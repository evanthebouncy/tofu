include("/home/evan/Documents/research/tofu/Distances.jl")
include("/home/evan/Documents/research/tofu/SimplePoly.jl")
include("/home/evan/Documents/research/tofu/Domain.jl")
using Optim

function get_max_value_exact(upper_bounds :: Function, dom_init, f_d_bounds = None)
  cnt = 0

  # negative of the maximum value as priority
  function get_priority(dom)
    -1.0 * upper_bounds(dom)
  end

  function interior(dom)
    for i in 1:length(dom)
      low, upp = dom[i]
      edge_low, edge_upp = dom_init[i]
      if (low == edge_low) | (upp == edge_upp)
        return false
      end
    end
    return true
  end

  # this constraint stats that if you are in the interior of the domain,
  # you must have the property that all your partial derivs can possibly be 0 for max to exist
  function constraint_partial(dom)
    if (interior(dom))
      all_bnds = [bnder(dom) for bnder in f_d_bounds]
      reduce(&, [bnd[1] - 1e-10 <=0.0<=bnd[2] + 1e-10 for bnd in all_bnds])
    else
      true
    end
  end

  doms = Collections.PriorityQueue()
  Collections.enqueue!(doms, dom_init, get_priority(dom_init))

  function sat_split_constraints!(dom)
    if dom in keys(doms)
      return false
    end
    if f_d_bounds == None
      true
    else
      constraint_partial(dom)
    end
  end

  while (cnt < 100)
    cnt += 1
    dom_cur = Collections.dequeue!(doms)
    @show(cnt, dom_cur, get_priority(dom_cur))
    half1, half2 = split_half(dom_cur)
    if sat_split_constraints!(half1)
      Collections.enqueue!(doms, half1, get_priority(half1))
    end
    if sat_split_constraints!(half2)
      Collections.enqueue!(doms, half2, get_priority(half2))
    end
  end
  dom_last = Collections.dequeue!(doms)
  upper_bounds(dom_last)
end

function get_max_value_approx(diff_func :: Function, dom, f_d = None)

  #provide the DifferentiableFunction object for the sampler
  function spp_enum(x::Vector)
    -1.0 * diff_func(x)
  end

  function spp_enum_grad!(x::Vector, storage::Vector)
    for i in 1:length(x)
      storage[i] = -1.0 * f_d[i](x)
    end
  end

  function spp_with_grad!(x::Vector, storage)
    for i in 1:length(x)
      storage[i] = -1.0 * f_d[i](x)
    end
    spp_enum(x)
  end

  ddf = (if f_d != None
    DifferentiableFunction(spp_enum, spp_enum_grad!, spp_with_grad!)
         else
           DifferentiableFunction(spp_enum)
         end
         )

  n_samples = length(dom) * 3

  function try_once()
    try
      xtol = max_length(dom) / 100
      ret = fminbox(ddf, get_single_sample(dom), Float64[x[1] for x in dom], Float64[x[2] for x in dom], xtol=xtol).f_minimum
      ret
    catch
      spp_enum(get_single_sample(dom))
    end
  end

  results = [try_once() for j in 1:n_samples]
  max([-1.0 * result for result in results]...)
end


# get the polynomial lower and upper bound over a domain for a potential f
# return a pair of SumPolyProdC object
function get_poly_lower_upper(f :: Potential, var_order, dom, degree, bnd_method="exact")
  spp_approx = get_m_projections_approx(f.potential_fun, var_order, degree, dom)

  function bound_diff(dom)
    Lf, Uf = f.potential_bnd(dom)
    Lp, Up = pbounds_phc(spp_approx, dom)
    Up-Lf, Uf - Lp
  end

  spp_enum(x) = peval(spp_approx, [x...])
  diff_fun_down(x) = spp_enum(x) - f.potential_fun(x...)
  diff_fun_up(x) = f.potential_fun(x...) - spp_enum(x)


  shift_down_bound(dom) = max(0, bound_diff(dom)[1])
  shift_up_bound(dom) = max(0, bound_diff(dom)[2])

  shift_down = (if bnd_method == "exact"
                  get_max_value_exact(shift_down_bound, dom)
                else
                  get_max_value_approx(diff_fun_down, dom)
                end
                )

  shift_up = (if bnd_method == "exact"
                get_max_value_exact(shift_up_bound, dom)
              else
                  get_max_value_approx(diff_fun_up, dom)
                end
              )

  SumPolyProdC(spp_approx, -1.0 * shift_down), SumPolyProdC(spp_approx, shift_up)
end
