include("/home/evan/Documents/research/tofu/Distances.jl")
include("/home/evan/Documents/research/tofu/SimplePoly.jl")
include("/home/evan/Documents/research/tofu/Domain.jl")
using Optim

function get_max_value(f_bounds :: Function, dom_init, f_d_bounds = None)
  cnt = 0

  # negative of the maximum value as priority
  function get_priority(dom)
    bound = f_bounds(dom)
    -1.0 * bound[2]
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

  while (cnt < 500)
    cnt += 1
    @show(cnt, length(doms))
    dom_cur = Collections.dequeue!(doms)
    @show(max_length(dom_cur), [f_d_bound(dom_cur) for f_d_bound in f_d_bounds], f_bounds(dom_cur))
    half1, half2 = split_half(dom_cur)
    if sat_split_constraints!(half1)
      @show("queued", half1)
      Collections.enqueue!(doms, half1, get_priority(half1))
    end
    if sat_split_constraints!(half2)
      @show("queued", half2)
      Collections.enqueue!(doms, half2, get_priority(half2))
    end
  end
  dom_last = Collections.dequeue!(doms)
  f_bounds(dom_last), doms
end

# get the polynomial abstraction over a domain for a potential f
function get_poly_abstraction(f :: Potential, var_order, dom, degree)
  spp_approx = get_m_projections_approx(f.potential_fun, var_order, degree, dom)

  function positive_err_bound(dom)
  end

  shift_down, shift_up = get_maximum_errors(f, spp_approx, dom)
  SumPolyProdC(spp_approx, -1.0 * shift_down), SumPolyProdC(spp_approx, shift_up)
end

# multiplication is tricky because of the negative numbers
# we need to explicitly force it to be >= 0
immutable SPPCMult
  sppcs :: Array{SumPolyProdC}
end

function peval(sppcm :: SPPCMult, x)
  reduce(*, [max(0, peval(sppc, x)) for sppc in sppcm.sppcs])
end
