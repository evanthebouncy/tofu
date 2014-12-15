include("/home/evan/Documents/research/tofu/Distances.jl")
include("/home/evan/Documents/research/tofu/SimplePoly.jl")
include("/home/evan/Documents/research/tofu/Domain.jl")
using Optim

# get maximum error
# given a Potential f and a polynomial p over a domain d,
# give the bounds btwn them, i.e. a and b such that
# p - a <= f <= p + b. i.e. we have a tube containing f
# i.e. a >= p - f  and  b >= f - p
function get_maximum_errors(f :: Potential, p :: SumPolyProd, dom_init)
  # book keepings
  cnt = 0
  # some constants
  dom_tol = max_length(dom_init) / 1000
  range_tol = 1.0e-10

  cur_a = -Inf
  cur_b = -Inf

  doms = typeof(dom_init)[]
  push!(doms, dom_init)

  function numer_a(x :: Vector)
    peval(p, x) - f.potential_fun(x...)
  end
  ddf_a = DifferentiableFunction(numer_a)
  function numer_b(x :: Vector)
    f.potential_fun(x...) - peval(p, x)
  end
  ddf_b = DifferentiableFunction(numer_b)

  p_bnd_obj = get_bound_object(p)

  function value_bound_ab(dom)
    p_low, p_high = p_bnd_obj(dom)
    f_low, f_high = f.potential_bnd(dom)
    a_bnd = p_high - f_low
    b_bnd = f_high - p_low
    a_bnd, b_bnd
  end

  function sat_split_constraints(dom)
    a_bnd, b_bnd = value_bound_ab(dom)
    (a_bnd > cur_a) | (b_bnd > cur_b)
  end

  while ( length(doms) > 0 )
    cnt += 1
    @show(cnt, length(doms), cur_a, cur_b)
    dom_cur = pop!(doms)

    # attempt to improve cur_a and cur_b, may fail with linsearch converge
    try
      cur_a = max(cur_a, grad_approx_max(ddf_a, doms[rand(1:length(doms))], 3))
      cur_b = max(cur_b, grad_approx_max(ddf_b, doms[rand(1:length(doms))], 3))
      range_tol = 0.01 * 0.5 * (cur_a + cur_b)
    catch
    end

    a_bnd, b_bnd = value_bound_ab(dom_cur)
    if a_bnd + b_bnd < range_tol
      cur_a = max(cur_a, a_bnd)
      cur_b = max(cur_b, b_bnd)
    end
    if sat_split_constraints(dom_cur)
      half1, half2 = split_half(dom_cur)
      push!(doms, half1)
      push!(doms, half2)
    end
  end
  cur_a + range_tol, cur_b + range_tol
end

# get the polynomial abstraction over a domain for a potential f
function get_poly_abstraction(f :: Potential, var_order, dom, degree)
  spp_approx = get_m_projections_approx(f.potential_fun, var_order, degree, dom)
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
