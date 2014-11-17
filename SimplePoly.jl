using Iterators

# ============= POLYNOMIALS ===============

# Univariate Polynomials
# to be expressed in standard basis, the simplest for now
# later maybe beneficial to extend to fancier basis
immutable Poly1
  coef :: Array{Float64}
end

# equality
== (p1 :: Poly1, p2 :: Poly1) = (p1.coef == p2.coef)

# evaluating an uni-variate polynomial
function peval(p::Poly1, x)
  ret::Float64 = 0.0
  for i in 0:length(p.coef)-1
    ret = ret + p.coef[i+1] * x^i
  end
  ret
end

# adding uni-variate polynomials
function +(p1::Poly1, p2::Poly1)
  coef1 = p1.coef
  coef2 = p2.coef
  l = max(length(coef1), length(coef2))
  ext_coef1 = [coef1, [0 for i in 1:l-length(coef1)]]
  ext_coef2 = [coef2, [0 for i in 1:l-length(coef2)]]
  Poly1(ext_coef1 + ext_coef2)
end

# multipling uni-variate polynomials
function * (p1::Poly1, p2::Poly1)
  coef1, coef2 = p1.coef, p2.coef
  expanded_coefs = [[[0 for i in 1:idx], [c1*coef2[idx+1] for c1 in coef1]] for idx in 0:length(coef2)-1]
  if length(expanded_coefs) == 1
    return Poly1(expanded_coefs[1])
  end
  reduce(+, [Poly1(cs) for cs in expanded_coefs])
end

# integrating polynomials
function ∫(p1::Poly1)
  coef1 = p1.coef
  int_coef = [0, [coef1[i] / i for i in 1:length(coef1)]]
  Poly1(int_coef)
end

# chebyshev points between [-1, 1]
function get_chebyshev_pts(n)
  [cos(j*pi/n) for j in 0:n]
end

# extend [-1, 1] to the range of [a, b]
function get_range_denormalizer(a,b)
  range_denomalizer(x) = (b - a) / 2 * (x + 1) + a
  range_denomalizer
end

# do the uni_variate interpolation at tchebychev points and give coefficients in std basis
# a simple algorithm, does not use chebyshev polynomials, but for low degrees it is ok
function get_univar_interpolant_coef(xs::Array{Float64}, ys::Array{Float64})
  n = length(xs)
  mat_vec = vcat([map(x->x^i, xs) for i in 0:n-1]...)
  mat = reshape(mat_vec, n, n)
  inv(mat) * ys
end

# do the uni_variate approximation by interpolation at chebyshev points
function get_univar_approximation(f, n, a, b)
  denorm = get_range_denormalizer(a,b)
  xs = map(denorm, get_chebyshev_pts(n))
  ys = Float64[f(x) for x in xs]
  Poly1(get_univar_interpolant_coef(xs,ys))
end


# Multivariate Polynomials
# To be expressed in a summation of products
# Σ Σ Σ c_ijk * Pi(x)Qj(y)Rk(z)


# Product of univar polynomials
# const * P(x) * Q(y) * R(z)
immutable PolyProd
  c :: Float64
  var_order :: Array{ASCIIString}
  polys :: Dict{ASCIIString, Poly1}
end

# evaluation
function peval(pp :: PolyProd, args :: Array{Float64})
  assert(length(args) == length(pp.var_order))
  ret = pp.c
  for i in 1:length(pp.polys)
    ret *= peval(pp.polys[pp.var_order[i]], args[i])
  end
  ret
end

function peval(pp :: PolyProd, asmts :: Dict{ASCIIString, Float64})
  ret = pp.c
  for var_name in keys(asmts)
    ret *= peval(pp.polys[var_name], asmts[var_name])
  end
  ret
end

# takes 2 var orders and merge them to 1, perserves some order
function merge_varorder(vo1, vo2)
  [[x for x in vo1], filter(y->!(y in vo1), vo2)]
end

function can_join (pp1 :: PolyProd, pp2 :: PolyProd)
  pp1.var_order == pp2.var_order &&
    pp1.polys == pp2.polys
end

# can only add if the terms match
function + (pp1 :: PolyProd, pp2 :: PolyProd)
  assert(can_join(pp1, pp2))
  c = pp1.c + pp2.c
  PolyProd(c, pp1.var_order, pp1.polys)
end

function * (pp1 :: PolyProd, pp2 :: PolyProd)
  c = pp1.c * pp2.c
  var_order = merge_varorder(pp1.var_order, pp2.var_order)
  polys = Dict{ASCIIString, Poly1}()
  for var_name in var_order
    polys[var_name] = Poly1([1])
    if var_name in pp1.var_order
      polys[var_name] = polys[var_name] * pp1.polys[var_name]
    end
    if var_name in pp2.var_order
      polys[var_name] = polys[var_name] * pp2.polys[var_name]
    end
  end
  PolyProd(c, var_order, polys)
end

function ∫ (pp :: PolyProd, x, a, b)
  assert(x in pp.var_order)
  poly1 = pp.polys[x]
  poly1_int = ∫(poly1)
  c = pp.c * (peval(poly1_int,b) - peval(poly1_int, a))
  var_order = filter(name->(name != x), pp.var_order)
  polys = [name => pp.polys[name] for name in var_order]
  PolyProd(c, var_order, polys)
end


# Sums of Products Of (univar)Polynomials
immutable SumPolyProd
  var_order :: Array{ASCIIString}
  polyprods :: Array{PolyProd}
end

# try to successively add to a list of PolyProd while reducing like-terms
function add_polyprod!(pp1 :: PolyProd, pplist :: Array{PolyProd})
  for i in 1:length(pplist)
    pp2 = pplist[i]
    if can_join(pp1, pp2)
      splice!(pplist, i)
      push!(pplist, pp1 + pp2)
      return pplist
    end
  end
  push!(pplist, pp1)
end


# evaluation
function peval(spp :: SumPolyProd, args :: Array{Float64})
  assert(length(args) == length(spp.var_order))
  ret = 0
  for pp in spp.polyprods
    ret += peval(pp, args)
  end
  ret
end

function peval(spp :: SumPolyProd, asmts :: Dict{ASCIIString, Float64})
  ret = 0
  for pp in spp.polyprods
    ret += peval(pp, asmts)
  end
  ret
end


function * (spp1 :: SumPolyProd, spp2 :: SumPolyProd)
  var_order = merge_varorder(spp1.var_order, spp2.var_order)
  polyprods = PolyProd[]
  for pp1 in spp1.polyprods
    for pp2 in spp2.polyprods
      push!(polyprods, pp1 * pp2)
    end
  end
  SumPolyProd(var_order, polyprods)
end

function ∫ (spp :: SumPolyProd, x :: ASCIIString, a, b)
  assert(x in spp.var_order)
  var_order = filter(name->(name != x), spp.var_order)
  polyprods = PolyProd[]
  for polyprod in spp.polyprods
    to_add = ∫(polyprod, x, a, b)
    add_polyprod!(to_add, polyprods)
  end
  SumPolyProd(var_order, polyprods)
end

# Do the approximation into a multi-variate polynomial as a SumPolyProd

# get a d dimensional tchebyshev grid of n sample in each dimension
# over the domain of the box
function get_chebyshev_grid(dim, n, box_domain)
  scaled_nodes = [map(get_range_denormalizer(box_domain[i]...), get_chebyshev_pts(n)) for i in 1:dim]
  product(scaled_nodes...)
end

# get the greatest error of a function over a grid
function sample_greatest_err(f, grid)
  max([(f(pt...),pt) for pt in grid]...)
end

# get the function on each lines radiating from a point evaluated
function get_univar_func(f, pt)
  ret = Function[]
  for idx in 1:length(pt)
    push!(ret, (x) -> (
            eval_pt = map(j->(j == idx ? x : pt[j]), 1:length(pt));
            f(eval_pt...)
          ))
  end
  return ret
end

# get the approximation by pivot at the greatest error points
# f the function
# var_order the list of variable names
# box_domain the box of region to do the approx over
# number of sample points (degree - 1)
# give a PolyProd output
function get_single_approx(f, var_order, n, box_domain)
  dim = length(var_order)
  # set up the grid to do the sampling
  grid = get_chebyshev_grid(dim, n, box_domain)
  # find by sampling the biggest error
  err, pivot = sample_greatest_err(f, grid)
  # construct univariate functions around the pivot
  univar_funs = get_univar_func(f, pivot)
  c = 1 / err ^ (dim - 1)
  univar_funs_approxs = Dict{ASCIIString, Poly1}()
  for i in 1:dim
    dom_lower, dom_upper = box_domain[i];
    univar_approx = get_univar_approximation(univar_funs[i], n, dom_lower, dom_upper)
    univar_funs_approxs[var_order[i]] = univar_approx
  end
  PolyProd(c, var_order, univar_funs_approxs)
end

# get (rep) number of projections onto the outer product space
# approximate a function as a SumPolyProd
function get_m_projections_approx(func, var_order, n, box_domain, rep=None)
  if rep == None
    rep = n
  end

  function rec_get_projection(err_fun :: Function, pp_list :: Array{PolyProd}, proj_togo :: Int64)
    if proj_togo == 0
      pp_list
    else
      pp_prox = get_single_approx(err_fun, var_order, n, box_domain)
      nxt_pp_list = [pp_list, pp_prox]
      nxt_err_fun = (x...) -> err_fun(x...) - peval(pp_prox, [x...])
      rec_get_projection(nxt_err_fun, nxt_pp_list, proj_togo - 1)
    end
  end
  SumPolyProd(var_order, rec_get_projection(func, PolyProd[], rep))
end