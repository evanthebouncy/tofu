using Iterators
using Optim

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
  ret = 0.0
  for i in 0:length(p.coef)-1
    ret = ret + p.coef[i+1] * x^i
  end
  ret
end

# The Durand_Kerner method for finding root of a uni-variate polynomial in complex
# The leading coefficient for p (the one with the highest degree) must be 1
function Durand_Kerner(p, deg)
  alph = 0.9 + 0.4im
  vect = [alph^i for i in 1:deg]
  for i in 1:30
    vect_new = []
    for j in 1:deg
      rootj = vect[j]
      others = [vect[i] for i in filter(x->x!=j, 1:deg)]
      rootj_new = rootj - p(rootj) / reduce(*, Complex64[rootj - rooti for rooti in others])
      vect_new = [vect_new, rootj_new]
    end
    vect = vect_new
  end
  vect
end

# find the real roots of a uni-variate polynomial
function real_roots(p::Poly1)
  coefs = p.coef
  if length(coefs) <= 1
    return Float64[]
  end
  last_term = last(coefs)
  normalized_coef = coefs / last_term
  normalized_poly = Poly1(normalized_coef)
  normalized_poly_numer(x) = peval(normalized_poly, x)
  roots = Durand_Kerner(normalized_poly_numer, length(coefs) - 1)
  Float64 [real(root) for root in filter(x-> imag(x) < 0.00001, roots)]
end

# find the minimum and maximum of a Poly1 over a domain
# assume roots of the derivitive are given (they maybe pre-computed)
function pbounds(p::Poly1, roots_of_deriv, dom)
  included_roots = filter(x->(dom[1]<x<dom[2]), roots_of_deriv)
  crit_pts = [dom[1], dom[2], included_roots]
  values = [peval(p, pt) for pt in crit_pts]
  min(values...), max(values...)
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

# negation
function -(p1::Poly1)
  coef = [-1.0 * x for x in p1.coef]
  Poly1(coef)
end

# subtraction
function -(p1::Poly1, p2::Poly1)
  p1 + (- p2)
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

# take derivitive of polynomial
function δ(p1::Poly1)
  coef1 = p1.coef
  derv_coef = [coef1[i] * (i-1) for i in 2:length(coef1)]
  Poly1(derv_coef)
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

# find bounds of a PolyProd object
function pbounds(pp :: PolyProd, roots_dict, rng)
  minmaxs = [pbounds(pp.polys[pp.var_order[i]],
                     roots_dict[pp.polys[pp.var_order[i]]],
                     rng[i])
             for i in 1:length(pp.polys)]
  cross = product(minmaxs...)
  values = pp.c * [reduce(*, tup) for tup in cross]
  min(values...), max(values...)
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

# negate
function - (pp1 :: PolyProd)
  c = -pp1.c
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

function δ (pp :: PolyProd, x)
  assert(x in pp.var_order)
  poly1 = pp.polys[x]
  poly1_der = δ(poly1)
  c = pp.c
  var_order = copy(pp.var_order)
  polys = copy(pp.polys)
  polys[x] = poly1_der
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

# give a spp, the root of its derivitives, and a domain
# provide a bound
function pbounds(spp::SumPolyProd, root_dict, dom)
  values = [pbounds(pp, root_dict, dom) for pp in spp.polyprods]
  minval = reduce(+, map(x->x[1], values))
  maxval = reduce(+, map(x->x[2], values))
  minval, maxval
end

# given a spp, create a bound object that, when taking in a domain
# give a lower / upper bound of it
function get_bound_object(spp::SumPolyProd)
  root_dict = Dict{Poly1, Array{Float64}}()
  for pp in spp.polyprods
    for v_name in pp.var_order
      p = pp.polys[v_name]
      if !(p in keys(root_dict))
        root_dict[p] = real_roots(δ(p))
      end
    end
  end

  function bound_object(dom)
    pbounds(spp, root_dict, dom)
  end

  bound_object
end

# arithmentic operations

# addition
function + (spp1 :: SumPolyProd, spp2 :: SumPolyProd)
  var_order = merge_varorder(spp1.var_order, spp2.var_order)
  polyprods = [spp1.polyprods, spp2.polyprods]
  SumPolyProd(var_order, polyprods)
end

# negation
function - (spp1 :: SumPolyProd)
  var_order = spp1.var_order
  polyprods = [PolyProd(-1.0 * pp.c, pp.var_order, pp.polys) for pp in spp1.polyprods]
  SumPolyProd(var_order, polyprods)
end

# subtraction
function - (spp1 :: SumPolyProd, spp2 :: SumPolyProd)
  spp1 + (- spp2)
end

# multiplication
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

function δ (spp :: SumPolyProd, der_var :: ASCIIString)
  var_order = copy(spp.var_order)
  polyprods = PolyProd[δ(pp, der_var) for pp in spp.polyprods]
  SumPolyProd(var_order, polyprods)
end

function get_single_sample(dom1)
  left_end = [x[1] for x in dom1]
  lengthz  = [x[2]-x[1] for x in dom1]
  [left_end[i] + rand() * lengthz[i] for i in 1:length(dom1)]
end

# use n random starting points combined with gradient descent to
# find a point of maximum within the domain
function grad_approx_max(ddf :: DifferentiableFunction, dom, n_samples)
  results = [fminbox(ddf, get_single_sample(dom), [x[1] for x in dom], [x[2] for x in dom])]
  max([-1.0 * result.f_minimum for result in results]..., -Inf)
end

# for a spp solve for the maximum of it
function solve_max(spp :: SumPolyProd, dom_init)
  # some constants and a counter
  dom_tol = 1.0e-10
  range_tol = 1.0e-10
  # for information only
  cnt = 0
  best_dom = (0.0, 0.0)
  best_bnd = (0.0, 0.0)

  # divide a domain in half along the greatest axis
  function split_half(dom)
    dom_w_length = [(d[2]-d[1],d) for d in dom]
    size, maxd = max(dom_w_length...)
    half1 = (maxd[1], (maxd[1] + maxd[2]) / 2)
    half2 = ((maxd[1] + maxd[2]) / 2, maxd[2])
    ret1 = [x for x in dom]
    ret1[findfirst(ret1, maxd)] = half1
    ret2 = [x for x in dom]
    ret2[findfirst(ret2, maxd)] = half2
    ret1, ret2
  end

  function max_length(dom)
    max([d[2]-d[1] for d in dom]...)
  end



  function interior(dom)
    return true
    for i in 1:length(dom)
      low, upp = dom[i]
      edge_low, edge_upp = dom_init[i]
      if (low == edge_low) | (upp == edge_upp)
        return false
      end
    end
    return true
  end

  # give the bounding object for the function itself
  value_bound_obj = get_bound_object(spp)

  # give the bounding objects for all the partial derivatives
  all_partials = SumPolyProd[δ(spp, x) for x in spp.var_order]
  all_partials_bnd_obj = [get_bound_object(spp1) for spp1 in all_partials]

  #provide the DifferentiableFunction object for the sampler
  function spp_enum(x::Vector)
    -1.0 * peval(spp, x)
  end
  function spp_enum_grad!(x::Vector, storage::Vector)
    for i in 1:length(x)
      storage[i] = -1.0 * peval(all_partials, x)
    end
  end
  function spp_with_grad!(x::Vector, storage)
    for i in 1:length(x)
      storage[i] = -1.0 * peval(all_partials[i], x)
    end
    spp_enum(x)
  end
  ddf = DifferentiableFunction(spp_enum, spp_enum_grad!, spp_with_grad!)


  # this constraint stats that if you are in the interior of the domain,
  # you must have the property that all your partial derivs can possibly be 0 for max to exist
  function constraint_partial(dom)
    if (interior(dom))
      all_bnds = [bnder(dom) for bnder in all_partials_bnd_obj]
      reduce(&, [bnd[1]<=0.0<=bnd[2] for bnd in all_bnds])
    else
      true
    end
  end

  # this constraint stats that if you want to have a shot of being the max
  # you better have an interval whose maximum is greater than the current max
  function constraint_max(dom)
    low, upp = value_bound_obj(dom)
    (upp - cur_max) > range_tol
  end

  # this constraint stats that if your dom_tol is too small
  # we should stop splitting you
  function constraint_dom_tol!(dom)
    if (max_length(dom) < dom_tol)
      println("giving up, dom too small, updating cur_max with whats provable on this region")
      cur_max = max(cur_max, value_bound_obj(dom)[2])
    end
    max_length(dom) > dom_tol
  end

  # if you satisify this you will be split and put back into the queue
  function sat_split_constraints(dom)
    constraint_partial(dom) &
      constraint_max(dom) &
      constraint_dom_tol!(dom)
  end

  dom_tol = max_length(dom_init) / 100
  doms = typeof(dom_init)[]
  push!(doms, dom_init)
  # the current maximum we want to keep track of, initialized at a random point's value
  # cur_max = peval(spp, get_single_sample(dom_init))
  cur_max = -Inf
  while ( length(doms) > 0 )
    cnt += 1
    @show(cnt, length(doms), cur_max)
    dom_cur = pop!(doms)

    # attempt to improve max by sampling, may fail with linsearch converge
    try
      cur_max = max(cur_max, grad_approx_max(ddf, dom_init, 3))
      range_tol = 0.01 * cur_max
    catch

    end

    low, upp = value_bound_obj(dom_cur)
    if upp - low < range_tol
      if (upp > cur_max)
        best_dom = dom_cur
        best_bnd = (low, upp)
      end
      cur_max = max(cur_max, upp)
    end
    if sat_split_constraints(dom_cur)
      half1, half2 = split_half(dom_cur)
      push!(doms, half1)
      push!(doms, half2)
    end
  end
  cur_max + range_tol
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

# ============ The STD basis for polynomials =======================
# ============ Useful for finding differences of polynomials =======

type STDPolyTerm
  c :: Float64
  var_order :: Array{ASCIIString}
  term :: Dict{ASCIIString, Int64}
end

function to_print(term :: STDPolyTerm)
  print(term.c)
  for var in term.var_order
    print(var)
    print("^")
    print(term.term[var])
  end
end

function can_join(spt1 :: STDPolyTerm, spt2 :: STDPolyTerm)
  spt1.term == spt2.term
end

function + (term1 :: STDPolyTerm, term2 :: STDPolyTerm)
  assert(can_join(term1, term2))
  STDPolyTerm(term1.c + term2.c, copy(term1.var_order), copy(term1.term))
end

function - (term1 :: STDPolyTerm)
  STDPolyTerm(-1.0 * term1.c, copy(term1.var_order), copy(term1.term))
end

function * (term1 :: STDPolyTerm, term2 :: STDPolyTerm)
  c = term1.c * term2.c
  var_order = merge_varorder(term1.var_order, term2.var_order)
  term = copy(term1.term)
  for var_name in keys(term2.term)
    if var_name in keys(term)
      term[var_name] = term[var_name] + term2.term[var_name]
    else
      term[var_name] = term2.term[var_name]
    end
  end
  STDPolyTerm(c, var_order, term)
end

function peval(term::STDPolyTerm, lst_vals::Array{Float64})
  assert(length(lst_vals) == length(term.var_order))
  ret = term.c
  for i in 1:length(lst_vals)
    var_name = term.var_order[i]
    var_exp = term.term[var_name]
    var_value = lst_vals[i]
    ret *= var_value ^ var_exp
  end
  ret
end

function peval(term::STDPolyTerm, lst_asmt::Dict{ASCIIString,Float64})
  # make sure all var_name are in the keys
  for var_name in term.var_order
    assert(var_name in keys(lst_asmt))
  end
  ret = term.c
  for var_name in term.var_order
    var_exp = term.term[var_name]
    var_value = lst_asmt[var_name]
    ret *= var_value ^ var_exp
  end
  ret
end

type STDPoly
  var_order :: Array{ASCIIString}
  terms :: Array{STDPolyTerm}
end

function to_print(pol :: STDPoly)
  for t in pol.terms
    to_print(t)
    print(" + ")
  end
end

function peval(stdpoly :: STDPoly, lst_vals :: Array{Float64})
  assert(length(stdpoly.var_order) == length(lst_vals))
  asmts = [stdpoly.var_order[i] => lst_vals[i] for i in 1:length(lst_vals)] :: Dict{ASCIIString,Float64}
  reduce(+, [peval(term, asmts) for term in stdpoly.terms])
end

function add1!(terms, other_term)
  for term in terms
    if can_join(term, other_term)
      c = term.c + other_term.c
      tterm = copy(term.term)
      splice!(terms, findfirst(terms, term))
      push!(terms, STDPolyTerm(c, term.var_order, tterm))
      return terms
    end
  end
  push!(terms, other_term)
end

function + (p1 :: STDPoly, p2 :: STDPoly)
  var_order = merge_varorder(p1.var_order, p2.var_order)
  terms = STDPolyTerm[term for term in p1.terms]

  for other_term in p2.terms
    add1!(terms, other_term)
  end
  STDPoly(var_order, terms)
end

function - (p1 :: STDPoly)
  var_order = copy(p1.var_order)
  terms = [-(term) for term in p1.terms]
  STDPoly(var_order, terms)
end

function - (p1 :: STDPoly, p2 :: STDPoly)
  p1 + (- p2)
end

function * (p1 :: STDPoly, p2 :: STDPoly)
  var_order = merge_varorder(p1.var_order, p2.var_order)
  all_prods = STDPolyTerm[t1 * t2 for (t1, t2) in product(p1.terms, p2.terms)]
  terms = STDPolyTerm[]
  for t in all_prods
    add1!(terms, t)
  end
  STDPoly(var_order, terms)
end

# convert SumPolyProd to STDPoly
# one at a time

function PolyProd_to_STDPoly(p1 :: PolyProd)
  var_order = p1.var_order
  c = p1.c
  stdpolys = STDPoly[ ]
  to_push =  STDPoly([], [STDPolyTerm(c, [], Dict{ASCIIString,Int64}())]) :: STDPoly
  push!(stdpolys,  to_push)
  for var_name in var_order
    coefs = p1.polys[var_name].coef
    terms = STDPolyTerm[]
    for i in 1:length(coefs)
      termc = coefs[i]
      term_term = Dict{ASCIIString, Int64}()
      term_var_order = ASCIIString[]
      if i > 1
        expon = i - 1
        push!(term_var_order, var_name)
        term_term[var_name] = expon
      end
      push!(terms, STDPolyTerm(termc, term_var_order, term_term))
    end
    push!(stdpolys, STDPoly([var_name], terms))
  end
  reduce(*, stdpolys)
end

function SumPolyProd_to_STDPoly(p1 :: SumPolyProd)
  std_polys = [PolyProd_to_STDPoly(p) for p in p1.polyprods]
  reduce(+, std_polys)
end

