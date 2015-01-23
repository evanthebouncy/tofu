include("/home/evan/Documents/research/tofu/Abstraction.jl")

# a patch is a region along with a polynomial lower and upper bound
type Patch
  dom :: Domain
  p_l :: SumPolyProdC
  p_u :: SumPolyProdC
end

# a factor consist of a variable ordering, a binary_space_tree
# which has leaves that forms a partition. It also associate each domain in the bsp
# with a Patch, dictating which lower/upper bound it must provide over that abstraction
type Factor
  f_name
  var_order
  bsp :: BSP
  partition :: Partition
  potential_bounds :: Dict{Domain, Patch}
end

# this keeps relationships between factors as operations Mult, and Inte
abstract FactorOp

type FactorMult <: FactorOp
  f1 :: Factor
  f2 :: Factor
  f_mult :: Factor
end

type FactorInte <: FactorOp
  f1 :: Factor
  inte_var :: ASCIIString
  f_inte :: Factor
end

# this is analegous relationship
abstract DomainOp

type DomainMult <: DomainOp
  fmultop :: FactorMult
  d1 :: Domain
  d2 :: Domain
  d_mult :: Domain
end

type DomainInte <: DomainOp
  finteop :: FactorInte
  d1s :: Array{Domain}
  d_inte :: Domain
end


# A collection of Factors and their relationships (and the domain relationships)
type FactorGraph
  factors :: Array{Factor}
  # a mapping between a factor and a FactorOp remembering its parents
  # e.g. f3 -> FactorMult(f1, f2, f3), i.e. f1 and f2 are parents of f3
  rel_factor_parents :: Dict{Factor, FactorOp}
  # a mapping between a factor and a FactorOp remember its child
  # e.g. f1 -> FactorMult(f1, f2, f3), i.e. f1 has child f3 via a mult with f2
  rel_factor_child :: Dict{Factor, FactorOp}
  # a mapping between the pair (Factor,Domain) to the parents whose contribution made
  # this particular domain in this factor.
  # e.g. (f3, dom3) -> DomainMult(FactorMult(f1,f2,f3), d1, d2, d3), i.e. (f3, dom3) has parent
  # d1 from f1, and d2 from f2. Each (fact,dom) can have only one, unique, domainOp associated with it
  rel_domain_parents :: Dict{(Factor, Domain), DomainOp}
  # similar as above, except remember its multiple childrens. it is possible for a single (Factor,Domain)
  # to register multiple DomainOp due to overlapping
  rel_domain_children :: Dict{(Factor, Domain), Array{DomainOp}}
end

# ================================ methods to modify the factor graph and its relations ==================================

# initialize a factor graph
function init_factor_graph()
  FactorGraph(Factor[], Dict{Factor,FactorOp}(), Dict{Factor,FactorOp}(),
              Dict{(Factor,Domain), DomainOp}(), Dict{(Factor,Domain),Array{DomainOp}}())
end

# adding a factor to the list of factors in the FG
function add_factor!(FG :: FactorGraph, f :: Factor)
  push!(FG.factors, f)
end

# add (factor, domain) -> domainOp to the FactorGraph.rel_domain_children
function add_rel_domain_children!(FG :: FactorGraph, f :: Factor, d :: Domain, dop :: DomainOp)
  key = (f, d)
  if !(key in keys(FG.rel_domain_children))
    FG.rel_domain_children[key] = DomainOp[]
  end
  push!(FG.rel_domain_children[key], dop)
end

# add a relation domainOP to the factor graph
function add_rel_domain!(FG :: FactorGraph, d_op :: DomainOp)
  if typeof(d_op) == DomainMult
    FG.rel_domain_parents[(d_op.fmultop.f_mult, d_op.d_mult)] = d_op
    add_rel_domain_children!(FG, d_op.fmultop.f1, d_op.d1, d_op)
    add_rel_domain_children!(FG, d_op.fmultop.f2, d_op.d2, d_op)
    return
  end
  if typeof(d_op) == DomainInte
    FG.rel_domain_parents[(d_op.finteop.f_inte, d_op.d_inte)] = d_op
    for p_domain in d_op.d1s
      add_rel_domain_children!(FG, d_op.finteop.f1, p_domain, d_op)
    end
    return
  end
  assert(false)
end

# add a relation FactorOp to the factorGraph
function add_rel_factor!(FG :: FactorGraph, f_op :: FactorOp)
  if typeof(f_op) == FactorMult
    FG.rel_factor_parents[f_op.f_mult] = f_op
    FG.rel_factor_child[f_op.f1] = f_op
    FG.rel_factor_child[f_op.f2] = f_op
    return
  end
  if typeof(f_op) == FactorInte
    FG.rel_factor_parents[f_op.f_inte] = f_op
    FG.rel_factor_child[f_op.f1] = f_op
    return
  end
  assert(false)
end


function feval_lower(f :: Factor, x)
  contained_dom = find_best_containing_domain(f.bsp, x)
  peval(f.potential_bounds[contained_dom].p_l, x)
end

function feval_upper(f :: Factor, x)
  contained_dom = find_best_containing_domain(f.bsp, x)
  peval(f.potential_bounds[contained_dom].p_u, x)
end

# shatter the domain into 2^k pieces where k is the dimension
# return the resulting domains
function shatter!(var_order, dom :: Domain, p :: Partition, bsp :: BSP)

  function split_once!(var_order, dom :: Domain, p :: Partition, bsp :: BSP)
    sp, d1, d2 = split_half(dom)
    setdiff!(p, Domain[dom])
    union!(p, Domain[d1, d2])
    split_var = var_order[sp]
    split_val = 0.5 * (dom[sp][1] + dom[sp][2])
    grow_bsp!(bsp, dom, d1, d2, split_var, split_val)
    Domain[d1, d2]
  end

  function shatter!_rec(var_order, doms :: Array{Domain}, p :: Partition, bsp :: BSP, n)
    if n == 0
      doms
    else
      new_doms = Domain[]
      for d in doms
        new_doms = [split_once!(var_order, d, p, bsp), new_doms]
      end
      shatter!_rec(var_order, new_doms, p, bsp, n - 1)
    end
  end

  shatter!_rec(var_order, Domain[dom], p, bsp, length(var_order))
end

# make a patch from a potential function
function patch_pot(FG :: FactorGraph, var_order, dom :: Domain, pot :: Potential, bnd_method, deg)
  lower_poly, upper_poly = get_poly_lower_upper(pot, var_order, dom, deg, bnd_method)
  lowest_of_lower = lowest_possible(lower_poly, dom, bnd_method)
  lower_poly = (if lowest_of_lower < -1e-6
                  make_zero(lower_poly)
                else
                  lower_poly
                end
                )
  Patch(dom, lower_poly, upper_poly)
end

# make a factor from a potential function
function f_pot(FG :: FactorGraph, f_name, pot :: Potential, var_order, dom :: Domain, bnd_method="exact", deg=2)
  bsp = new_leaf(dom)
  partition = Set{Domain}(Domain[dom])
  potential_bounds = Dict{Domain, Patch}()
  potential_bounds[dom] = patch_pot(FG, var_order, dom, pot, bnd_method, deg)
  ret = Factor(f_name, var_order, bsp, partition, potential_bounds)
  add_factor!(FG, ret)
  ret
end

# grow a factor from a potential function
function f_pot_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain, pot :: Potential, bnd_method="exact", deg=2)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  # add the patches for the shattered domains
  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_pot(FG, f.var_order, d1, pot, bnd_method, deg)
  end
end

function rand_f_pot_grow!(FG :: FactorGraph, f :: Factor, pot :: Potential, bnd_method="exact", deg=2)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_pot_grow!(FG, f, dom_to_grow, pot, bnd_method, deg)
end

# give a domain, and 2 factors f1 and f2
function patch_mult (FG :: FactorGraph, var_order, dom :: Domain, f1 :: Factor, f2 :: Factor, f_mult :: Factor)
  shrink1 = diminish_dom_dim(var_order, f1.var_order, dom)
  shrink2 = diminish_dom_dim(var_order, f2.var_order, dom)
  best_cover1 = best_covering(shrink1, keys(f1.potential_bounds))
  best_cover2 = best_covering(shrink2, keys(f2.potential_bounds))
  patch1, patch2 = f1.potential_bounds[best_cover1], f2.potential_bounds[best_cover2]
  p_l = patch1.p_l * patch2.p_l
  p_u = patch1.p_u * patch2.p_u
  # register relationship in FG
  fact_mult = FactorMult(f1, f2, f_mult)
  domain_mult = DomainMult(fact_mult, best_cover1, best_cover2, dom)
  add_rel_domain!(FG, domain_mult)
  Patch(dom, p_l, p_u)
end

function f_mult(FG :: FactorGraph, f_name, f1 :: Factor, f2 :: Factor, cover_dom :: Domain)
  var_order = merge_varorder(f1.var_order, f2.var_order)
  bsp = new_leaf(cover_dom)
  partition = Set{Domain}(Domain[cover_dom])
  potential_bounds = Dict{Domain, Patch}()
  ret = Factor(f_name, var_order, bsp, partition, potential_bounds)
  # register fact relationship to FG
  push!(FG.factors, ret)
  add_rel_factor!(FG, FactorMult(f1, f2, ret))
  ret.potential_bounds[cover_dom] = patch_mult(FG, var_order, cover_dom, f1, f2, ret)
  ret
end

# grow a factor from a multiplication
function f_mult_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain, f1 :: Factor, f2 :: Factor)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_mult(FG, f.var_order, d1, f1, f2, f)
  end
end

function rand_f_mult_grow!(FG :: FactorGraph, f :: Factor, f1 :: Factor, f2 :: Factor)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_mult_grow!(FG, f, dom_to_grow, f1, f2)
end


# give a domain, and the factor it integrates from
function patch_integrate (FG :: FactorGraph, var_order, dom :: Domain, f1 :: Factor, inte_var_name, f_inte :: Factor)
  smallest_cover = find_smallest_cover(f1.bsp, dom, var_order, f1.var_order, inte_var_name)
  integrated_pots = (SumPolyProdC,SumPolyProdC)[]
  for s_dom in smallest_cover
    inte_patch = f1.potential_bounds[s_dom]
    lower, upper = inte_patch.p_l, inte_patch.p_u
    bound_left, bound_right = s_dom[findfirst(f1.var_order, inte_var_name)]
    lower_poly_bnd = ∫(lower, inte_var_name, bound_left, bound_right)
    upper_poly_bnd = ∫(upper, inte_var_name, bound_left, bound_right)
    push!(integrated_pots, (lower_poly_bnd, upper_poly_bnd))
  end
  lower_pol, upper_pol = (if length(integrated_pots) == 1
                            integrated_pots[1]
                          else
                            reduce((x,y)-> (x[1]+y[1],x[2]+y[2]), integrated_pots)
                          end
                          )
  # register relationship in FG
  fact_inte = FactorInte(f1, inte_var_name, f_inte)
  domain_inte = DomainInte(fact_inte, [s_dom for s_dom in smallest_cover], dom)
  add_rel_domain!(FG, domain_inte)
  Patch(dom, lower_pol, upper_pol)
end

function f_inte(FG :: FactorGraph, f_name, f1 :: Factor, inte_var_name :: ASCIIString, cover_dom :: Domain)
  var_order = filter(x->x!=inte_var_name, f1.var_order)
  bsp = new_leaf(cover_dom)
  partition = Set{Domain}(Domain[cover_dom])
  potential_bounds = Dict{Domain, Patch}()
  ret = Factor(f_name, var_order, bsp, partition, potential_bounds)
  ret.potential_bounds[cover_dom] = patch_integrate(FG, var_order, cover_dom, f1, inte_var_name, ret)
  # register fact relationship to FG
  push!(FG.factors, ret)
  add_rel_factor!(FG, FactorInte(f1, inte_var_name, ret))
  ret
end

# grow a factor from a potential function
function f_inte_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain, f1 :: Factor, inte_var_name :: ASCIIString)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_integrate(FG, f.var_order, d1, f1, inte_var_name, f)
  end
end

function rand_f_inte_grow!(FG :: FactorGraph, f :: Factor, f1 :: Factor, inte_var_name :: ASCIIString)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_inte_grow!(FG, f, dom_to_grow, f1, inte_var_name)
end

# measurement of errors ===========================

# given a factor, and a set of domains in that factor
# if this factor is a multiplication relation, return a list of length 2 of [[(f1, dom_f1)...], [(f2, dom_f2)...]]
# if this factor is an integration relation, return a list of length 1 of [[f1, dom_f1)...]
# the domain whom contributed to the list of domains in the argument
function one_step_upstream_contribution(factor :: Factor, doms :: Set{Domain})
  if !(factor in keys(FG.rel_factor_parents))
    return
  end
  # grab the f_type of this factor
  f_type = FG.rel_factor_parents[factor]
  # if it is multiplication...
  if typeof(f_type) == FactorMult
    f1 = f_type.f1
    f2 = f_type.f2
    dom1s = Set{Domain}()
    dom2s = Set{Domain}()
    for dom in doms
      rel_multi = FG.rel_domain_parents[(factor, dom)]
      union!(dom1s, Domain[rel_multi.d1])
      union!(dom2s, Domain[rel_multi.d2])
    end
    return [(f1, dom1s), (f2, dom2s)]
  end
  # if it is integration...
  if typeof(f_type) == FactorInte
    f1 = f_type.f1
    dom1s = Set{Domain}()
    for dom in doms
      rel_inte = FG.rel_domain_parents[(factor, dom)]
      union!(dom1s, rel_inte.d1s)
    end
    return [(f1, dom1s)]
  end
end

function all_upstream_contributions_rec(factor :: Factor, doms :: Set{Domain})
  attempt_upstream = one_step_upstream_contribution(factor, doms)
  if (attempt_upstream == Nothing())
    [(factor, doms)]
  else
    ret = (Factor, Array{Domain})[]
    for factor_doms in attempt_upstream
      ret = [ret, all_upstream_contributions_rec(factor_doms...)]
    end
    push!(ret, (factor, doms))
    ret
  end
end

# find out all the domains that could be further shattered with effect
function get_valid_split_domains(FG :: FactorGraph)
  result_factor = last(FG.factors)
  @show([f.f_name for f in FG.factors])
  upstreams = all_upstream_contributions_rec(result_factor, result_factor.partition)
  @show(length(upstreams))
  ret = (Factor, Domain)[]
  for fact_dom in upstreams
    fact, doms = fact_dom
    filtered_doms = filter(d->d in fact.partition, doms)
    @show(filtered_doms)
    for filt_dom in filtered_doms
      push!(ret, (fact, filt_dom))
    end
  end
  ret
end

# pretend the other bounds are prefect i.e. lower == upper
# we propagate the imperfect lower bound
function propagate_bad_lower_bnd_mult(f_mult :: Factor, f_src_bad :: Factor, bad_bounds :: Dict{Domain, SumPolyProdC})
  all_containmanated = find_all_mult_containmenates(f_src_bad.var_order, f_mult.var_order, Set(keys(bad_bounds)), f_mult.bsp)

  good_parent = (if f_src_bad == f_mult.parent1
                   f_mult.parent2
                 else
                   f_mult.parent1
                 end)

  # return the upper bound of a patch's good parent
  function get_good_parent_upper_bound(mult_dom :: Domain)
    shrink = diminish_dom_dim(f_mult.var_order, good_parent.var_order, mult_dom)
    best_cover = best_covering(shrink, keys(good_parent.potential_bounds))
    good_parent.potential_bounds[best_cover].p_u
  end

  ret = Dict{Domain, SumPolyProdC}()
  for bad_mult_dom in all_containmanated
    ret[bad_mult_dom] = get_good_parent_upper_bound(bad_mult_dom) * bad_bounds
  end

end

# ============== DUMPSTER ===============
# make a factor by having a different partition (maybe coarser)
function abstract_factor(f :: Factor, p :: Partition, bnd_method="exact", deg=2)
  lower_fun(x...) = feval_lower(f, [x...])
  upper_fun(x...) = feval_upper(f, [x...])
  var_order = f.var_order

  potential_bounds = Dict{Domain, (SumPolyProdC, SumPolyProdC)}()
  for dom in p
    # the two spp approximations
    spp_approx_lower = get_m_projections_approx(lower_fun, var_order, deg, dom)
    spp_approx_upper = get_m_projections_approx(upper_fun, var_order, deg, dom)

    # try to figure out how to shift those 2 approximations to get sound results
    intersected_doms = filter(x->has_intersect(x,dom), keys(f.potential_bounds))
    shift_down = 0.0
    shift_up = 0.0

    # for all the original factor's dom that intersect with this dom in p
    for x in intersected_doms
      # get the intersection dom of the 2 overlapping doms
      common_dom = get_intersect(x, dom)
      # compute the diff of the 2 polynomial bounds for each lower / upper bound
      shift_down_poly = spp_approx_lower - f.potential_bounds[x][1]
      shift_up_poly = f.potential_bounds[x][2] - spp_approx_upper
      # and turn them into their numerical versions
      shift_down_poly_enum(x) = peval(shift_down_poly, [x...])
      shift_up_poly_enum(x) = eval(shift_up_poly, [x...])
      # compute now how much needs to be shifted to perserve soundness
      # compute it over the common dom where both polynomials share
      cur_shift_down = (if bnd_method == "exact"
                          get_max_value_exact(shift_down_poly, common_dom)
                        else
                          get_max_value_approx(shift_down_poly_enum, common_dom)
                        end
                        )
      shift_down = max(shift_down, cur_shift_down)
      cur_shift_up = (if bnd_method == "exact"
                        get_max_value_exact(shift_down_poly, common_dom)
                      else
                        get_max_value_approx(shift_down_poly_enum, common_dom)
                      end
                      )
      shift_up = max(shift_up, cur_shift_up)
    end
    lower_poly = SumPolyProdC(spp_approx_lower, -1.0 * shift_down)
    upper_poly = SumPolyProdC(spp_approx_upper, shift_up)

    lowest_of_lower = lowest_possible(lower_poly, dom, bnd_method)
    lower_poly = (if lowest_of_lower < 0.0
                    make_zero(lower_poly)
                  else
                    lower_poly
                  end
                  )

    potential_bounds[dom] = lower_poly, upper_poly
  end
  Factor(var_order, potential_bounds)
end


function * (f1 :: Factor, f2 :: Factor)
  var_order = merge_varorder(f1.var_order, f2.var_order)
  potential_bounds = Dict{Domain, (SumPolyProdC, SumPolyProdC)}()
  for dom1 in keys(f1.potential_bounds)
    for dom2 in keys(f2.potential_bounds)
      dom1_enlarged = enlarge_dom_dim(f1.var_order, var_order, dom1)
      dom2_enlarged = enlarge_dom_dim(f2.var_order, var_order, dom2)
      if has_intersect(dom1_enlarged, dom2_enlarged)
        dom_common = get_intersect(dom1_enlarged, dom2_enlarged)
        lower1, upper1 = f1.potential_bounds[dom1]
        lower2, upper2 = f2.potential_bounds[dom2]
        potential_bounds[dom_common] = (lower1*lower2, upper1*upper2)
      end
    end
  end
  Factor(var_order, potential_bounds)
end

function ∫ (f :: Factor, inte_var_name :: ASCIIString)
  var_order = filter(x->x!=inte_var_name, f.var_order)
  # project the partition on f into a lower dimension while keeping track of where each new domain in the partition
  # is an intersection of what in the old place
  projected_partition_map = partition_projection(f.var_order, var_order, inte_var_name, keys(f.potential_bounds))
  # we attempt to integrate the lower/upper bounds from all the original partitions, and cache the results
  cached_int_potential_bounds = Dict{Domain, (SumPolyProdC, SumPolyProdC)}()
  for s_dom in keys(f.potential_bounds)
    lower, upper = f.potential_bounds[s_dom]
    bound_left, bound_right = s_dom[findfirst(f.var_order, inte_var_name)]
    lower_poly_bnd = ∫(lower, inte_var_name, bound_left, bound_right)
    upper_poly_bnd = ∫(upper, inte_var_name, bound_left, bound_right)
    cached_int_potential_bounds[s_dom] = (lower_poly_bnd, upper_poly_bnd)
  end

  potential_bounds = Dict{Domain, (SumPolyProdC, SumPolyProdC)}()
  for intersected_dom in keys(projected_partition_map)
    dom_sources = Domain[projected_partition_map[intersected_dom]...]
    if (length(dom_sources)) == 1
      potential_bounds[intersected_dom] = cached_int_potential_bounds[dom_sources[1]]
    else
      inted_lower_uppers = (SumPolyProdC, SumPolyProdC)[cached_int_potential_bounds[dom_source] for dom_source in dom_sources]
      comebine_pair(l1_u1::(SumPolyProdC,SumPolyProdC), l2_u2::(SumPolyProdC,SumPolyProdC)) = (l1_u1[1]+l2_u2[1], l1_u1[2]+l2_u2[2])
      potential_bounds[intersected_dom] = reduce(comebine_pair, inted_lower_uppers)
    end
  end
  Factor(var_order, potential_bounds)
end
