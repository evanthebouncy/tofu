include("/home/evan/Documents/research/tofu/Abstraction.jl")
import Base.isless
PROP_AMOUNT = 0.9
Gdict = Dict{Any, Any}()

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

# factors are un-ordered by default ...
isless(f1::Factor, f2::Factor) = true

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
  # remember the potentials with the source factors
  potentials :: Dict{Factor, Potential}
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
  # memoize the cost of the partitions
  memoized_cost :: Dict{(Factor, Domain), Float64}
  # keep track of if memozied_cost has been recently computed
  cost_recent :: Set{(Factor, Domain)}
  # keep track of a list of integrals recently computed
  integral_recent :: Dict{(SumPolyProdC, ASCIIString, Float64, Float64), SumPolyProdC}
end

# ================================ methods to modify the factor graph and its relations ==================================

# initialize a factor graph
function init_factor_graph()
  FactorGraph(Dict{Factor, Potential}(), Factor[], Dict{Factor,FactorOp}(), Dict{Factor,FactorOp}(),
              Dict{(Factor,Domain), DomainOp}(), Dict{(Factor,Domain),Array{DomainOp}}(),
              Dict{(Factor, Domain), Float64}(), Set{(Factor, Domain)}(),
              Dict{(SumPolyProdC, ASCIIString, Float64, Float64), SumPolyProdC}() )
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
function patch_pot(FG :: FactorGraph, var_order, dom :: Domain, pot :: Potential, bnd_method="approx", deg=2)
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
  # book keeping
  add_factor!(FG, ret)
  FG.potentials[ret] = pot
  ret
end

# grow a factor from a potential function
function f_pot_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)
  pot = FG.potentials[f]
  # add the patches for the shattered domains
  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_pot(FG, f.var_order, d1, pot, "approx", 2)
  end
end

function rand_f_pot_grow!(FG :: FactorGraph)
  len = length(FG.factors)
  rand_id = convert(Int64, round((len-1) * rand())) + 1

  f = FG.factors[rand_id]
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  factor_grow!(FG, f, dom_to_grow)
end

# give a domain, and 2 factors f1 and f2
function patch_mult (FG :: FactorGraph, var_order, dom :: Domain, f1 :: Factor, f2 :: Factor, f_mult :: Factor)
  shrink1 = diminish_dom_dim(var_order, f1.var_order, dom)
  shrink2 = diminish_dom_dim(var_order, f2.var_order, dom)
  best_cover1 = best_covering(shrink1, keys(f1.potential_bounds))
  best_cover2 = best_covering(shrink2, keys(f2.potential_bounds))
  patch1, patch2 = f1.potential_bounds[best_cover1], f2.potential_bounds[best_cover2]
  p_l = sppc_mult(patch1.p_l, patch2.p_l)
  p_u = sppc_mult(patch1.p_u, patch2.p_u)
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
  # register fact relationship to FG and various book keepings
  push!(FG.factors, ret)
  add_rel_factor!(FG, FactorMult(f1, f2, ret))
  ret.potential_bounds[cover_dom] = patch_mult(FG, var_order, cover_dom, f1, f2, ret)
  ret
end

# grow a factor from a multiplication
function f_mult_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain)
  parents_rel = FG.rel_factor_parents[f]
  f1, f2 = parents_rel.f1, parents_rel.f2
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_mult(FG, f.var_order, d1, f1, f2, f)
  end
end

function rand_f_mult_grow!(FG :: FactorGraph, f :: Factor)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_mult_grow!(FG, f, dom_to_grow)
end

function ∫_memo (FG::FactorGraph, sppc :: SumPolyProdC, x :: ASCIIString, a, b)
  if (sppc, x, a, b) in FG.integral_recent
    FG.integral_recent[(sppc, x, a, b)]
  else
    ret = ∫(sppc, x, a, b)
    FG.integral_recent[(sppc, x, a, b)] = ret
    ret
  end
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
    upper_poly_bnd = ∫_memo(FG, upper, inte_var_name, bound_left, bound_right)
    push!(integrated_pots, (lower_poly_bnd, upper_poly_bnd))
  end
  lower_pol, upper_pol = (if length(integrated_pots) == 1
                            integrated_pots[1]
                          else
                            reduce((x,y)-> (x[1]+y[1],x[2]+y[2]), integrated_pots)
                          end
                          )
  # register relationship in FG and book keepings
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
function f_inte_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain)
  parents_rel = FG.rel_factor_parents[f]
  f1, inte_var_name = parents_rel.f1, parents_rel.inte_var
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_integrate(FG, f.var_order, d1, f1, inte_var_name, f)
  end
end

function rand_f_inte_grow!(FG :: FactorGraph, f :: Factor)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_inte_grow!(FG, f, dom_to_grow)
end

function factor_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain)
  # if it is a source potential factor
  if f in keys(FG.potentials)
    f_pot_grow!(FG, f, dom)
  # if it is not a source, it has parents
  else
    f_parent_rel = FG.rel_factor_parents[f]
    # if it is a multiplication ...
    if typeof(f_parent_rel) == FactorMult
      f_mult_grow!(FG, f, dom)
      # otherwise it is integration
    else
      f_inte_grow!(FG, f, dom)
    end
  end
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
    return [f1 => dom1s, f2 => dom2s]
  end
  # if it is integration...
  if typeof(f_type) == FactorInte
    f1 = f_type.f1
    dom1s = Set{Domain}()
    for dom in doms
      rel_inte = FG.rel_domain_parents[(factor, dom)]
      union!(dom1s, rel_inte.d1s)
    end
    return [f1 => dom1s]
  end
end

function all_upstream_contributions_rec(factor :: Factor, doms :: Set{Domain})
  # println("upstream contribution of ", factor.f_name)
  attempt_upstream = one_step_upstream_contribution(factor, doms)
  if (attempt_upstream == Nothing())
    [factor => doms]
  else
    ret = Dict{Factor, Set{Domain}}()
    for factor1 in keys(attempt_upstream)
      ret = merge(ret, all_upstream_contributions_rec(factor1, attempt_upstream[factor1]))
    end
    ret = merge(ret, [factor => doms])
    ret
  end
end

# find out all the domains that could be further shattered with effect from upstream contributions
function get_valid_split_domains(FG :: FactorGraph, upstreams_contribution :: Dict{Factor, Set{Domain}})
  ret = Set{(Factor, Domain)}()
  for fact in keys(upstreams_contribution)
    doms = upstreams_contribution[fact]
    filtered_doms = filter(d->d in fact.partition, doms)
    for filt_dom in filtered_doms
      push!(ret, (fact, filt_dom))
    end
  end
  ret
end

# get imprecision on a single factor and a single domain
# if the result is cached, return it, if not, call the expensife funcion
# and cache the result
function get_imprecision!(FG :: FactorGraph ,f :: Factor, dom :: Domain, contributors :: Dict{Factor, Set{Domain}})
  if (f, dom) in FG.cost_recent
    FG.memoized_cost[(f, dom)]
  else
    ret_imprecision = get_set_imprecisions!(FG, f, [dom => f.potential_bounds[dom].p_l], contributors, 1.0)
    FG.memoized_cost[(f, dom)] = ret_imprecision
    ret_imprecision
  end
end

# given a factor and a set of domains in it, find out the imprecision
# of them all the way to the bottom
function get_set_imprecisions!(FG::FactorGraph, f::Factor, dom_lower_bnds::Dict{Domain, SumPolyProdC}, contributors::Dict{Factor, Set{Domain}}, relative_weight::Float64)
  if !(f in keys(FG.rel_factor_child))
    final_bnd, final_f = dom_lower_bnds, f
    dom_sppcs = [(final_dom, sppc_sub(final_f.potential_bounds[final_dom].p_u, final_bnd[final_dom])) for final_dom in keys(final_bnd)]
    volumes = [poly_volume(final_f.var_order, d_sppc[1], d_sppc[2]) for d_sppc in dom_sppcs]
    if length(volumes) == 0
      return 0.0
    end
    return reduce((x,y)->x+y, volumes) * relative_weight
  end

  factor_child_rel = FG.rel_factor_child[f]
  if typeof(factor_child_rel) == FactorMult
    nxt_f, nxt_bnd, nxt_relative_weight = propagate_mult_imprecise_lower(FG, f, factor_child_rel, dom_lower_bnds, contributors)
    get_set_imprecisions!(FG, nxt_f, nxt_bnd, contributors, relative_weight * nxt_relative_weight)
  else
    assert(typeof(factor_child_rel) == FactorInte)
    nxt_f, nxt_bnd, nxt_relative_weight = propagate_inte_imprecise_lower(FG, f, factor_child_rel, dom_lower_bnds, contributors)
    get_set_imprecisions!(FG, nxt_f, nxt_bnd, contributors, relative_weight * nxt_relative_weight)
  end
end

function propagate_mult_imprecise_lower(FG::FactorGraph, f::Factor, factor_child_rel::FactorMult, dom_lower_bnds::Dict{Domain, SumPolyProdC}, contributors::Dict{Factor, Set{Domain}})
  # the next factor is the f_mult field of the relation
  nxt_fact = factor_child_rel.f_mult
  # what domain in the next factor should we even consider?
  valid_children = Set{Domain}()
  # for all the domain in our input, that has imprecise lower bound
  for f_dom in keys(dom_lower_bnds)
    # how is it related to its children in the nxt_fact ? there are multiple relationships possible
    dom_children_rel = FG.rel_domain_children[(f, f_dom)]
    # for each of the possible relationship
    for dom_child_rel in dom_children_rel
      # get the domain of the child
      child_dom = dom_child_rel.d_mult
      # consider the child domain only if it is a contributor
      if child_dom in contributors[nxt_fact]
        push!(valid_children, child_dom)
      end
    end
  end
  # now we want to filter the valid children to include only a percentage of the imprecisions they occur!!
  filtered_children, relative_weight = filter_doms_by_imprecision(FG, nxt_fact, valid_children, contributors, PROP_AMOUNT)
  nxt_imprecise = Dict{Domain, SumPolyProdC}()
  # now we have a set of valid children dom that we know we want to consider... for each of them
  for c_dom in filtered_children
    parents_rel = FG.rel_domain_parents[nxt_fact, c_dom]
    factor_parent1 = parents_rel.fmultop.f1
    factor_parent2 = parents_rel.fmultop.f2
    poly1 = (if factor_parent1 == f
               dom_lower_bnds[parents_rel.d1]
             else
               factor_parent1.potential_bounds[parents_rel.d1].p_u
             end)
    poly2 = (if factor_parent2 == f
               dom_lower_bnds[parents_rel.d2]
             else
               factor_parent2.potential_bounds[parents_rel.d2].p_l
             end)
    nxt_imprecise[c_dom] = sppc_mult(poly1, poly2)
  end
  nxt_fact, nxt_imprecise, relative_weight
end

function propagate_inte_imprecise_lower(FG::FactorGraph, f::Factor, factor_child_rel::FactorInte, dom_lower_bnds::Dict{Domain, SumPolyProdC}, contributors::Dict{Factor, Set{Domain}})
  # the next factor is the f_inte field of the relation
  nxt_fact = factor_child_rel.f_inte
  # what domain in the next factor should we even consider?
  valid_children = Set{Domain}()
  # for all the domain in our input, that has imprecise lower bound
  for f_dom in keys(dom_lower_bnds)
    # how is it related to its children in the nxt_fact ? there are multiple relationships possible
    dom_children_rel = FG.rel_domain_children[(f, f_dom)]
    # for each of the possible relationship
    for dom_child_rel in dom_children_rel
      # get the domain of the child
      child_dom = dom_child_rel.d_inte
      # consider the child domain only if it is a contributor
      if child_dom in contributors[nxt_fact]
        push!(valid_children, child_dom)
      end
    end
  end
  # now we want to filter the valid children to include only a percentage of the imprecisions they occur!!
  filtered_children, relative_weight = filter_doms_by_imprecision(FG, nxt_fact, valid_children, contributors, PROP_AMOUNT)

  nxt_imprecise = Dict{Domain, SumPolyProdC}()
  # now we have a set of valid children dom that we know we want to consider... for each of them
  for c_dom in filtered_children
    # get the parents
    parents_rel = FG.rel_domain_parents[nxt_fact, c_dom]
    inte_var = parents_rel.finteop.inte_var
    # this is all the parent domain that gets integrated down
    dom1s = parents_rel.d1s
    integrated_pots = SumPolyProdC[]
    # for all the potential parents
    for s_dom in dom1s
      bound_left, bound_right = s_dom[findfirst(f.var_order, inte_var)]
      # if it is a parent that has been made imprecise...
      if s_dom in keys(dom_lower_bnds)
        push!(integrated_pots, ∫(dom_lower_bnds[s_dom], inte_var, bound_left, bound_right))
        # if not, use the upper bound
      else
        push!(integrated_pots, ∫_memo(FG, f.potential_bounds[s_dom].p_u, inte_var, bound_left, bound_right))
      end
    end
    imprecise_lower = (if length(integrated_pots) == 1
                         integrated_pots[1]
                       else
                         reduce((x,y)->x+y, integrated_pots)
                       end
                       )
    nxt_imprecise[c_dom] = imprecise_lower
  end
  nxt_fact, nxt_imprecise, relative_weight
end

# given a set of domain, put them in queue and pull fraction of them out based on their imprecisions
function filter_doms_by_imprecision(FG::FactorGraph, f::Factor, doms::Set{Domain}, contributors::Dict{Factor, Set{Domain}}, fraction::Float64)
  # the queue contains entries of the form ((factor, domain), cost) : cost
  cost_queue = Collections.PriorityQueue{((Factor, Domain), Float64), Float64}()
  total_cost =  0.0
  for dom in doms
    cost = get_imprecision!(FG, f, dom, contributors)
#     cost = (if (f, dom) in keys(FG.memoized_cost)
#               FG.memoized_cost[(f, dom)]
#             else
#               get_imprecision!(FG, f, dom, contributors)
#             end
#             )
    Collections.enqueue!(cost_queue, ((f,dom),cost), -cost)

    total_cost += cost
  end
  ret = Set{Domain}()
  cur_cost = 0.0
  while (cur_cost < (total_cost * fraction - 1e-6)) & (length(cost_queue) > 0)
    f_dom, d_cost = Collections.dequeue!(cost_queue)
    push!(ret, f_dom[2])
    cur_cost = cur_cost + d_cost
  end
  ret, cur_cost / total_cost
end



# pretend the other bounds are prefect i.e. lower == upper
# we propagate the imperfect lower bound
function propagate_imprecise_lower(FG :: FactorGraph, f :: Factor, dom_lower_bnds :: Dict{Domain, SumPolyProdC}, contributors :: Dict{Factor, Set{Domain}})
  if !(f in keys(FG.rel_factor_child))
    return (f, dom_lower_bnds)
  end
  # how is this factor f related to its children?
  factor_child_rel = FG.rel_factor_child[f]

  # in the case of multiplication
  if typeof(factor_child_rel) == FactorMult
    # the next factor is the f_mult field of the relation
    nxt_fact = factor_child_rel.f_mult
    # what domain in the next factor should we even consider?
    valid_children = Set{Domain}()
    # for all the domain in our input, that has imprecise lower bound
    for f_dom in keys(dom_lower_bnds)
      # how is it related to its children in the nxt_fact ? there are multiple relationships possible
      dom_children_rel = FG.rel_domain_children[(f, f_dom)]
      # for each of the possible relationship
      for dom_child_rel in dom_children_rel
        # get the domain of the child
        child_dom = dom_child_rel.d_mult
        # consider the child domain only if it is a contributor
        if child_dom in contributors[nxt_fact]
          push!(valid_children, child_dom)
        end
      end
    end
    nxt_imprecise = Dict{Domain, SumPolyProdC}()
    # now we have a set of valid children dom that we know we want to consider... for each of them
    for c_dom in valid_children
      parents_rel = FG.rel_domain_parents[nxt_fact, c_dom]
      factor_parent1 = parents_rel.fmultop.f1
      factor_parent2 = parents_rel.fmultop.f2
      poly1 = (if factor_parent1 == f
                 dom_lower_bnds[parents_rel.d1]
               else
                 factor_parent1.potential_bounds[parents_rel.d1].p_u
               end)
      poly2 = (if factor_parent2 == f
                 dom_lower_bnds[parents_rel.d2]
               else
                 factor_parent2.potential_bounds[parents_rel.d2].p_l
               end)
      nxt_imprecise[c_dom] = sppc_mult(poly1, poly2)
    end
    return nxt_fact, nxt_imprecise
  end

  # in the case of integration
  if typeof(factor_child_rel) == FactorInte
    # the next factor is the f_mult field of the relation
    nxt_fact = factor_child_rel.f_inte
    # what domain in the next factor should we even consider?
    valid_children = Set{Domain}()
    # for all the domain in our input, that has imprecise lower bound
    for f_dom in keys(dom_lower_bnds)
      # how is it related to its children in the nxt_fact ? there are multiple relationships possible
      dom_children_rel = FG.rel_domain_children[(f, f_dom)]
      # for each of the possible relationship
      for dom_child_rel in dom_children_rel
        # get the domain of the child
        child_dom = dom_child_rel.d_inte
        # consider the child domain only if it is a contributor
        if child_dom in contributors[nxt_fact]
          push!(valid_children, child_dom)
        end
      end
    end
    nxt_imprecise = Dict{Domain, SumPolyProdC}()
    # now we have a set of valid children dom that we know we want to consider... for each of them
    for c_dom in valid_children
      # get the parents
      parents_rel = FG.rel_domain_parents[nxt_fact, c_dom]
      inte_var = parents_rel.finteop.inte_var
      # this is all the parent domain that gets integrated down
      dom1s = parents_rel.d1s
      integrated_pots = SumPolyProdC[]
      # for all the potential parents
      for s_dom in dom1s
        bound_left, bound_right = s_dom[findfirst(f.var_order, inte_var)]
        # if it is a parent that has been made imprecise...
        if s_dom in keys(dom_lower_bnds)
          push!(integrated_pots, ∫(dom_lower_bnds[s_dom], inte_var, bound_left, bound_right))
        # if not, use the upper bound
        else
          push!(integrated_pots, ∫(f.potential_bounds[s_dom].p_u, inte_var, bound_left, bound_right))
        end
      end
      imprecise_lower = (if length(integrated_pots) == 1
                           integrated_pots[1]
                         else
                           reduce((x,y)->x+y, integrated_pots)
                         end
                         )
      nxt_imprecise[c_dom] = imprecise_lower
    end
    return nxt_fact, nxt_imprecise
  end

end

function propagate_imprecise_lower_rec(FG :: FactorGraph, f :: Factor, dom_lower_bnds :: Dict{Domain, SumPolyProdC}, contributors :: Dict{Factor, Set{Domain}})
  nxt_f, nxt_bnd = propagate_imprecise_lower(FG, f, dom_lower_bnds, contributors)
  if nxt_f == f
    nxt_f, nxt_bnd
  else
    propagate_imprecise_lower_rec(FG , nxt_f, nxt_bnd, contributors)
  end
end

function find_imprecise_error(FG :: FactorGraph, f :: Factor, dom_lower_bnds :: Dict{Domain, SumPolyProdC}, contributors :: Dict{Factor, Set{Domain}})
  # use the recursion thing to find the final factor and the imprecise lower bound of it
  final_f, final_bnd = propagate_imprecise_lower_rec(FG , f, dom_lower_bnds, contributors)
  # for all the domain that has imprecise lower bound, compute p_u - p_l on it...
  dom_sppcs = [(final_dom, sppc_sub(final_f.potential_bounds[final_dom].p_u, final_bnd[final_dom])) for final_dom in keys(final_bnd)]
  volumes = [poly_volume(final_f.var_order, d_sppc[1], d_sppc[2]) for d_sppc in dom_sppcs]
  reduce((x,y)->x+y, volumes)
end

function find_best_split!(FG :: FactorGraph)
  last_factor = last(FG.factors)
  upstream_contribution = all_upstream_contributions_rec(last_factor, last_factor.partition)
  all_split_doms = get_valid_split_domains(FG, upstream_contribution)

  # the queue contains entries of the form ((factor, domain), cost) : cost
  cost_queue = Collections.PriorityQueue{((Factor, Domain), Float64), Float64}()
  for fact_dom in all_split_doms
    fact, dom = fact_dom[1], fact_dom[2]
    if fact_dom in keys(FG.memoized_cost)
      Collections.enqueue!(cost_queue, ((fact, dom), FG.memoized_cost[fact_dom]), -FG.memoized_cost[fact_dom])
    else
      imprecise_err = find_imprecise_error(FG, fact, [dom => fact.potential_bounds[dom].p_l], upstream_contribution)
      Collections.enqueue!(cost_queue, ((fact, dom), imprecise_err), -imprecise_err)
    end
  end

  # pop the largest element, and re-computes its cost
  cur_f_dom, cur_cost = Collections.dequeue!(cost_queue)
  cur_f, cur_dom = cur_f_dom
  recomputed_cost = find_imprecise_error(FG, cur_f, [cur_dom => cur_f.potential_bounds[cur_dom].p_l], upstream_contribution)
  @show(cur_cost, recomputed_cost)
  # if the recomputed cost is the same, for SURE we want to split this one!! so we exit the loop
  # however if the recomputed cost is smaller... then we have to update
  stuck = 0
  while recomputed_cost < (cur_cost - 1e-6)
    stuck = stuck + 1
    @show(stuck)
    FG.memoized_cost[(cur_f, cur_dom)] = recomputed_cost
    Collections.enqueue!(cost_queue, ((cur_f, cur_dom), recomputed_cost), -recomputed_cost)
  end
  cur_f, cur_dom
end

function find_best_split_lazy!(FG :: FactorGraph)
  last_factor = last(FG.factors)
  upstream_contribution = all_upstream_contributions_rec(last_factor, last_factor.partition)
  Gdict["upstream_controbution"] = upstream_contribution
  all_split_doms = get_valid_split_domains(FG, upstream_contribution)

  Gdict["all_split_doms"] = all_split_doms

  # the queue contains entries of the form ((factor, domain), cost) : cost
  cost_queue = Collections.PriorityQueue{((Factor, Domain), Float64), Float64}()
  for fact_dom in all_split_doms
    fact, dom = fact_dom[1], fact_dom[2]
    if fact_dom in keys(FG.memoized_cost)
      Collections.enqueue!(cost_queue, ((fact, dom), FG.memoized_cost[fact_dom]), -FG.memoized_cost[fact_dom])
    else
      imprecise_err = get_imprecision!(FG, fact, dom, upstream_contribution)
      Collections.enqueue!(cost_queue, ((fact, dom), imprecise_err), -imprecise_err)
    end
  end

  # pop the largest element, and re-computes its cost
  cur_f_dom, cur_cost = Collections.dequeue!(cost_queue)
  cur_f, cur_dom = cur_f_dom
  recomputed_cost = get_imprecision!(FG, cur_f, cur_dom, upstream_contribution)
  # if the recomputed cost is the same, for SURE we want to split this one!! so we exit the loop
  # however if the recomputed cost is smaller... then we have to update
  stuck = 0
  while recomputed_cost < (cur_cost - 1e-6)
    stuck = stuck + 1
    @show(stuck)
    FG.memoized_cost[(cur_f, cur_dom)] = recomputed_cost
    Collections.enqueue!(cost_queue, ((cur_f, cur_dom), recomputed_cost), -recomputed_cost)

    # pop the largest element, and re-computes its cost
    cur_f_dom, cur_cost = Collections.dequeue!(cost_queue)
    cur_f, cur_dom = cur_f_dom
    recomputed_cost = get_imprecision!(FG, cur_f, cur_dom, upstream_contribution)
  end
  println("cost is ", recomputed_cost)
  cur_f, cur_dom
end

function heuristic_grow!(FG :: FactorGraph)
  best_fact, best_dom = find_best_split_lazy!(FG)
  println(best_fact.f_name, "    ", best_dom)
  FG.cost_recent = Set{(Factor, Domain)}()
  FG.integral_recent = Dict{(SumPolyProdC, ASCIIString, Float64, Float64), SumPolyProdC}()

  factor_grow!(FG, best_fact, best_dom)
end
