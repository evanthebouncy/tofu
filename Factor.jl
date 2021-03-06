include("/home/evan/Documents/research/tofu/Abstraction.jl")
import Base.isless
import Base.==
import Base.show

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
function show(io::IO, f::Factor)
    print(io,f.f_name)
end

# this keeps relationships between factors as operations Mult, and Inte
abstract FactorOp

immutable FactorMult <: FactorOp
  f1 :: Factor
  f2 :: Factor
  f_mult :: Factor
end
function show(io::IO, fM::FactorMult)
  print(io,"FactorMult($(fM.f1.f_name), $(fM.f2.f_name), $(fM.f_mult.f_name))")
end

immutable FactorInte <: FactorOp
  f1 :: Factor
  inte_var :: ASCIIString
  f_inte :: Factor
end

# this is analegous relationship
abstract DomainOp

immutable DomainMult <: DomainOp
  fmultop :: FactorMult
  d1 :: Domain
  d2 :: Domain
  d_mult :: Domain
end

function == (dop1::DomainMult, dop2::DomainMult)
  (dop1.fmultop == dop2.fmultop) &
    (dop1.d1 == dop2.d1) &
    (dop1.d2 == dop2.d2) &
    (dop1.d_mult == dop2.d_mult)
end

immutable DomainInte <: DomainOp
  finteop :: FactorInte
  d1s :: Array{Domain}
  d_inte :: Domain
end

function == (dop1::DomainInte, dop2::DomainInte)
  (dop1.finteop == dop2.finteop) &
    (dop1.d1s == dop2.d1s) &
    (dop1.d_inte == dop2.d_inte)
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

function remove_rel_domain_children!(FG::FactorGraph, f::Factor, d::Domain, dop::DomainOp)
  key = (f, d)
  old_children_list = FG.rel_domain_children[key]
  FG.rel_domain_children[key] = filter(x->x!=dop, old_children_list)
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

# remove a relation domainOP from the factor graph
function remove_rel_domain!(FG::FactorGraph, d_op::DomainOp)
  if typeof(d_op) == DomainMult
    delete!(FG.rel_domain_parents, (d_op.fmultop.f_mult, d_op.d_mult))
    remove_rel_domain_children!(FG, d_op.fmultop.f1, d_op.d1, d_op)
    remove_rel_domain_children!(FG, d_op.fmultop.f2, d_op.d2, d_op)
    return
  end
  if typeof(d_op) == DomainInte
    delete!(FG.rel_domain_parents, (d_op.finteop.f_inte, d_op.d_inte))
    for p_domain in d_op.d1s
      remove_rel_domain_children!(FG, d_op.finteop.f1, p_domain, d_op)
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
function shatter(dom::Domain)
    function split_once(dom :: Domain)
    sp, d1, d2 = split_half(dom)
    Domain[d1, d2]
  end

  function shatter_rec(doms :: Array{Domain}, n)
    if n == 0
      doms
    else
      new_doms = Domain[]
      for d in doms
        new_doms = [split_once(d), new_doms]
      end
      shatter_rec(new_doms, n - 1)
    end
  end

  shatter_rec(Domain[dom], length(dom))
end


# shattering while modifying some data structures
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

  #shatter!_rec(var_order, Domain[dom], p, bsp, length(var_order))
  shatter!_rec(var_order, Domain[dom], p, bsp, 1)
end

# rewire the children parent relationship. Given a parent that's been shattered, rewire its children for better fit
# also computes a set of better bounds for the children's patches, and give a set of domains (which had altered patches) upon return
function rewire!(FG::FactorGraph, p_fact::Factor, p_dom::Domain, shattered_doms::Array{Domain})
  old_children_rels = FG.rel_domain_children[(p_fact, p_dom)]
  # erase all the old childrens relation from FG
  for old_child_rel in old_children_rels
    remove_rel_domain!(FG, old_child_rel)
  end
  # repopulate the child_relation by recomputing the patch
  type1 = old_children_rels[1]
  if typeof(type1) == DomainMult
    ret = Set{Domain}()
    fact_rel_type = type1.fmultop
    f_mult, f1, f2 = fact_rel_type.f_mult, fact_rel_type.f1, fact_rel_type.f2
    for old_child_rel in old_children_rels
      dom_mult = old_child_rel.d_mult
      new_patch = patch_mult(FG, f_mult.var_order, dom_mult, f1, f2, f_mult)
      f_mult.potential_bounds[dom_mult] = new_patch
      push!(ret, dom_mult)
    end
    return f_mult, ret
  end
  if typeof(type1) == DomainInte
    ret = Set{Domain}()
    fact_rel_type = type1.finteop
    f1, inte_var, f_inte = fact_rel_type.f1, fact_rel_type.inte_var, fact_rel_type.f_inte
    for old_child_rel in old_children_rels
      dom_inte = old_child_rel.d_inte
      new_patch = patch_integrate(FG, f_inte.var_order, dom_inte, f1, inte_var, f_inte)
      f_inte.potential_bounds[dom_inte] = new_patch
      push!(ret, dom_inte)
    end
    return f_inte, ret
  end
end

# recursively propagate a set of (better) bounds by re-evaluating the patches of the domains
# and get the next set of domains too.
function propagate_bounds!(FG::FactorGraph, prev_fact::Factor, prev_changed_doms::Set{Domain})
  cur_rel = FG.rel_factor_child[prev_fact]
  if typeof(cur_rel) == FactorMult
    f1, f2, f_cur = cur_rel.f1, cur_rel.f2, cur_rel.f_mult
    changed_doms = Set{Domain}()
    for prev_dom in prev_changed_doms
      if (prev_fact, prev_dom) in keys(FG.rel_domain_children)
        children_rels = FG.rel_domain_children[(prev_fact, prev_dom)]
        for c in children_rels
          push!(changed_doms, c.d_mult)
        end
      end
    end
    for changed_dom in changed_doms
      new_patch = patch_mult(FG, f_cur.var_order, changed_dom, f1, f2, f_cur, false)
      f_cur.potential_bounds[changed_dom] = new_patch
    end
    return f_cur, changed_doms
  end

  if typeof(cur_rel) == FactorInte
    f1, inte_var, f_cur = cur_rel.f1, cur_rel.inte_var, cur_rel.f_inte
    changed_doms = Set{Domain}()
    for prev_dom in prev_changed_doms
      if (prev_fact, prev_dom) in keys(FG.rel_domain_children)
        children_rels = FG.rel_domain_children[(prev_fact, prev_dom)]
        for c in children_rels
          push!(changed_doms, c.d_inte)
        end
      end
    end
    for changed_dom in changed_doms
      new_patch = patch_integrate(FG, f_cur.var_order, changed_dom, f1, inte_var, f_cur, false)
      f_cur.potential_bounds[changed_dom] = new_patch
    end
    return f_cur, changed_doms
  end
end

function propagate_bounds_rec!(FG::FactorGraph, fact::Factor, changed_doms::Set{Domain})
  if fact == last(FG.factors)
    return
  else
    nxt_f, nxt_changed_doms = propagate_bounds!(FG, fact, changed_doms)
    return propagate_bounds_rec!(FG, nxt_f, nxt_changed_doms)
  end
end

# make a patch from a potential function
function patch_pot(FG :: FactorGraph, var_order, dom :: Domain, pot :: Potential, bnd_method="approx", deg=2)
  lower_poly, upper_poly = get_poly_lower_upper(pot, var_order, dom, deg, bnd_method)
  lowest_of_lower = lowest_possible(lower_poly, dom, bnd_method)
  lower_poly = (if lowest_of_lower < 0.0
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
  shattered_doms
end

# give a domain, and 2 factors f1 and f2
function patch_mult (FG :: FactorGraph, var_order, dom :: Domain, f1 :: Factor, f2 :: Factor, f_mult :: Factor, mutate=true)
  shrink1 = diminish_dom_dim(var_order, f1.var_order, dom)
  shrink2 = diminish_dom_dim(var_order, f2.var_order, dom)
  best_cover1 = best_covering(shrink1, f1.bsp)
  best_cover2 = best_covering(shrink2, f2.bsp)
  patch1, patch2 = f1.potential_bounds[best_cover1], f2.potential_bounds[best_cover2]
  p_l = sppc_mult(patch1.p_l, patch2.p_l)
  p_u = sppc_mult(patch1.p_u, patch2.p_u)
  # register relationship in FG
  if mutate
    fact_mult = FactorMult(f1, f2, f_mult)
    domain_mult = DomainMult(fact_mult, best_cover1, best_cover2, dom)
    add_rel_domain!(FG, domain_mult)
  end
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
  shattered_doms
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
function patch_integrate (FG :: FactorGraph, var_order, dom :: Domain, f1 :: Factor, inte_var_name, f_inte :: Factor, mutate=true)
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
  if mutate
    fact_inte = FactorInte(f1, inte_var_name, f_inte)
    domain_inte = DomainInte(fact_inte, [s_dom for s_dom in smallest_cover], dom)
    add_rel_domain!(FG, domain_inte)
  end
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
  shattered_doms
end

function factor_grow!(FG :: FactorGraph, f :: Factor, dom :: Domain)
  shattered_doms = Domain[]
  # if it is a source potential factor
  if f in keys(FG.potentials)
    shattered_doms = f_pot_grow!(FG, f, dom)
  # if it is not a source, it has parents
  else
    f_parent_rel = FG.rel_factor_parents[f]
    # if it is a multiplication ...
    if typeof(f_parent_rel) == FactorMult
      shattered_doms = f_mult_grow!(FG, f, dom)
      # otherwise it is integration
    else
      shattered_doms = f_inte_grow!(FG, f, dom)
    end
  end
  # if it is not the result
  if f != last(FG.factors)
    nxt_f, new_doms = rewire!(FG, f, dom, shattered_doms)
    propagate_bounds_rec!(FG, nxt_f, new_doms)
  end
end


