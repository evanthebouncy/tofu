include("/home/evan/Documents/research/tofu/Factor.jl")
PROP_AMOUNT = 0.9
Gdict = Dict{Any, Any}()

# the baseline random strategy for growth
function rand_grow!(FG :: FactorGraph)
  len = length(FG.factors)
  rand_id = convert(Int64, round((len-1) * rand())) + 1

  f = FG.factors[rand_id]
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  factor_grow!(FG, f, dom_to_grow)
end


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
function get_contributing_partitions(FG :: FactorGraph, upstreams_contribution :: Dict{Factor, Set{Domain}})
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

function get_valid_split_domains(FG::FactorGraph, upstreams_contribution::Dict{Factor, Set{Domain}})

  function get_parent_doms(dom::Domain, f_type::FactorOp)
    # if it is multiplication...
    if typeof(f_type) == FactorMult
      f1 = f_type.f1
      f2 = f_type.f2
      f_mult = f_type.f_mult
      shrink1 = diminish_dom_dim(f_mult.var_order, f1.var_order, dom)
      shrink2 = diminish_dom_dim(f_mult.var_order, f2.var_order, dom)
      best_cover1 = best_covering(shrink1, keys(f1.potential_bounds))
      best_cover2 = best_covering(shrink2, keys(f2.potential_bounds))
      return [f1=>Set{Domain}(Domain[best_cover1]),
       f2=>Set{Domain}(Domain[best_cover2])]
    end
    # if it is integration...
    if typeof(f_type) == FactorInte
      f1 = f_type.f1
      inte_var_name = f_type.inte_var
      factor = f_type.f_inte
      smallest_cover = find_smallest_cover(f1.bsp, dom, factor.var_order, f1.var_order, inte_var_name)
      return [f1=>Set(smallest_cover)]
    end
  end

  # first get the contributing partitions...
  contrib_partition = get_contributing_partitions(FG, upstreams_contribution)

  ret = Set{(Factor,Domain)}()
  for f_dom in contrib_partition
    factor, dom = f_dom
    # println("my name is ", factor.f_name)
    # if it has no parents, it is a potential, shattering it will for sure
    # increase precision, so we add it already
    if !(factor in keys(FG.rel_factor_parents))
      push!(ret, f_dom)
    else
      # otherwise, we attempt to shatter it first...
      shattered_doms = shatter(dom)
      # get the parents for the factor...
      f_type = FG.rel_factor_parents[factor]
      shattered_parents = Dict{Factor, Set{Domain}}[get_parent_doms(shat_d, f_type) for shat_d in shattered_doms]
      # for debugging
      # Gdict["sp$(factor.f_name)"] = shattered_parents
      # println("showing parents ", shattered_parents)
      #end debugging
      is_sames = [shattered_parents[i] == shattered_parents[i+1] for i in 1:(length(shattered_parents)-1)]
      # if there exists a difference...
      if false in is_sames
        push!(ret, f_dom)
      end
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


# =================== Diagnostics ======================
function diagnostic_f_imprecision(FG, f, dom)
  println("fuck")
  @show(f.f_name, dom)

end

##### drawing a 2d domain
function draw_dom2d(doms)
  doms = [d for d in doms]
  color = [rand() * 1000 for i in doms]
  x_min = Float64[x[1][1] for x in doms]
  x_max = Float64[x[1][2] for x in doms]
  y_min = Float64[x[2][1] for x in doms]
  y_max = Float64[x[2][2] for x in doms]
  plot(x_max = x_max, x_min = x_min, y_max = y_max, y_min = y_min,
       color=color,
       Geom.rectbin)
end

function draw_dom1d(doms)
  doms = [d for d in doms]
  color = [rand() * 1000 for i in doms]
  x_min = Float64[x[1][1] for x in doms]
  x_max = Float64[x[1][2] for x in doms]
  y_min = Float64[0.0 for x in doms]
  y_max = Float64[1.0 for x in doms]
  plot(x_max = x_max, x_min = x_min, y_max = y_max, y_min = y_min,
       color=color,
       Geom.rectbin)
end

function draw_f_imprecision2d(FG, f)
  doms = [d for d in f.partition]
  color = [
   (   if (f, dom) in keys(FG.memoized_cost)
    FG.memoized_cost[(f, dom)]
  else
    0.0
  end)

    for dom in doms]
  x_min = Float64[x[1][1] for x in doms]
  x_max = Float64[x[1][2] for x in doms]
  y_min = Float64[x[2][1] for x in doms]
  y_max = Float64[x[2][2] for x in doms]
  plot(x_max = x_max, x_min = x_min, y_max = y_max, y_min = y_min,
       color=color,
       Geom.rectbin)
end

function draw_f_imprecision1d(FG, f)
  doms = [d for d in f.partition]
  color = [
   (   if (f, dom) in keys(FG.memoized_cost)
    FG.memoized_cost[(f, dom)]
  else
    0.0
  end)

    for dom in doms]
  x_min = Float64[x[1][1] for x in doms]
  x_max = Float64[x[1][2] for x in doms]
  y_min = Float64[0.0 for x in doms]
  y_max = Float64[1.0 for x in doms]
  plot(x_max = x_max, x_min = x_min, y_max = y_max, y_min = y_min,
       color=color,
       Geom.rectbin)
end
