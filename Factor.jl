include("/home/evan/Documents/research/tofu/Abstraction.jl")
type Patch
  dom :: Domain
  p_l :: SumPolyProdC
  p_u :: SumPolyProdC
  parents :: Array{Domain}
end

type Factor
  var_order
  cover_dom :: Domain
  partition :: Partition
  potential_bounds :: Dict{Domain, Patch}
end

function feval_lower(f :: Factor, x)
  for dom in f.partition
    if dom_contains(dom, x)
      return peval(f.potential_bounds[dom].p_l, x)
    end
  end
  assert(false)
end

function feval_upper(f :: Factor, x)
  for dom in f.partition
    if dom_contains(dom, x)
      return peval(f.potential_bounds[dom].p_u, x)
    end
  end
  assert(false)
end

# make a patch from a potential function
function patch_pot(var_order, dom :: Domain, pot :: Potential, bnd_method, deg)
  lower_poly, upper_poly = get_poly_lower_upper(pot, var_order, dom, deg, bnd_method)
  lowest_of_lower = lowest_possible(lower_poly, dom, bnd_method)
  lower_poly = (if lowest_of_lower < -1e-6
                  make_zero(lower_poly)
                else
                  lower_poly
                end
                )
  Patch(dom, lower_poly, upper_poly, Domain[])
end

# make a factor from a potential function
function make_factor(pot :: Potential, cover_dom :: Domain, p :: Partition, var_order, bnd_method="exact", deg=2)
  potential_bounds = Dict{Domain, Patch}()
  for dom in p
    potential_bounds[dom] = patch_pot(var_order, dom, pot, bnd_method, deg)
  end
  Factor(var_order, cover_dom, p, potential_bounds)
end

# grow a factor from a potential function
function pot_grow!(pot :: Potential, f :: Factor, dom :: Domain, bnd_method="exact", deg=2)
  assert(dom in f.partition)
  # shatter the domain to grow on
  shattered_doms = rec_split_half(dom, length(f.var_order))
  # remove the domain from the partition
  setdiff!(f.partition, Domain[dom])
  # add the shattered domain to the partition
  union!(f.partition, shattered_doms)
  # add the patches for the shattered domains
  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_pot(f.var_order, d1, pot, bnd_method, deg)
  end
end

function rand_pot_grow!(pot :: Potential, f :: Factor, bnd_method="exact", deg=2)
  dom_to_grow = find_split_dom(f.cover_dom, f.partition)
  pot_grow!(pot, f, dom_to_grow, bnd_method, deg)
end


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

# give a domain, and 2 factors f1 and f2
function patch_mult (var_order, dom :: Domain, f1 :: Factor, f2 :: Factor)
  shrink1 = diminish_dom_dim(var_order, f1.var_order, dom)
  shrink2 = diminish_dom_dim(var_order, f2.var_order, dom)
  @show(dom, shrink1, shrink2)
  best_cover1 = best_covering(shrink1, keys(f1.potential_bounds))
  best_cover2 = best_covering(shrink2, keys(f2.potential_bounds))
  @show(best_cover1, best_cover2)
  p_l = f1.potential_bounds[best_cover1].p_l * f2.potential_bounds[best_cover2].p_l
  p_u = f1.potential_bounds[best_cover1].p_u * f2.potential_bounds[best_cover2].p_u
  Patch(dom, p_l, p_u, Domain[best_cover1, best_cover2])
end

function f_mult(f1 :: Factor, f2 :: Factor, p :: Partition)
  var_order = merge_varorder(f1.var_order, f2.var_order)
  potential_bounds = Dict{Domain, Patch}()
  for dom in p
    potential_bounds[dom] = patch_mult(var_order, dom, f1, f2)
  end
  Factor(var_order, p, potential_bounds)
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
