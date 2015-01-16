include("/home/evan/Documents/research/tofu/Abstraction.jl")
type Patch
  dom :: Domain
  p_l :: SumPolyProdC
  p_u :: SumPolyProdC
  parents :: Array{Domain}
end

type Factor
  var_order
  bsp :: BSP
  partition :: Partition
  potential_bounds :: Dict{Domain, Patch}
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
function f_pot(pot :: Potential, var_order, dom :: Domain, bnd_method="exact", deg=2)
  bsp = new_leaf(dom)
  partition = Set{Domain}(Domain[dom])
  potential_bounds = Dict{Domain, Patch}()
  potential_bounds[dom] = patch_pot(var_order, dom, pot, bnd_method, deg)
  Factor(var_order, bsp, partition, potential_bounds)
end

# grow a factor from a potential function
function f_pot_grow!(f :: Factor, dom :: Domain, pot :: Potential, bnd_method="exact", deg=2)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  # add the patches for the shattered domains
  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_pot(f.var_order, d1, pot, bnd_method, deg)
  end
end

function rand_f_pot_grow!(f :: Factor, pot :: Potential, bnd_method="exact", deg=2)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_pot_grow!(f, dom_to_grow, pot, bnd_method, deg)
end

# give a domain, and 2 factors f1 and f2
function patch_mult (var_order, dom :: Domain, f1 :: Factor, f2 :: Factor)
  shrink1 = diminish_dom_dim(var_order, f1.var_order, dom)
  shrink2 = diminish_dom_dim(var_order, f2.var_order, dom)
  best_cover1 = best_covering(shrink1, keys(f1.potential_bounds))
  best_cover2 = best_covering(shrink2, keys(f2.potential_bounds))
  p_l = f1.potential_bounds[best_cover1].p_l * f2.potential_bounds[best_cover2].p_l
  p_u = f1.potential_bounds[best_cover1].p_u * f2.potential_bounds[best_cover2].p_u
  Patch(dom, p_l, p_u, Domain[best_cover1, best_cover2])
end

function f_mult(f1 :: Factor, f2 :: Factor, cover_dom :: Domain)
  var_order = merge_varorder(f1.var_order, f2.var_order)
  bsp = new_leaf(cover_dom)
  partition = Set{Domain}(Domain[cover_dom])
  potential_bounds = Dict{Domain, Patch}()
  potential_bounds[cover_dom] = patch_mult(var_order, cover_dom, f1, f2)
  Factor(var_order, bsp, partition, potential_bounds)
end

# grow a factor from a multiplication
function f_mult_grow!(f :: Factor, dom :: Domain, f1 :: Factor, f2 :: Factor)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_mult(f.var_order, d1, f1, f2)
  end
end

function rand_f_mult_grow!(f :: Factor, f1 :: Factor, f2 :: Factor)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_mult_grow!(f, dom_to_grow, f1, f2)
end


# give a domain, and the factor it integrates from
function patch_integrate (var_order, dom :: Domain, f1 :: Factor, inte_var_name)
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
  Patch(dom, lower_pol, upper_pol, smallest_cover)
end

function f_inte(f1 :: Factor, inte_var_name :: ASCIIString, cover_dom :: Domain)
  var_order = filter(x->x!=inte_var_name, f1.var_order)
  bsp = new_leaf(cover_dom)
  partition = Set{Domain}(Domain[cover_dom])
  potential_bounds = Dict{Domain, Patch}()
  potential_bounds[cover_dom] = patch_integrate(var_order, cover_dom, f1, inte_var_name)
  Factor(var_order, bsp, partition, potential_bounds)
end

# grow a factor from a potential function
function f_inte_grow!(f :: Factor, dom :: Domain, f1 :: Factor, inte_var_name :: ASCIIString)
  shattered_doms = shatter!(f.var_order, dom, f.partition, f.bsp)

  for d1 in shattered_doms
    f.potential_bounds[d1] = patch_integrate(f.var_order, d1, f1, inte_var_name)
  end
end

function rand_f_inte_grow!(f :: Factor, f1 :: Factor, inte_var_name :: ASCIIString)
  dom_to_grow = find_split_dom(f.bsp.cover_domain, f.partition)
  f_inte_grow!(f, dom_to_grow, f1, inte_var_name)
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
