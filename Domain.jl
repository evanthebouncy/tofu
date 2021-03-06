# define the basic domain operations
# a domain is simply a list of tuples denoting lower/upper bound
# over a box

typealias Domain Array{(Float64,Float64),1}
typealias Partition Set{Domain}

abstract BSP
type Split <: BSP
  cover_domain :: Domain
  split_var :: ASCIIString
  split_val :: Float64
  left :: BSP
  right :: BSP
end
type Leaf <: BSP
end

#### FUNCTION CONCERNING DOMAINS ============================


# divide a domain in half along the greatest axis
# also gives the index of where the split happened
function split_half(dom :: Domain)
  function split_at_index(dom :: Domain, idx)
    to_split = dom[idx]
    mid_pt = 0.5 * (to_split[1] + to_split[2])
    half1, half2 = (to_split[1], mid_pt), (mid_pt, to_split[2])
    ret1,ret2 = [x for x in dom], [x for x in dom]
    ret1[idx] = half1
    ret2[idx] = half2
    ret1, ret2
  end
  function find_split_idx(dom :: Domain)
    length_idx = [(dom[i][2]-dom[i][1], i) for i in 1:length(dom)]
    max_len, max_idx = max(length_idx..., (0.0, 0))
    max_idx
  end
  idx = find_split_idx(dom)
  ret1, ret2 = split_at_index(dom, idx)
  idx, ret1, ret2
end

function has_intersect(dom1, dom2)
  function interval_intersect(i1, i2)
    a, b = i1
    c, d = i2
    stronger_lower_bnd = max(a, c)
    stronger_upper_bnd = min(b, d)
    stronger_lower_bnd < stronger_upper_bnd
  end
  assert(length(dom1) == length(dom2))
  for i in 1:length(dom1)
    if !(interval_intersect(dom1[i], dom2[i]))
      return false
    end
  end
  true
end

function get_intersect(dom1, dom2)
  function get_interval_intersect(i1, i2)
    a, b = i1
    c, d = i2
    stronger_lower_bnd = max(a, c)
    stronger_upper_bnd = min(b, d)
    (stronger_lower_bnd, stronger_upper_bnd)
  end
  (Float64, Float64)[get_interval_intersect(dom1[i], dom2[i]) for i in 1:length(dom1)]
end

# returns true if dom1 is a subset of dom2
function dom_subset(dom1, dom2)
  get_intersect(dom1, dom2) == dom1
end

# find the center of a domain
function center(dom)
  [0.5*(d[1] + d[2]) for d in dom]
end

# get the largest axis length of the domain
function max_length(dom)
  max([d[2]-d[1] for d in dom]..., 0.0)
end

# get the euclidian diagnal dist
function diag_radius_length(dom)
  sum_of_sq = reduce(+, [(d[2]-d[1])^2 for d in dom])
  sum_of_sq ^ 0.5
end

# give a single random sample of the domain
function get_single_sample(dom1)
  left_end = [x[1] for x in dom1]
  lengthz  = [x[2]-x[1] for x in dom1]
  Float64[left_end[i] + rand() * lengthz[i] for i in 1:length(dom1)]
end

# check if dom contains a point
function dom_contains(dom, pt)
  assert(length(dom) == length(pt))
  for i in 1:length(pt)
    x = pt[i]
    x_bnd_low, x_bnd_high = dom[i]
    if !(x_bnd_low <= x <= x_bnd_high)
      return false
    end
  end
  return true
end

# enlarge a dom to encompass a new dimention
function enlarge_dom_dim(old_var_order, new_var_order, old_dom)
  ret = typeof(old_dom[1])[]
  for v_n in new_var_order
    idx = findfirst(old_var_order, v_n)
    if idx == 0
      push!(ret, (-Inf, Inf))
    else
      push!(ret, old_dom[idx])
    end
  end
  ret
end

# diminish a dom to remove a dimension
function diminish_dom_dim(old_var_order, new_var_order, old_dom)
  ret = typeof(old_dom[1])[]
  for v_n in new_var_order
    idx = findfirst(old_var_order, v_n)
    push!(ret, old_dom[idx])
  end
  ret
end


# FUNCTION CONCERNING PARTITIONS =====================

# p1 and p2 are both partitions of some domain
# give the intersection/overlaps of these 2 partitions
function partition_intersect(p1, p2)
  ret = Set{Domain}()
  for d1 in p1
    for d2 in p2
      if has_intersect(d1, d2)
        union!(ret, get_intersect(d1, d2))
      end
    end
  end
  ret
end

# same as above but we track a set of unflatted domains
function partition_intersect_remember(dic1, dic2)
  ret = Dict{Domain, Set{Domain}}()
  for d1 in keys(dic1)
    for d2 in keys(dic2)
      if has_intersect(d1, d2)
        ret[get_intersect(d1, d2)] = union(dic1[d1], dic2[d2])
      end
    end
  end
  ret
end

function find_split_dom(dom_init :: Domain, p :: Partition)
  rnd_pt = get_single_sample(dom_init)
  for dom in p
    if dom_contains(dom, rnd_pt)
      return dom
    end
  end
  dom_init
end

# generarte a semi-random partition from an initial domain
function gen_test_partition(dom_init::Domain, n_splits)
  ret = Set{Domain}(Domain[dom_init])
  for i in 1:n_splits
    #assert(typeof(ret) == Partition)
    to_split = find_split_dom(dom_init, ret)
    setdiff!(ret, Domain[to_split])
    splt1, splt2 = split_half(to_split)
    union!(ret, Domain[splt1])
    union!(ret, Domain[splt2])
  end
  ret
end

# shatter a partition at a set domain, modify the partition and return the
# shattered doms
function shatter_partition!(dom :: Domain, partition :: Partition)
  assert(dom in partition)
  # shatter the domain to grow on
  shattered_doms = rec_split_half(dom, length(dom))
  # remove the domain from the partition
  setdiff!(partition, Domain[dom])
  # add the shattered domain to the partition
  union!(partition, shattered_doms)
  shattered_doms
end


# check if a value slices a domain
function value_slices(var_order, var_name, value, dom)
  a,b = dom[findfirst(var_order, var_name)]
  a <= value <= b
end

# projection of a partition
function partition_projection(old_var_order, new_var_order, projected_variable, p)
  squashed_layers = Dict{Domain, Set{Domain}}[]

  left_value(dom) = dom[findfirst(old_var_order, projected_variable)][1]
  right_value(dom) = dom[findfirst(old_var_order, projected_variable)][2]

  # a simple sorted values for all the left ends of domains in p
  all_left_values = sort(Float64[x for x in Set([left_value(dom) for dom in p])])
  all_slices = [0.5*(all_left_values[i+1] + all_left_values[i]) for i in 1:length(all_left_values)-1]
  for y in all_slices
    sliced_layer = filter(dom->value_slices(old_var_order, projected_variable, y, dom), p)
    to_add = Dict{Domain, Set{Domain}}()
    for dom in sliced_layer
      add_key = diminish_dom_dim(old_var_order, new_var_order, dom) :: Domain
      add_value = Set{Domain}(Domain[dom])
      to_add[add_key] = add_value
    end
    push!(squashed_layers, to_add)
  end
  reduce((x,y)->partition_intersect_remember(x,y), squashed_layers)
end


###### FUNCTIONS CONCERNING a BSP

function new_leaf(cover_domain :: Domain)
  Split (cover_domain, "", Inf, Leaf(), Leaf())
end

function find_best_containing_domain(bsp :: BSP, x)
  if (typeof(bsp.left) == Leaf) & (dom_contains(bsp.cover_domain, x))
    bsp.cover_domain
  else
    if dom_contains(bsp.left.cover_domain, x)
      find_best_containing_domain(bsp.left, x)
    else
      find_best_containing_domain(bsp.right, x)
    end
  end
end

# find the leaf that contains the domain
function find_leaf(bsp :: BSP, target_domain :: Domain)
  if (typeof(bsp.left) == Leaf)
    bsp
  else
    if dom_subset(target_domain, bsp.left.cover_domain)
      find_leaf(bsp.left, target_domain)
    else
      find_leaf(bsp.right, target_domain)
    end
  end
end

function grow_bsp!(bsp :: BSP, old_dom :: Domain, new_dom1 :: Domain, new_dom2 :: Domain, split_var :: ASCIIString, split_val :: Float64)
  to_grow = find_leaf(bsp, old_dom)
  left = new_leaf(new_dom1)
  right = new_leaf(new_dom2)
  to_grow.split_var = split_var
  to_grow.split_val = split_val
  to_grow.left = left
  to_grow.right = right
end

# for domain of smaller dimension: target_domain
# find the smallest covering of that domain in the bsp, as a lst of domains
function find_smallest_cover(bsp :: BSP, target_domain :: Domain, small_var_order, bsp_var_order, inte_var_name)
  if typeof(bsp.left) == Leaf
    Domain[bsp.cover_domain]
  else
    if bsp.split_var == inte_var_name
      rec_left = find_smallest_cover(bsp.left, target_domain, small_var_order, bsp_var_order, inte_var_name)
      rec_right = find_smallest_cover(bsp.right, target_domain, small_var_order, bsp_var_order, inte_var_name)
      ret = Domain[d for d in rec_left]
      for d in rec_right
        push!(ret, d)
      end
      ret
    else
      left_squished = diminish_dom_dim(bsp_var_order, small_var_order, bsp.left.cover_domain)
      right_squished = diminish_dom_dim(bsp_var_order, small_var_order, bsp.right.cover_domain)
      # if left side contains, go down left
      if dom_subset(target_domain, left_squished)
        find_smallest_cover(bsp.left, target_domain, small_var_order, bsp_var_order, inte_var_name)
      else
        # if right side contains, go down right
        if dom_subset(target_domain, right_squished)
          find_smallest_cover(bsp.right, target_domain, small_var_order, bsp_var_order, inte_var_name)
        # if neither contains, return as it is
        else
          Domain[bsp.cover_domain]
        end
      end
    end
  end
end


# find the most suitable covering of a dom from a set of Doms
function best_covering(dom :: Domain, bsp :: BSP)
  # if base case, it is what it is
  if typeof(bsp.left) == Leaf
    bsp.cover_domain
  else
    # otherwise, try to walk down either left or right, but if neither, return self domain
    if dom_subset(dom, bsp.left.cover_domain)
      best_covering(dom, bsp.left)
    else
      if dom_subset(dom, bsp.right.cover_domain)
        best_covering(dom, bsp.right)
      else
        bsp.cover_domain
      end
    end
  end
end

# find the contamination of the multiplication
# a domain is containminated by the source domain if
# its projection touches the source
function find_mult_containmenate(src_var_order, mult_var_order, src_dom :: Domain, bsp :: BSP)
  if typeof(bsp.left) == Leaf
    Domain[bsp.cover_domain]
  else
    # if the split is not in the src var order, we need to include both end
    if !(bsp.split_var in src_var_order)
      rec_left = find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.left)
      rec_right = find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.right)
      ret = Domain[d for d in rec_left]
      for d in rec_right
        push!(ret, d)
      end
      ret
    else
      left_squished = diminish_dom_dim(mult_var_order, src_var_order, bsp.left.cover_domain)
      if dom_subset(left_squished, src_dom)
        find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.left)
      else
        find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.right)
      end
    end
  end
end

# find all the containminates and track whom contaiminated whom in a dictionary
# src_doms is a subset of the partition, and they cannot overlap, as a result...
# we have the assertion that .. ?
function find_all_mult_containmenates(src_var_order, mult_var_order, src_doms :: Set{Domain}, bsp :: BSP)
  ret = Dict(Domain, Set{Domain})()
  for bad_src in src_doms
    containminated = find_mult_containmenate(src_var_order, mult_var_order, bad_src, bsp)
    for cont_dom in containminated
      if cont_dom in keys(ret)
        union!(ret[cont_dom], bad_src)
      else
        ret[cont_dom] = Set{Domain}(Domain[bad_src])
      end
    end
  end
  ret
end

# find the contamination of integration, same idea as above
function find_inte_containmenate(src_var_order, inte_var_order, inte_var_name, src_dom :: Domain, bsp :: BSP)
  if typeof(bsp.left) == Leaf
    Domain[bsp.cover_domain]
  else
    # if the split is not in the src var order, we need to include both end
    if !(bsp.split_var in src_var_order)
      rec_left = find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.left)
      rec_right = find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.right)
      ret = Domain[d for d in rec_left]
      for d in rec_right
        push!(ret, d)
      end
      ret
    # otherwise, we squish the source domain because that is the one that gets squished in integration
    else
      source_squished = diminish_dom_dim(src_var_order, inte_var_order, src_dom)
      if dom_subset(bsp.left.cover_domain, source_squished)
        find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.left)
      else
        find_mult_containmenate(src_var_order, mult_var_order, src_dom, bsp.right)
      end
    end
  end
end
