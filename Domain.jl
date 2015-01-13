# define the basic domain operations
# a domain is simply a list of tuples denoting lower/upper bound
# over a box

typealias Domain Array{(Float64,Float64),1}
typealias Partition Array{Domain}

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

# divide a domain in half along the greatest axis
function split_half(dom)
  dom_w_length = [(d[2]-d[1],d) for d in dom]
  size, maxd = max(dom_w_length..., (0.0, dom[1]))
  half1 = (maxd[1], (maxd[1] + maxd[2]) / 2)
  half2 = ((maxd[1] + maxd[2]) / 2, maxd[2])
  ret1 = [x for x in dom]
  ret1[findfirst(ret1, maxd)] = half1
  ret2 = [x for x in dom]
  ret2[findfirst(ret2, maxd)] = half2
  ret1, ret2
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

# p1 and p2 are both partitions of some domain
# give the intersection/overlaps of these 2 partitions
function partition_intersect(p1, p2)
  ret = Domain[]
  for d1 in p1
    for d2 in p2
      if has_intersect(d1, d2)
        push!(ret, get_intersect(d1, d2))
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

# find the center of a domain
function center(dom)
  [0.5*(d[1] + d[2]) for d in dom]
end

# get the largest axis length of the domain
function max_length(dom)
  max([d[2]-d[1] for d in dom]...)
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

# generarte a semi-random partition from an initial domain
function gen_test_partition(dom_init, n_splits)
  ret = Domain[dom_init]
  function find_split_idx()
    rnd_pt = get_single_sample(dom_init)
    for j in 1:length(ret)
      if dom_contains(ret[j], rnd_pt)
        return j
      end
    end
    return 1
  end
  for i in 1:n_splits
    to_split = splice!(ret, find_split_idx())
    splt1, splt2 = split_half(to_split)
    push!(ret, splt1)
    push!(ret, splt2)
  end
  ret
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
