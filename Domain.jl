# define the basic domain operations
# a domain is simply a list of tuples denoting lower/upper bound
# over a box

typealias Domain Array{(Float64,Float64),1}

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


# function intersect_var(d1 :: DicDomain, d2 :: DicDomain)
#   var1, var2 = keys(d1), keys(d2)
#   filter(x->x in var2, var1)
# end

# function has_intersect(d1 :: DicDomain, d2 :: DicDomain)
#   length(intersect_var(d1, d2)) > 0
# end

# # check if dom contains a point
# function dom_contains(dom :: DicDomain, pt)
#   assert(length(dom) == length(pt))
#   for i in 1:length(pt)
#     x = pt[i]
#     x_bnd_low, x_bnd_high = dom[i]
#     if !(x_bnd_low <= x <= x_bnd_high)
#       return false
#     end
#   end
#   return true
# end
