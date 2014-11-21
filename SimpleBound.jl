# divide a domain in half along the greatest axis
function halves(dom)
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

# given 2 bound objects f_bnd_obj, g_bnd_obj
# each of them able of giving lower/upper bound over a domain
# for function f and g respectively
# given also a domain dom
# compute a, b
# where a >= max(f - g) over the domain
# where b >= max(g - f) over the domain
function get_max_dist(f_bnd_obj, g_bnd_obj, dom)
  f_min, f_max = f_bnd_obj(dom)
  g_min, g_max = g_bnd_obj(dom)
  f_max - g_min, g_max - f_min
end

# recursively apply get_max_dist on smaller and smaller domains
# returns A where A >= max(f - g) over the domain
function get_sound_max_dist(f_bnd_obj, g_bnd_obj, init_domain, rep=1000)
  pq = Collections.PriorityQueue()
  init_cost = -1.0 * get_max_dist(f_bnd_obj, g_bnd_obj, init_domain)[1]
  Collections.enqueue!(pq, (init_domain, init_cost), init_cost)
  for i in 1:rep
    dom, cost = Collections.dequeue!(pq)
    if i % 50 == 0
      @show(dom, cost)
    end
    half1, half2 = halves(dom)
    cost1 = -1.0 * get_max_dist(f_bnd_obj, g_bnd_obj, half1)[1]
    cost2 = -1.0 * get_max_dist(f_bnd_obj, g_bnd_obj, half2)[1]
    Collections.enqueue!(pq, (half1,cost1), cost1)
    Collections.enqueue!(pq, (half2,cost2), cost2)
  end
  final_dom, final_cost = Collections.dequeue!(pq)
  get_max_dist(f_bnd_obj, g_bnd_obj, final_dom)[1]
end

function get_sound_max_dist_spp(f :: SumPolyProd, g :: SumPolyProd, init_domain, rep=1000)
  diff_func = f - g
  diff_bnd = get_bound_object(diff_func)
  zero_bnd(x) = (0.0,0.0)
  get_sound_max_dist(diff_bnd, zero_bnd, init_domain, rep)
end

function get_sample_max_dist_spp(f::SumPolyProd, g::SumPolyProd, domain, n_sample = 1000)
  left_end = [x[1] for x in domain]
  lengthz  = [x[2]-x[1] for x in domain]
  function get_single_sample()
    [left_end[i] + rand() * lengthz[i] for i in 1:length(domain)]
  end
  max_location = left_end
  max_val = -Inf
  for i in 1:n_sample
    fresh_sample = get_single_sample()
    cur_val = peval(f, fresh_sample) - peval(g, fresh_sample)
    if cur_val > max_val
      max_location = fresh_sample
      max_val = cur_val
    end
  end
  max_val, max_location
end
