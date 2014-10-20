using SymPy
using Iterators

###########################################################################
# 1 dimensional polynomial approximation via evaluation at chebyshev points.
###########################################################################

# contract [a,b] to the range of [1-,1]
function get_range_normalizer(a,b)
  range_normalizer(x) = 2 / (b - a) * (x - a) - 1
  range_normalizer
end

# extend [-1, 1] to the range of [a, b]
function get_range_denormalizer(a,b)
  range_denomalizer(x) = (b - a) / 2 * (x + 1) + a
  range_denomalizer
end

# contract the function from [a, b] to the [-1, 1] domain
function get_normalized_function(f, a, b)
  denormalizer = get_range_denormalizer(a,b)
  normalized_function(t) = f(denormalizer(t))
  normalized_function
end

# extend the function from the [-1, 1] domain to [a, b]
function get_denormalized_function(g, a, b)
  normalizer = get_range_normalizer(a, b)
  denormalized_function(x) = g(normalizer(x))
  denormalized_function
end

# chebyshev points
function get_chebyshev_pts(n)
  [cos(j*pi/n) for j in 0:n]
end

# legrange interpolation at (x,f(x)) samples
function get_legrange_interpolation(samples, name)
  x = Sym(name)
  all_xj = [sample[1] for sample in samples]
  all_products = reduce(*, [x-xj for xj in all_xj])
  ret = 0
  for (xi, yi) in samples
    ret += yi * ( all_products / (x - xi) ) / ( (all_products / (x - xi)) |> replace(x, xi) |> float )
  end
  x, ret
end

# give a nth degree polynomial approximation of a function by
# interpolation at chebyshev points over domain [a, b]
function get_polynomial_approx(f, n, a, b, name="x")
  denormalizer = get_range_denormalizer(a, b)
  interpolant_pts = map(denormalizer, get_chebyshev_pts(n))
  samples = [(x, f(x)) for x in interpolant_pts]
  get_legrange_interpolation(samples, name)
end


################################################################################
# multi-dimensional polynomial approximation by projection onto the outerproduct
################################################################################

# get a d dimensional tchebyshev grid of n sample in each dimension
# over the domain of the box
function get_chebyshev_grid(d, n, box_domain)
  scaled_nodes = [map(get_range_denormalizer(box_domain[i]...), get_chebyshev_pts(n)) for i in 1:d]
  product(scaled_nodes...)
end

# get the greatest error of a function over a grid
function sample_greatest_err(f, grid)
  max([(f(pt...),pt) for pt in grid]...)
end

# get the function on each lines radiating from a point evaluated
function get_univar_func(f, pt)
  ret = []
  for idx in 1:length(pt)
    ret = [ret, (x) -> (
            eval_pt = map(j->(j == idx ? x : pt[j]), 1:length(pt));
            f(eval_pt...)
          )]
  end
  return ret
end

# get the approximation by pivot at the greatest error points
function get_single_approx(f, d, n, box_domain, name="x")
  # set up the grid to do the sampling
  grid = get_chebyshev_grid(d, n, box_domain)
  # find by sampling the biggest error
  err, pivot = sample_greatest_err(f, grid)
  # construct univariate functions around the pivot
  univar_funs = get_univar_func(f, pivot)

  univar_funs_approx = []
  for i in 1:d
    dom_lower, dom_upper = box_domain[i];
    univar_funs_approx = [univar_funs_approx,
                          get_polynomial_approx(univar_funs[i], n, dom_lower, dom_upper, string(name,string(i)) )
                         ]
  end
  sym_vars = map(x->x[1], univar_funs_approx)
  sym_expr = reduce(*, map(x->x[2], univar_funs_approx))
  sym_expr = simplify(sym_expr) / (err ^ (d-1))
  sym_vars, sym_expr
end

# turn a symbolic function to an executable one
function sym_to_numer(f, sym_vars)
  function cur_approx_numer(x...)
    result = f
    for j in 1:length(sym_vars)
      result = result |> replace(sym_vars[j], x[j])
    end;
    result |> float
  end
  cur_approx_numer
end

# get n projections onto the outer product space
function get_m_projections_approx(f, d, n, box_domain, name, rep=None)

  if rep == None
    rep = n
  end

  function get_1_projection_chain(err_f, sym_f, sym_vars)
    sym_vars, cur_approx = get_single_approx(err_f,d,n,box_domain,name)
    sym_f_new = sym_f + cur_approx
    cur_approx_numer = sym_to_numer(cur_approx, sym_vars)
    err_f_new(x...) = err_f(x...) - cur_approx_numer(x...)
    err_f_new, sym_f_new, sym_vars
  end

  function m_projection_rec(m)
    if m == 0
      f, Sym(0), Any
    else
      get_1_projection_chain(m_projection_rec(m-1)...)
    end
  end
  last_err_f, last_sym_f, sym_vars = m_projection_rec(rep)
  sym_vars, last_sym_f
end
