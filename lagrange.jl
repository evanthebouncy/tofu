using SymPy

# 1 dimensional polynomial approximation via evaluation at chebyshev points.

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

# legrange interpolation at (x,y) coordinates
function get_legrange_interpolation(samples)
  x = Sym("x")
  all_xj = [sample[1] for sample in samples]
  all_products = reduce(*, [x-xj for xj in all_xj])
  ret = 0
  for (xi, yi) in samples
    ret += yi * ( all_products / (x - xi) ) / ( (all_products / (x - xi)) |> replace(x, xi) |> float )
  end
  ret
end

# give a nth degree polynomial approximation of a function by
# interpolation at chebyshev points over domain [a, b]
function get_polynomial_approx(f, n, a, b)
  denormalizer = get_range_denormalizer(a, b)
  interpolant_pts = map(denormalizer, get_chebyshev_pts(n))
  samples = [(x, f(x)) for x in interpolant_pts]
  println(samples)
  get_legrange_interpolation(samples)
end
