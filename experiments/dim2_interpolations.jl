using PolyInterpolate
using Gadfly
using SymPy
using Calculus

size = 100
x = linspace(-2*pi, 2*pi, size)
y = linspace(-2*pi, 2*pi, size)

function g(x,y)
    1 / (1 + (x - y)^2)
end

function g22(x,y)
  if x < 0
    x + abs(y)
  else
    abs(x - y)
  end
end


integrate
integrate(x -> -sin(x), 0.0, float64(pi))

function g22_integrated(y)
  integrate(x -> g(y,x), float64(-2*pi), float64(2*pi))
end

g22_integrated(0.4)
plot(g22_integrated, -2*pi, 2*pi)

# approximating it by 2 degree 7 polynomials (total of coefficients)|
symvar, symfun = get_m_projections_approx(g, 2, 7, [(-2*pi, 2*pi), (-2*pi, 2*pi)], "x", None)
numerfun = sym_to_numer(symfun, symvar)

symx1 = Sym("x1")
asdf = integrate(symfun, Sym("x1"))
definite_asdf = (asdf |> replace(symx1, 2*pi)) - (asdf |> replace(symx1, -2*pi))

numerfun_asdf = (x) -> (definite_asdf |> replace(Sym("x2"), x) |> float)

definite_asdf |> (Sym("x2"), 2)

numerfun_asdf(2)

plot(numerfun_asdf, -2*pi, 2*pi)

# the accurate function
plot(z=g, x=x, y=y, Geom.contour)
# the approximate polynomial of the function
plot(z=numerfun, x=x, y=y, Geom.contour)

