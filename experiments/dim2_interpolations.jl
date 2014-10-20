using PolyInterpolate
using Gadfly

size = 100
x = linspace(-2*pi, 2*pi, size)
y = linspace(-2*pi, 2*pi, size)

function g(x,y)
    1 / (1 + (x - y)^2)
end

# approximating it by 2 degree 7 polynomials (total of coefficients)
symvar, symfun = get_m_projections_approx(g, 2, 7, [(-2*pi, 2*pi), (-2*pi, 2*pi)], "x", None)
numerfun = sym_to_numer(symfun, symvar)
# the accurate function
plot(z=g, x=x, y=y, Geom.contour)
# the approximate polynomial of the function
plot(z=numerfun, x=x, y=y, Geom.contour)
