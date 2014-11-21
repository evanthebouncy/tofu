include("/home/evan/Documents/research/tofu/SimplePoly.jl")
using Gadfly

ls_size = 100
x = linspace(-5, 5, ls_size)
y = linspace(-5, 5, ls_size)

f11(x,y) = 1 / (1 +(x - y)^2)
f22(x,y) = 1 / (1 + (x + y)^2)
f_sq(x,y) = 1 / (1 + (x^2 - y - 2)^2)
f_y(y) = 1 / (1 + (y-2)^2)

f11_approx = get_m_projections_approx(f11, ["x", "y"], 11, [(-5, 5), (-5, 5)])
numer_f11(x,y) = peval(f11_approx, [x,y])
plot(z=numer_f11, x=x, y=y, Geom.contour)
plot(z=f11, x=x, y=y, Geom.contour)
plot(f_y, -5, 5)

f_sq_approx = get_m_projections_approx(f_sq, ["x", "y"], 15, [(-5, 5), (-5, 5)])
f_sq_numer(x,y) = peval(f_sq_approx, [x,y])
plot(z=f_sq_numer, x=x, y=y, Geom.contour)
plot(z=f_sq, x=x, y=y, Geom.contour)

two_const = f11_approx * f_sq_approx
two_const_numer(x,y) = peval(two_const, [x,y])
two_const_numer(1.1, 2.2)
consts_numer(x,y) = peval(two_const, [x,y])

plot(z=consts_numer, x=x, y=y, Geom.contour)

prob_x = ∫(two_const, "x", -5, 5)
prob_x_numer(x) = peval(prob_x, [x])
plot(prob_x_numer, -5, 5)

f_y_approx = get_m_projections_approx(f_y, ["y"], 5, [(-5, 5)])3
f_prod = f_sq_approx * f_y_approx
f_prod_numer(x,y) = peval(f_prod, [x,y])
plot(z=f_prod_numer, x=x, y=y, Geom.contour)

f_x_approx = ∫(f_sq_approx*f_y_approx, "y", -5, 5)
f_x_numer(x) = peval(f_x_approx, [x])
plot(f_x_numer, -5, 5)

f22_approx = get_m_projections_approx(f22, ["x", "y"], 5, [(-5, 5), (-5, 5)])
numer_f22(x,y) = peval(f22_approx, [x,y])
plot(z=numer_f22, x=x, y=y, Geom.contour)

f33_approx = f11_approx * f22_approx
numer_f33(x,y) = peval(f33_approx, [x,y])
plot(z=numer_f33, x=x, y=y, Geom.contour)

f33_int = ∫(f33_approx, "y", -4, 4)
numer_f33_int(x) = peval(f33_int, [x])
numer_f33_int(3.0)
plot(numer_f33_int, -4, 4)


# ========================================testing stuff only...

P = Poly1([1,1])
Q = Poly1([1,1])
R = Poly1([1,1])
S = Poly1([2,1])
T = Poly1([2,1])

Prod1 = PolyProd(2, ["x", "y", "z"], ["x" => P, "y" => Q, "z" => R])
peval(Prod1, [1.0,0.0,0.0])
peval(Prod1, ["x"=>1.0, "y"=>1.0, "z"=>1.0])
Prod2 = PolyProd(3, ["x", "y"], ["x" => S, "y" => T])
Prod3 = ∫(Prod1 * Prod2, "x", 0.0, 1.0)
Prod4 = PolyProd(3, ["x", "y"], ["x" => P, "y" => T])
Prod5 = ∫(Prod1 * Prod4, "x", 4.0, 8.0)

Prod11 = PolyProd(2, ["x", "y", "z"], ["x" => P, "y" => Q, "z" => R])
Prod12 = PolyProd(1, ["x", "y", "z"], ["x" => Q, "y" => P, "z" => R])
spp1 = SumPolyProd(["x", "y", "z"], [Prod11, Prod12])
spp1_sq = spp1 * spp1
spp1_sq.var_order
map(x-> spp1_sq.polyprods
spp2 = ∫(spp1_sq, "x", 3, 5)
spp3 = ∫(spp2, "y", 4, 5)

peval(Prod11, [1.0, 1.0, 0.0])
peval(Prod12, [1.0, 1.0, 0.0])
peval(spp1, [1.0, 1.0, 0.0])
can_join(Prod3, Prod3)
can_join(Prod3, Prod5)
Prod3 + Prod5
pp1 = Poly1([1,2,3])
∫(pp1)
pp4 = Poly1([1,3,4])
pp1 * pp4
pp1 * Poly1([1])
pp1 + pp1 + pp1
pp2 = Poly1([1])
pp3 = pp1 + pp2
peval(pp3, 2)
peval(pp1, 0)
peval(pp1, 4)
ffff(x) = 1 / (1 + x)
using Gadfly

l1 = [pp1]
reduce(+, [l1,l1])
xxx = [1,2,3]
yyy = [convert(Float64,ffff(x)) for x in xxx]

l1 = layer(ffff, 0, 7)
poly1 = get_univar_approximation(ffff, 4, 2, 5)

function func11(x)
  peval(poly1, x)
end
l2 = layer(func11, 0, 7)

ffff(2.5)
func11(2.5)
plot(l1, l2)
get_univar_interpolant([1,2,3],[1,3,2])
A=rand(1:3,3,3)
B = inv(A)
C = A * B

A = reshape(1:9, 3, 3)
B = [1,2,3]
A * B

A = reshape([1,1,1,1,2,3,1,4,9], 3, 3)
A * [1,1,1]
Ainv = inv(A)
Ainv * [1,3,2]

P(x) = -4.0 + 6.5 * x + -1.5 * x^2

P(1)
P(2)
P(3)

map(x->x^2,[1,2,3,4])

typeof([1.1,1.2])

[1,1,1,1] + [1,1,1]

kkl = PolyProd[]

kaka = (x...) -> (reduce(+, x))
kaka(1,2,3)

