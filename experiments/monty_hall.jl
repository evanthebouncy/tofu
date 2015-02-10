include("/home/evan/Documents/research/tofu/Refinement.jl")

bigdom = (-1.0, 3.0) :: (Float64, Float64)
smalldom = (-1.0, 2.0)
dom_x = linspace(bigdom..., 1000)
dom_y = linspace(bigdom..., 1000)

p_door_pot = get_discrete_uniform([0.0, 1.0, 2.0])
neq_pot = get_potential_from_dist(unequal_dist)
g_door_pot = get_discrete_uniform([0.0])
eq_eq_pot = get_potential_from_dist(eq_eq_dist)

FG = init_factor_graph()

f_p_door = f_pot(FG, "pdoor", p_door_pot, ["p"], (Float64,Float64)[bigdom for i in 1:1], "approx", 2)
f_g_door = f_pot(FG, "gdoor", g_door_pot, ["g"], (Float64,Float64)[bigdom for i in 1:1], "approx", 2)
f_neq1 = f_pot(FG, "neq1", neq_pot, ["p", "h"], (Float64,Float64)[bigdom for i in 1:2], "approx", 2)
f_neq2 = f_pot(FG, "neq2", neq_pot, ["g", "h"], (Float64,Float64)[bigdom for i in 1:2], "approx", 2)
f_eqeq = f_pot(FG, "eqeq", eq_eq_pot, ["ans","p","g"], (Float64,Float64)[bigdom for i in 1:3], "approx", 3)

mult1 = f_mult(FG, "mult1", f_neq1, f_neq2, (Float64,Float64)[bigdom for i in 1:3])
inte1 = f_inte(FG, "inte1", mult1, "h", (Float64,Float64)[bigdom for i in 1:2])
mult2 = f_mult(FG, "mult2", f_p_door, f_g_door, (Float64,Float64)[bigdom for i in 1:2])
mult3 = f_mult(FG, "mult3", mult2, inte1, (Float64,Float64)[bigdom for i in 1:2])
mult4 = f_mult(FG, "mult4", mult3, f_eqeq, (Float64,Float64)[bigdom for i in 1:3])
inte2 = f_inte(FG, "inte2", mult4, "p", (Float64,Float64)[bigdom for i in 1:2])
inte3 = f_inte(FG, "inte3", inte2, "g", (Float64,Float64)[bigdom for i in 1:1])

for i in 1:200
  @show(i)
  heuristic_grow!(FG)
end

layer1 = layer((x)->feval_lower(inte3,[x]), bigdom...)
layer2 = layer((x)->feval_upper(inte3,[x]), bigdom...)
plot(layer1, layer2)


unequal_dist(2.0, 0.0)
plot(z=(x,y)->neq_pot1.potential_fun(x,y), x=dom_x, y=dom_y, Geom.contour)

layer1 = layer((x)->p_door_pot.potential_fun(x), bigdom...)
plot(layer1)

