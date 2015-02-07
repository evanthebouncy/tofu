include("/home/evan/Documents/research/tofu/Refinement.jl")

bigdom = (0.0, 16.0) :: (Float64, Float64)
dom_x = linspace(bigdom..., 1000)
dom_y = linspace(bigdom..., 1000)

const_pot = get_potential_from_dist(get_const_dist(5.0))
eq_pot = get_potential_from_dist(equal_dist)

FG = init_factor_graph()
f_const = f_pot(FG, "f_c", const_pot, ["x"], (Float64,Float64)[bigdom for i in 1:1], "approx", 2)
f_eq = f_pot(FG, "f_eq", eq_pot,["y", "x"],(Float64,Float64)[bigdom for i in 1:2],  "approx", 2)
f_mult1 = f_mult(FG, "f_mult1", f_const, f_eq, (Float64,Float64)[bigdom for i in 1:2])
f_inte1 = f_inte(FG, "f_inte1", f_mult1, "x", (Float64,Float64)[bigdom for i in 1:1])

for i in 1:1
  @show(i)
  heuristic_grow!(FG)
end

using Gadfly

draw_dom1d(f_const.partition)
draw_f_imprecision1d(FG, f_const)

draw_dom2d(f_eq.partition)

draw_f_imprecision2d(FG, f_eq)

draw_dom2d(f_mult1.partition)

draw_f_imprecision2d(FG, f_mult1)

draw_dom1d(f_inte1.partition)

draw_f_imprecision1d(FG, f_inte1)


layer1 = layer((x)->feval_lower(f_inte1,[x]), bigdom...)
layer2 = layer((x)->feval_upper(f_inte1,[x]), bigdom...)
plot(layer1, layer2)

