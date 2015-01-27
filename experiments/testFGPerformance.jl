using ProfileView
include("/home/evan/Documents/research/tofu/Factor.jl")

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

heuristic_grow!(FG)

function profile_test()
  for i in 1:20
    heuristic_grow!(FG)
  end
end

Profile.init(10^8, 0.001)
Profile.clear()
@profile profile_test()
ProfileView.view()



# sanity check to see if we ruined it by caching
using Gadfly

draw_dom2d(f_mult1.partition)
draw_dom2d(f_eq.partition)
layer1 = layer((x)->feval_lower(f_inte1,[x]), bigdom...)
layer2 = layer((x)->feval_upper(f_inte1,[x]), bigdom...)
plot(layer1, layer2)

