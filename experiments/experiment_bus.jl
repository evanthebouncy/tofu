include("/home/evan/Documents/research/tofu/Refinement.jl")

bigdom = (0.0, 60.0) :: (Float64, Float64)
dom_x = linspace(bigdom..., 1000)
dom_y = linspace(bigdom..., 1000)

function initialize_test ()
  bigdom = (0.0, 60.0) :: (Float64, Float64)
  dom_x = linspace(bigdom..., 1000)
  dom_y = linspace(bigdom..., 1000)

  const_pot = get_potential_from_dist(get_const_dist(10.0))
  ob_pot = get_potential_from_dist(get_const_dist(20.0))
  plus_pot = get_potential_from_dist(plus_dist)

  FG = init_factor_graph()

  # our constant c is 20
  f_const = f_pot(FG, "f_c", const_pot, ["c"], (Float64,Float64)[bigdom for i in 1:1], "approx", 1)
  # our observation is 10
  f_ob = f_pot(FG, "f_ob", ob_pot, ["ob"], (Float64,Float64)[bigdom for i in 1:1], "approx", 1)
  # end = start + 20
  f_plus1 = f_pot(FG, "f_plus", plus_pot, ["end", "start", "c"], (Float64,Float64)[bigdom for i in 1:3], "approx", 1)
  # the observation is drawn from uniform
  f_unif1 = f_pot(FG, "f_unif1", uniform_pot, ["ob", "start", "end"],(Float64,Float64)[bigdom for i in 1:3],  "approx", 1)
  # the prediction is also drawn from uniform
  f_unif2 = f_pot(FG, "f_unif2", uniform_pot, ["pred", "start", "end"], (Float64,Float64)[bigdom for i in 1:3],  "approx", 1)

  # we first multiply the plus by const
  f_plus2 = f_mult(FG, "f_plus2", f_plus1, f_const, (Float64,Float64)[bigdom for i in 1:3])
  # then rid of c...
  f_start_end = f_inte(FG, "f_start_end", f_plus2, "c", (Float64,Float64)[bigdom for i in 1:2])
  # we then multiply start_end with observation
  f_join_obs = f_mult(FG, "f_join_obs", f_start_end, f_unif1, (Float64,Float64)[bigdom for i in 1:3])
  f_join_obs2 = f_mult(FG, "f_join_obs2", f_join_obs, f_ob, (Float64,Float64)[bigdom for i in 1:3])

  # we then integrate away the constant 10
  f_unif_update = f_inte(FG, "f_unif_update", f_join_obs2, "ob", (Float64,Float64)[bigdom for i in 1:2])
  # we then multiply by our final uniform
  f_unif_pred = f_mult(FG, "f_unif_pred", f_unif_update, f_unif2, (Float64,Float64)[bigdom for i in 1:3])
  # integrate away the 2 unused
  f_unif_nostart = f_inte(FG, "f_unif_nostart", f_unif_pred, "start", (Float64,Float64)[bigdom for i in 1:2])
  f_answer = f_inte(FG, "f_answer", f_unif_nostart, "end", (Float64,Float64)[bigdom for i in 1:1])
  FG
end

FG = initialize_test()

FG.factors[9].potential_bounds

function profile_test()
  for i in 1:1000
    println(i)
    heuristic_grow!(FG)
  end
end

profile_test()

yada = give_sorted_partition_imprecision(FG)

displayshit(Collections.dequeue!(yada))

# Profile.init(10^8, 0.1)
# Profile.clear()
# @profile profile_test()
# using ProfileView
# ProfileView.view()

using Gadfly

layer1 = layer((x)->feval_lower(last(FG.factors),[x]), bigdom...)
layer2 = layer((x)->feval_upper(last(FG.factors),[x]), bigdom...)
plot(layer1, layer2)


draw_dom1d(last(FG.factors).partition)

for f in FG.factors
  println(f.f_name, "    ", length(f.partition))
end

