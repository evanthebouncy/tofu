include("/home/evan/Documents/research/tofu/Factor.jl")
using Gadfly

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

# check the factor relationships are correct
length(FG.factors) == 4
FG.rel_factor_child[f_const].f1 == f_const
FG.rel_factor_child[f_const].f2 == f_eq
FG.rel_factor_child[f_const].f_mult == f_mult1
FG.rel_factor_child[f_const] == FG.rel_factor_child[f_eq]
FG.rel_factor_child[f_const] == FG.rel_factor_parents[f_mult1]
FG.rel_factor_child[f_mult1].f1 == f_mult1
FG.rel_factor_child[f_mult1].f_inte == f_inte1
FG.rel_factor_child[f_mult1] == FG.rel_factor_parents[f_inte1]

# check the domain relationships are correct
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][1].d1 == [(0.0, 16.0)]
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][1].d2 == [(0.0, 16.0), (0.0, 16.0)]
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][1].d_mult == [(0.0, 16.0), (0.0, 16.0)]
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][1].fmultop.f1 == f_const
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][1].fmultop.f2 == f_eq
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][1].fmultop.f_mult == f_mult1

FG.rel_domain_children[(f_eq, [(0.0, 16.0),(0.0, 16.0)])]
FG.rel_domain_children[(f_eq, [(0.0, 16.0),(0.0, 16.0)])][1].d2 == [(0.0, 16.0), (0.0, 16.0)]
FG.rel_domain_children[(f_eq, [(0.0, 16.0),(0.0, 16.0)])][1].d_mult == [(0.0, 16.0), (0.0, 16.0)]
FG.rel_domain_children[(f_eq, [(0.0, 16.0),(0.0, 16.0)])][1].fmultop.f1 == f_const
FG.rel_domain_children[(f_eq, [(0.0, 16.0),(0.0, 16.0)])][1].fmultop.f2 == f_eq
FG.rel_domain_children[(f_eq, [(0.0, 16.0),(0.0, 16.0)])][1].fmultop.f_mult == f_mult1

FG.rel_domain_parents[(f_mult1, [(0.0, 16.0),(0.0, 16.0)])].d_mult
FG.rel_domain_children[(f_mult1, [(0.0, 16.0),(0.0, 16.0)])]

# try to go back to see ...should have 4 domains total.
upstreams1 = all_upstream_contributions_rec(f_inte1, f_inte1.partition)
[(f_d.f_name, upstreams1[f_d]) for f_d in keys(upstreams1)]
all_splits = get_valid_split_domains(FG)

rand_f_mult_grow!(FG, f_mult1)
# should have 4 domains still...
upstreams2 = all_upstream_contributions_rec(f_inte1, f_inte1.partition)
[(f_d[1].f_name, f_d[2]) for f_d in upstreams2]
valid_splits = get_valid_split_domains(FG)
[(f_d[1].f_name, f_d[2]) for f_d in valid_splits]
FG.rel_domain_children[(f_const, [(0.0, 16.0)])][5].d_mult
FG.rel_domain_children[(f_eq, [(0.0, 16.0), (0.0, 16.0)])][3].d_mult

# this should be key not foud error
FG.rel_domain_children[(f_mult1, [(0.0, 8.0), (0.0, 8.0)])]
rand_f_inte_grow!(FG, f_inte1)
rand_f_pot_grow!(FG, f_const)
rand_f_pot_grow!(FG, f_eq)
# now it should have 8 domains now valid for splitting
upstreams3 = all_upstream_contributions_rec(f_inte1, f_inte1.partition)
[(f_d.f_name, upstreams1[f_d]) for f_d in keys(upstreams3)]
valid_splits = get_valid_split_domains(FG)
[(f_d[1].f_name, f_d[2]) for f_d in valid_splits]
hahaf, hahadic = propagate_imprecise_lower_rec(FG, f_const, [[(0.0, 16.0)]=>f_const.potential_bounds[[(0.0, 16.0)]].p_l], upstreams3)
hahaf.f_name
keys(hahadic)

hahadic = find_imprecise_error(FG, f_eq, [[(0.0, 16.0), (0.0, 16.0)]=>f_eq.potential_bounds[[(0.0, 16.0), (0.0, 16.0)]].p_l], upstreams3)
shitty = find_all_split_error(FG)
crappy = [(shit_key[1].f_name, shit_key[2], shitty[shit_key]) for shit_key in keys(shitty)]
filter(x->(x[1] == "f_mult1"), crappy)
for i in 1:100
  heuristic_grow!(FG)
end

draw_dom2d(f_mult1.partition)
draw_dom2d(f_eq.partition)

layer1 = layer((x)->feval_lower(f_inte1,[x]), bigdom...)
layer2 = layer((x)->feval_upper(f_inte1,[x]), bigdom...)
plot(layer1, layer2)

