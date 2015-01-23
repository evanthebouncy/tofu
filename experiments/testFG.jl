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
[(f_d[1].f_name, f_d[2]) for f_d in upstreams1]

rand_f_mult_grow!(FG, f_mult1, f_const, f_eq)
# should have 4 domains still...
upstreams2 = all_upstream_contributions_rec(f_inte1, f_inte1.partition)
[(f_d[1].f_name, f_d[2]) for f_d in upstreams2]

FG.rel_domain_children[(f_const, [(0.0, 16.0)])][5].d_mult
FG.rel_domain_children[(f_eq, [(0.0, 16.0), (0.0, 16.0)])][3].d_mult

# this should be key not foud error
FG.rel_domain_children[(f_mult1, [(0.0, 8.0), (0.0, 8.0)])]
rand_f_inte_grow!(FG, f_inte1, f_mult1, "x")
# now it should have 8 domains now...
upstreams3 = all_upstream_contributions_rec(f_inte1, f_inte1.partition)
[(f_d[1].f_name, f_d[2]) for f_d in upstreams3]

valid_splits = get_valid_split_domains(FG)
[(f_d[1].f_name, f_d[2]) for f_d in valid_splits]
rand_f_mult_grow!(FG, f_mult1, f_const, f_eq)

valid_splits2 = get_valid_split_domains(FG)
[(f_d[1].f_name, f_d[2]) for f_d in valid_splits2]
rand_f_inte_grow!(FG, f_inte1, f_mult1, "x")
valid_splits2 = get_valid_split_domains(FG)
[(f_d[1].f_name, f_d[2]) for f_d in valid_splits2]

# this should now have a child
FG.rel_domain_children[(f_mult1, [(0.0, 8.0), (0.0, 8.0)])]

# this should also have 2 parent domains
FG.rel_domain_parents[(f_inte1, [(0.0, 8.0)])].d1s
# this should have 2 parent domains as well ...
FG.rel_domain_parents[(f_inte1, [(8.0, 16.0)])].d1s

