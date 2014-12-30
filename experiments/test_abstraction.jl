using Gadfly
include("/home/evan/Documents/research/tofu/Abstraction.jl")

function drawStuff(doms)
  doms = [d for d in doms]
  color = [rand() * 1000 for i in doms]
  x_min = Float64[x[1][1] for x in doms]
  x_max = Float64[x[1][2] for x in doms]
  y_min = Float64[x[2][1] for x in doms]
  y_max = Float64[x[2][2] for x in doms]
  plot(x_max = x_max, x_min = x_min, y_max = y_max, y_min = y_min,
       color=color,
       Geom.rectbin)
end

dom_small = [(-1.0, 0), (0, 1.0)]
dom_large = [(-2.0, 0), (0, 2.0)]
dom_x = linspace(dom_small[1][1], dom_small[1][2], 1000)
dom_y = linspace(dom_small[2][1], dom_small[2][2], 1000)

neg_potential = get_potential_from_dist(neg_dist)
poly_apx = get_m_projections_approx(neg_potential.potential_fun, ["x", "y"], 5, dom_large)
poly_apx2 = get_m_projections_approx(neg_potential.potential_fun, ["x", "y"], 5, dom_small)

poly_diff = poly_apx - poly_apx2
all_partials = SumPolyProd[δ(poly_diff, x) for x in poly_diff.var_order]
all_partials_bnd_obj = [get_bound_object(spp1) for spp1 in all_partials]
diff_bound_poly = get_bound_object(poly_diff)
asdf1, bcd1 = get_max_value(diff_bound_poly, dom_small)
asdf2, bcd2 = get_max_value(diff_bound_poly, dom_small, all_partials_bnd_obj)
drawStuff(keys(bcd1))

asdf
poly_bnd = get_bound_object(poly_apx)
poly_bnd(dom_small)
diff_fun(x,y) = peval(poly_apx, [x,y]) - neg_potential.potential_fun(x,y)

function diff_bound(dom)
  @show(dom)
  aa,bb = neg_potential.potential_bnd(dom)
  a, b = poly_bnd(dom)
  (b - aa, bb - a)
end


plot(z=diff_fun, x = dom_x, y = dom_y, Geom.contour)

asdf, bcd = get_max_value(diff_bound, dom_small)
bcd = keys(bcd)
drawStuff(bcd)
poly_abstr_neg = get_poly_abstraction(neg_potential, ["x", "y"], dom_small, 5)
poly_abstr_neg[1].c, poly_abstr_neg[2].c

plot(z= (asd...) -> peval(poly_abstr_neg[1], [asd...]), x = dom_x, y = dom_y, Geom.contour)
plot(z= (asd...) -> peval(poly_abstr_neg[2], [asd...]), x = dom_x, y = dom_y, Geom.contour)

le_potential = get_potential_from_dist(le_dist)
poly_abstr_le = get_poly_abstraction(le_potential, ["x", "y"], dom_small, 3)
lowerz = poly_abstr_le[1]
upperz = poly_abstr_le[2]

plot(z= (asd...) -> peval(poly_abstr_le[1], [asd...]), x = dom_x, y = dom_y, Geom.contour)
plot(z= (asd...) -> peval(poly_abstr_le[2], [asd...]), x = dom_x, y = dom_y, Geom.contour)

mult_lower = SPPCMult([poly_abstr_neg[1], poly_abstr_le[1]])
mult_upper = SPPCMult([poly_abstr_neg[2], poly_abstr_le[2]])

plot(z= (asd...) -> peval(mult_lower, [asd...]), x = dom_x, y = dom_y, Geom.contour)
plot(z= (asd...) -> peval(mult_upper, [asd...]), x = dom_x, y = dom_y, Geom.contour)

best_thing = -Inf
for xx in linspace(1.0, 2.0, 1000)
  for yy in linspace(1.0, 2.0, 1000)
    best_thing = max(diff_thing(1.5, 1.5), best_thing)
  end
end
best_thing

diff_thing(1.0287537825579567,1.062499999999964)
diff_thing(1.124999999999801,1.2103772041450853)

plot(z= (asd...) -> peval(poly_abstr[1], [asd...]), x = dom_x, y = dom_y, Geom.contour)
