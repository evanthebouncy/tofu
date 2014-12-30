include("/home/evan/Documents/research/tofu/SimplePoly.jl")

feq(x,y) = 1 / (1 + (x - y)^2)
large_lower, large_upper = 1.0, 2.0
small_lower, small_upper = 1.5, 2.0

feq_apx = get_m_projections_approx(feq, ["x","y"], 3, [(large_lower, large_upper),(large_lower, large_upper)])
feq_apx2 = get_m_projections_approx(feq, ["x","y"], 3, [(small_lower, small_upper), (small_lower, small_upper)])
diff_poly = feq_apx2 - feq_apx
solve_max(diff_poly, [(small_lower, small_upper), (small_lower, small_upper)])

using Gadfly
diff_numer(x,y) = peval(diff_poly, [x,y])
plot(x = linspace(small_lower, small_upper, 1000), y = linspace(small_lower, small_upper, 1000), z = diff_numer, Geom.contour)


asdd = Poly1([1,2,3])
asdd2 = Poly1([1,2,3])
asdd == asdd2

f_poly = diff_poly
function diff_enum(x::Vector)
  -1.0 * peval(f_poly, x)
end


linspx = linspace(-2.0, 7.0, 1000)
plot(z=(aa,bb)->diff_enum([aa,bb]), x=linspx, y=linspx, Geom.contour)


partialx = δ(f_poly, "x")
partialy = δ(f_poly, "y")

function diff_enum_grad!(x::Vector, storage::Vector)
    storage[1] = -1.0 * peval(partialx, x)
    storage[2] = -1.0 * peval(partialy, x)
end

function crappy!(x::Vector, storage)
      storage[1] = -1.0 * peval(partialx, x)
    storage[2] = -1.0 * peval(partialy, x)
  diff_enum(x)
end

d4 = DifferentiableFunction(diff_enum,
                            diff_enum_grad!,
                            crappy!)

results = fminbox(d4, [2.001, 0.001], [2.0, 0.0], [3.0, 1.0])

results.f_minimum
res = optimize(diff_enum, diff_enum_grad!,
               [2.0, 0.0],
               method = :gradient_descent,
               grtol = 1e-12,
               iterations = 10,
               store_trace = false,
               show_trace = false)

optimize(diff_enum, diff_enum_grad!, [2.0, 0.0], method = :l_bfgs)


critz = find_all_crit_pt(diff_poly, [(small_lower, small_upper), (small_lower, small_upper)])
pt11 = critz[1]

bounder = get_bound_object(diff_poly)
bounder(pt11)

using Gadfly
partialx = δ(diff_poly, "x")
partialy = δ(diff_poly, "y")

peval(partialx, [3.99479, 3.99479])
peval(partialy, [3.99479, 3.99479])

diff_std = SumPolyProd_to_STDPoly(feq_apx) - SumPolyProd_to_STDPoly(feq_apx2)

peval(diff_std, [4.0, 4.0])
to_print(diff_std)

max_val = 0.0
for xx in linspace(small_lower, small_upper, 1000)
  for yy in linspace(small_lower, small_upper, 1000)
    max_val = max(max_val, peval(diff_std, [xx,yy]))
  end
end
max_val

@show(diff_std)

peval(diff_std, [3.0, 3.44427])

diff_std.terms[2]

function get_estimate_err(poly_idx, var_idx)
  fx1 = feq_apx.polyprods[poly_idx].polys[var_idx]
  fx2 = feq_apx2.polyprods[poly_idx].polys[var_idx]

  fx3 = fx1 - fx2
  ls = linspace(small_lower, small_upper, 1000)
  value = 0.0
  for x in ls
    value = max(abs(peval(fx3, x)), value)
  end
  value
end

[get_estimate_err(i,"x") for i in 1:7]

diff_numer(x,y) = peval(feq_apx,[x,y]) - peval(feq_apx2, [x,y])
diff_numer_apx = get_m_projections_approx(diff_numer, ["x","y"], 17, [(-2.0,2.0),(-2.0,2.0)])

feq_apx_bnd = get_bound_object(feq_apx)
feq_apx_bnd2 = get_bound_object(feq_apx2)
diff_numer_bnd = get_bound_object(diff_numer_apx)

feq_apx_bnd([(1.0,2.0),(1.0,2.0)])
feq_apx_bnd2([(-1.0,0.0),(1.5,2.0)])

get_sound_max_dist(feq_apx_bnd,feq_apx_bnd2,[(small_lower,small_lower),(small_lower,small_upper)],1000)
get_sound_max_dist_spp(feq_apx, feq_apx2, [(-2.0,2.0),(-2.0,2.0)], 1000)
get_sample_max_dist_spp(feq_apx, feq_apx2, [(small_lower,small_lower),(small_lower,small_upper)], 1000)
get_sound_max_dist(diff_numer_bnd, x->(0.0,0.0),[(-2.0,2.0),(-2.0,2.0)],1000)

