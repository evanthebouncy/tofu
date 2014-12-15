include("/home/evan/Documents/research/tofu/Abstraction.jl")
dom_small = [(-5.0, 5.0), (-5.0, 5.0)]

equal_potential = get_potential_from_dist(equal_dist)
poly_abstr = get_poly_abstraction(equal_potential, ["x", "y"], [(d_low, d_high), (d_low, d_high)], 3)

using Gadfly
plot(z= (asd...) -> peval(poly_abstr[1], [asd...]), x=linspace(1.0, 2.0, 1000), y=linspace(1.0, 2.0, 1000), Geom.contour)

plot(z= (asd...) -> peval(poly_abstr[2], [asd...]), x=linspace(1.0, 2.0, 1000), y=linspace(1.0, 2.0, 1000), Geom.contour)


best_thing = -Inf
for xx in linspace(1.0, 2.0, 1000)
  for yy in linspace(1.0, 2.0, 1000)
    best_thing = max(diff_thing(1.5, 1.5), best_thing)
  end
end
best_thing

diff_thing(1.0287537825579567,1.062499999999964)
diff_thing(1.124999999999801,1.2103772041450853)

