# The Durand_Kerner method for finding root of a uni-variate polynomial
# The leading coefficient for p (the one with the highest degree) must be 1
function Durand_Kerner(p, deg)
  α = 0.9 + 0.4im
  vect = [α^i for i in 1:deg]
  for i in 1:30
    vect_new = []
    for j in 1:deg
      rootj = vect[j]
      others = [vect[i] for i in filter(x->x!=j, 1:deg)]
      rootj_new = rootj - p(rootj) / reduce(*, [rootj - rooti for rooti in others])
      vect_new = [vect_new, rootj_new]
    end
    vect = vect_new
  end
  vect
end

test_fun(x) = (5*x^6 + 8x^3 - 7) / 5

answer = Durand_Kerner(test_fun, 6)
test_fun(answer[1])
test_fun(answer[2])
test_fun(answer[3])
test_fun(answer[4])


reduce(*,([1,[2,3]]))

[i+i for i in 1:4[(x) -> (x != 2)]]

