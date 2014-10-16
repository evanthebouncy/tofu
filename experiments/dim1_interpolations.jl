using PolyInterpolate
using Gadfly

function test_fun1(x)
  if x < 0
    1.0
  else
    1.0 / (1 + x)
  end
end

x, test_fun1_approx = get_polynomial_approx(test_fun1, 5, -2, 2)
test_fun1_approx_numer(s) = test_fun1_approx |> replace(x, s) |> float

x, test_fun1_approx_crazy = get_polynomial_approx(test_fun1, 100, -2, 2)
test_fun1_approx_crazy_numer(s) = test_fun1_approx_crazy |> replace(x, s) |> float

layer_real = layer(test_fun1, -2, 2)
layer_approx = layer(test_fun1_approx_numer, -2, 2)
layer_approx_crazy = layer(test_fun1_approx_crazy_numer, -2, 2)
plot(layer_real, layer_approx, layer_approx_crazy)

function test_fun2(x)
  1 / (abs(x - 2) + 1)
end

x, test_fun2_approx = get_polynomial_approx(test_fun2, 8, -8, 10)
test_fun2_approx_numer(s) = test_fun2_approx |> replace(x, s) |> float

plot(test_fun2_approx_numer, -8, 10)

