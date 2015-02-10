include("/home/evan/Documents/research/tofu/Domain.jl")

# ===========================================================================
# the distance functions of the potentials
# these are the basic constructs of the inference framework
# all distances are measured in euclidian distances

# constant x == c for some c
function get_const_dist(c)
  function const_dist(x)
    abs(x - c)
  end
end

# equality x1 == x2
function equal_dist(x1, x2)
  (0.5 * (x1 - x2)^2) ^ 0.5
end

# negation x1 == -x2
function neg_dist(x1, x2)
  (0.5 * (x1 + x2)^2) ^ 0.5
end

# less or equalto x1 <= x2
function le_dist(x1, x2)
  y = 1.0
  function lt_case1(y, x1, x2)
    if (x1 >= x2)
      y^2
    else
      y^2 + 2 * ((x1-x2)/2)^2
    end
  end
  function lt_case2(y, x1, x2)
    if (x1 <= x2)
      (y-1)^2
    else
      (y-1)^2 + 2 * ((x1-x2)/2)^2
    end
  end
  (min(lt_case1(y,x1,x2),lt_case2(y,x1,x2))) ^ 0.5
end

# addition y == x1 + x2
function plus_dist(y, x1, x2)
  (3 * (1/3 * (x1+x2-y))^2) ^ 0.5
end



# ======================================================================
# convert distances to actual potential functions!

# convert distance into potential
function dist_to_potential_sharp(dist)
  delta = 1.5
  delta_01 = 1 / (1 + 2 ^ (-1*delta))
  amped_err = (1-delta_01 + delta_01*dist) ^ delta
  1.0 / amped_err
end

# a gaussian version of it
function dist_to_potential(dist)
  2.0^(-1.0 * dist^2) + 0.001
end

# convert distance into potential bounder
function get_bounder_from_dist(dist_fun)
  function bounder(dom)
    mid_pt = center(dom)
    radius = diag_radius_length(dom)
    min_dist = max(0.0, dist_fun(mid_pt...) - radius)
    max_dist = dist_fun(mid_pt...) + radius
    pot_close, pot_far = dist_to_potential(min_dist), dist_to_potential(max_dist)
    # the lower bound actually comes from the far point when converted to potential
    pot_far, pot_close
  end
  bounder
end

# the potential object consist of a potential function
# and a potential bounder function
immutable Potential
  potential_fun :: Function
  potential_bnd :: Function
end

function get_potential_from_dist(dist)
  pot_fun(x...) = dist_to_potential(dist(x...))
  pot_bnd_fun = get_bounder_from_dist(dist)
  Potential(pot_fun, pot_bnd_fun)
end

# uniform potential
uniform_pot = Potential(
  ((x, a, b) ->
   (if (a < x < b) & ((b - a) > 1e-4)
      1.0 / (b - a)
    else
      0.0
    end
    )), (x, a, b) -> 0.0)

