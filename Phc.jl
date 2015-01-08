function write_and_run_phc(lst_polys)
  try
    run(`rm foo.phc`)
    run(`rm foo_ans`)
  catch
  end
  fd = open("foo.phc", "w")
  println(fd, string(length(lst_polys)))
  for poly in lst_polys
    println(fd, string(to_string(poly),";"))
  end
  close(fd)
  run(`/home/evan/Apps/phcpack/phc -b foo.phc foo_ans`)
end

function read_roots(num_var)
  fd_read = open("foo.phc")
  lines = readlines(fd_read)
  sol_idxs = filter(i->contains(lines[i],"the solution for t"), 1:length(lines))
  sols_strs = [[lines[sol_idx + offset] for offset in 1:num_var] for sol_idx in sol_idxs]
  sols_strs
  function get_1_sol(sol_str)
    ret = Dict()
    for sol in sol_str
      splitted = split(sol, r" |\n")
      splitted = filter(x->length(x)>0, splitted)
      try
        rl_part = float64(splitted[3])
        im_part = float64(splitted[4])
        ret[splitted[1]] = (rl_part, im_part)
      catch
      end

    end
    ret
  end
  close(fd_read)
  [get_1_sol(solz) for solz in sols_strs]
end

function solve_system(lst_polys)
  lenz = length(lst_polys)
  write_and_run_phc(lst_polys)
  read_roots(lenz)
end
