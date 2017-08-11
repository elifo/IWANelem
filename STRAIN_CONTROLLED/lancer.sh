mkdir -p outputfiles
gfortran nonlinear3c.f90 water.f90 main.f90 lu_solver.f90 -o huhu
[[ $? -eq 0 ]] && time ./huhu
rm huhu *~ fort.*

