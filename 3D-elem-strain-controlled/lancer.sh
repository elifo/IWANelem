mkdir -p outputfiles
rm -f *~ fort.*
gfortran nonlinear3c.f90 water.f90 iwan3c.f90 lu_solver.f90 -o huhu
[[ $? -eq 0 ]] && time ./huhu

