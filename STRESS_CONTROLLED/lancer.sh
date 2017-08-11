rm -rf outputfiles
mkdir -p outputfiles
rm -f *~ fort.*
gfortran nonlinear3c.f90 water.f90 main.f90  -o huhu
[[ $? -eq 0 ]] && time ./huhu
python plot.py

