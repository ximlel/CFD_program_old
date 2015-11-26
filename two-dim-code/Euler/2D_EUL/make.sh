#!/bin/bash

ulimit -s  102400


INITIAL_PATH=$PWD

#Compile

cd ../../lib/file_io/
gcc -c ./comman.c
gcc -c ./file_io.c -I ../
ar crv file_io.a comman.o file_io.o

cd ../custom/
gcc -c ./initialize_memory.c
gcc -c ./rinv.c
gcc -c ./rhd.c
ar crv custom.a initialize_memory.o rinv.o rhd.o

cd ../Riemann_solver/
gcc -c ./Riemann_solver_exact.c -g
gcc -c ./ROE_solver.c -I ../
gcc -c ./Roe_Goundov_solver.c -I ../
gcc -c ./HLL_solver.c -I ../
ar crv Riemann_solver.a Riemann_solver_exact.o ROE_solver.o Roe_Goundov_solver.o HLL_solver.o

cd $INITIAL_PATH/meshing/
gcc -c ./Sod_mesh.c
gcc -c ./Free_mesh.c
gcc -c ./odd_even_mesh.c
gcc -c ./Shear_mesh.c
gcc -c ./Cylinder_mesh.c
gcc -c ./RMI_mesh.c
ar crv meshing.a Sod_mesh.o Free_mesh.o odd_even_mesh.o Shear_mesh.o Cylinder_mesh.o RMI_mesh.o

cd ../cell_centered_scheme/
gcc -c ./linear_GRP_solver_Edir.c -I ../ -g
gcc -c ./first_order_solver.c -I ../
gcc -c ./first_order_two_species_solver.c -I ../
ar crv cell_centered_scheme.a linear_GRP_solver_Edir.o first_order_solver.o first_order_two_species_solver.o
ranlib cell_centered_scheme.a

cd ../
gcc -c ./EUL_source.c -g
gcc -o EUL_source.out ./EUL_source.o ../../lib/file_io/file_io.a ./meshing/meshing.a ./cell_centered_scheme/cell_centered_scheme.a ../../lib/Riemann_solver/Riemann_solver.a ../../lib/custom/custom.a -lm


#Initilize

cd ../../../data_in/two-dim/
./data_initilize.sh
cd $INITIAL_PATH

#Run

#./EUL_source.out Sod_10 Sod_10_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out Sod_10 Sod_10 Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out Sod_10 Sod_10_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out Sod Sod_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out Sod Sod_HLL Sod_mesh HLL -1 0.4
#./EUL_source.out Sod Sod Sod_mesh Riemann_exact -1 0.4

#./EUL_source.out Test1 Test1_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out Test1 Test1 Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out Test1 Test1_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out steady_shock steady_shock_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out steady_shock steady_shock Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out steady_shock steady_shock_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out steady_shear steady_shear_ROE Shear_mesh ROE -1 0.5
#./EUL_source.out steady_shear steady_shear Shear_mesh Riemann_exact -1 0.5
#./EUL_source.out steady_shear steady_shear_Roe_Goundov Shear_mesh Roe_Goundov -1 0.5

#./EUL_source.out contact_only contact_only_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out contact_only contact_only Sod_mesh Riemann_exact -1 0.4

#./EUL_source.out odd_even odd_even_ROE odd_even_mesh ROE -1 0.4
#./EUL_source.out odd_even odd_even_HLL odd_even_mesh HLL -1 0.4
#./EUL_source.out odd_even odd_even odd_even_mesh Riemann_exact -1 0.4
#./EUL_source.out odd_even odd_even_Roe_Goundov odd_even_mesh Roe_Goundov -1 0.4

#./EUL_source.out one_line_du one_line_du_ROE Free_mesh ROE -1 0.4
#./EUL_source.out one_line_du one_line_du_HLL Free_mesh HLL -1 0.4
#./EUL_source.out one_line_du one_line_du Free_mesh Riemann_exact -1 0.4

#./EUL_source.out one_line_dp one_line_dp_ROE Free_mesh ROE -1 0.4
#./EUL_source.out one_line_dp one_line_dp Free_mesh Riemann_exact -1 0.4

#./EUL_source.out Cylinder Cylinder_ROE Cylinder_mesh ROE -1 0.4
#./EUL_source.out Cylinder Cylinder Cylinder_mesh Riemann_exact -1 0.4
#./EUL_source.out Cylinder Cylinder_Roe_Goundov Cylinder_mesh Roe_Goundov -1 0.4

#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad_ROE Free_mesh ROE -1 0.4
#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad Free_mesh Riemann_exact -1 0.4


./EUL_source.out RMI RMI RMI_mesh Riemann_exact -1 0.4 Two_species


exit 0
