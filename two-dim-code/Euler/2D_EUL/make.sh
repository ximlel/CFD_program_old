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
gcc -c ./CreateDir.c
ar crv custom.a initialize_memory.o rinv.o rhd.o CreateDir.o

cd ../Riemann_solver/
gcc -c ./Riemann_solver_exact.c
gcc -c ./ROE_solver.c -I ../
gcc -c ./Roe_Goundov_solver.c -I ../
gcc -c ./Roe_HLL_solver.c -I ../
gcc -c ./HLL_solver.c -I ../
ar crv Riemann_solver.a Riemann_solver_exact.o ROE_solver.o Roe_Goundov_solver.o Roe_HLL_solver.o HLL_solver.o

cd $INITIAL_PATH/meshing/
gcc -c ./Sod_mesh.c
gcc -c ./Free_mesh.c
gcc -c ./odd_even_mesh.c
gcc -c ./odd_even_EW_mesh.c
gcc -c ./odd_even_EW_upstream_mesh.c
gcc -c ./Shear_mesh.c
gcc -c ./Cylinder_mesh.c
gcc -c ./RMI_mesh.c
ar crv meshing.a Sod_mesh.o Free_mesh.o odd_even_mesh.o odd_even_EW_mesh.o odd_even_EW_upstream_mesh.o Shear_mesh.o Cylinder_mesh.o RMI_mesh.o

cd ../cell_centered_scheme/
gcc -c ./linear_GRP_solver_Edir.c -I ../
gcc -c ./linear_GRP_solver_Edir_2D.c -I ../
gcc -c ./slope_limiter_Ven.c -I ../
gcc -c ./first_order_solver.c -I ../
gcc -c ./first_order_two_species_solver.c -I ../
gcc -c ./second_order_solver.c -I ../
gcc -c ./second_order_two_species_solver.c -I ../
ar crv cell_centered_scheme.a linear_GRP_solver_Edir.o linear_GRP_solver_Edir_2D.o slope_limiter_Ven.o first_order_solver.o first_order_two_species_solver.o second_order_solver.o second_order_two_species_solver.o
ranlib cell_centered_scheme.a

cd ../
gcc -c ./EUL_source.c
gcc -o EUL_source.out ./EUL_source.o ../../lib/file_io/file_io.a ./meshing/meshing.a ./cell_centered_scheme/cell_centered_scheme.a ../../lib/Riemann_solver/Riemann_solver.a ../../lib/custom/custom.a -lm


#Initilize

# cd ../../../data_in/two-dim/
# ./data_initilize.sh
cd $INITIAL_PATH

#Run


## first order

#./EUL_source.out Sod_10 Sod_10/Sod_10_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out Sod_10 Sod_10/Sod_10 Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out Sod_10 Sod_10/Sod_10_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out Sod Sod/Sod_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out Sod Sod/Sod_HLL Sod_mesh HLL -1 0.4
#./EUL_source.out Sod Sod/Sod Sod_mesh Riemann_exact -1 0.4

#./EUL_source.out Test1 Test1/Test1_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out Test1 Test1/Test1 Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out Test1 Test1/Test1_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out steady_shock steady_shock/steady_shock_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out steady_shock steady_shock/steady_shock Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out steady_shock steady_shock/steady_shock_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out steady_shear steady_shear/steady_shear_ROE Shear_mesh ROE -1 0.5
#./EUL_source.out steady_shear steady_shear/steady_shear Shear_mesh Riemann_exact -1 0.5
#./EUL_source.out steady_shear steady_shear/steady_shear_Roe_Goundov Shear_mesh Roe_Goundov -1 0.5

#./EUL_source.out contact_only contact_only/contact_only_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out contact_only contact_only/contact_only Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out contact_only contact_only/contact_only_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out odd_even odd_even/odd_even_ROE odd_even_mesh ROE -1 0.4
#./EUL_source.out odd_even odd_even/odd_even_HLL odd_even_mesh HLL -1 0.4
#./EUL_source.out odd_even odd_even/odd_even odd_even_mesh Riemann_exact -1 0.4
#./EUL_source.out odd_even odd_even/odd_even_Roe_Goundov odd_even_mesh Roe_Goundov -1 0.4
#./EUL_source.out odd_even odd_even/odd_even_Roe_HLL odd_even_mesh Roe_HLL -1 0.4

#./EUL_source.out odd_even_entropy_wave odd_even_EW/odd_even_EW_ROE odd_even_mesh ROE -1 0.4

#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream_ROE odd_even_mesh ROE -1 0.4
#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream_Roe_Goundov odd_even_mesh Roe_Goundov -1 0.4
#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream odd_even_mesh Riemann_exact -1 0.4

#./EUL_source.out drho_upstream drho_upstream/drho_upstream_ROE odd_even_EW_upstream_mesh ROE -1 0.4

#./EUL_source.out du_upstream du_upstream/du_upstream_ROE odd_even_EW_upstream_mesh ROE -1 0.4

#./EUL_source.out stationary_shock_OE stationary_shock_OE/stationary_shock_OE_ROE odd_even_EW_upstream_mesh ROE -1 0.4

./EUL_source.out stationary_shock_OE stationary_shock_OE/stationary_shock_OE_Roe_Goundov odd_even_EW_upstream_mesh Roe_Goundov -1 0.4

#./EUL_source.out no_stationary_shock_OE no_stationary_shock_OE/no_stationary_shock_OE_ROE odd_even_EW_upstream_mesh ROE -1 0.4

#./EUL_source.out entropy_wave_upstream entropy_wave_upstream/entropy_wave_upstream_ROE Sod_mesh ROE -1 0.4
#./EUL_source.out entropy_wave_upstream entropy_wave_upstream/entropy_wave_upstream Sod_mesh Riemann_exact -1 0.4
#./EUL_source.out entropy_wave_upstream entropy_wave_upstream/entropy_wave_upstream_Roe_Goundov Sod_mesh Roe_Goundov -1 0.4

#./EUL_source.out shear_wave_upstream shear_wave_upstream/shear_wave_upstream_ROE Sod_mesh ROE -1 0.4

#./EUL_source.out one_line_du one_line_du/one_line_du_ROE Free_mesh ROE -1 0.4
#./EUL_source.out one_line_du one_line_du/one_line_du_HLL Free_mesh HLL -1 0.4
#./EUL_source.out one_line_du one_line_du/one_line_du Free_mesh Riemann_exact -1 0.4
#./EUL_source.out one_line_du one_line_du/one_line_du_Roe_Goundov Free_mesh Roe_Goundov -1 0.4

#./EUL_source.out one_line_dp one_line_dp/one_line_dp_ROE Free_mesh ROE -1 0.4
#./EUL_source.out one_line_dp one_line_dp/one_line_dp Free_mesh Riemann_exact -1 0.4
#./EUL_source.out one_line_dp one_line_dp/one_line_dp_Roe_Goundov Free_mesh Roe_Goundov -1 0.4

#./EUL_source.out Cylinder Cylinder/Cylinder_ROE Cylinder_mesh ROE -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder Cylinder_mesh Riemann_exact -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder_Roe_Goundov Cylinder_mesh Roe_Goundov -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder_Roe_HLL Cylinder_mesh Roe_HLL -1 0.4

#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad_ROE Free_mesh ROE -1 0.4
#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad Free_mesh Riemann_exact -1 0.4

#./EUL_source.out shock-entropy_wave shock-entropy_wave/shock-entropy_wave_ROE Sod_mesh ROE -1 0.4

#./EUL_source.out steady_entropy_wave steady_entropy_wave/steady_entropy_wave_ROE RMI_mesh ROE -1 0.4

### Two_species

#./EUL_source.out RMI/RMI_origin RMI/RMI_Roe_Goundov RMI_mesh Roe_Goundov 1000 0.4
#./EUL_source.out RMI/RMI_origin RMI/RMI RMI_mesh Riemann_exact -1 0.4 Two_species
#./EUL_source.out RMI/RMI_origin RMI/RMI_single RMI_mesh Riemann_exact -1 0.4
#./EUL_source.out RMI/RMI_321 RMI/RMI_321 RMI_mesh Riemann_exact -1 0.4 Two_species


## second order

#./EUL_source.out Sod_10 Sod_10 Sod_mesh GRP -1 0.4 second_order

#./EUL_source.out RMI/RMI_origin RMI/RMI RMI_mesh GRP -1 0.4 second_order Two_species
#./EUL_source.out RMI/RMI_origin RMI/RMI_single RMI_mesh GRP -1 0.4 second_order
#./EUL_source.out RMI/RMI_641 RMI/RMI_641 RMI_mesh GRP -1 0.4 second_order Two_species
#./EUL_source.out RMI/RMI_641 RMI/RMI_641_single RMI_mesh GRP -1 0.4 second_order

#./EUL_source.out RMI/RMI_161 RMI_bp/RMI_161 RMI_mesh GRP -1 0.4 second_order Two_species
#./EUL_source.out RMI/RMI_K1A10 RMI_bp/RMI_K1A10 RMI_mesh GRP -1 0.4 second_order Two_species

./RMI_test.sh


exit 0
