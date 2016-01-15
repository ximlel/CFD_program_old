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
gcc -c ./HLL_solver.c -I ../
ar crv Riemann_solver.a Riemann_solver_exact.o ROE_solver.o Roe_Goundov_solver.o HLL_solver.o

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

#./EUL_source.out odd_even_entropy_wave odd_even_EW/odd_even_EW_ROE odd_even_mesh ROE -1 0.4
#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream_ROE odd_even_mesh ROE -1 0.4

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
#=============================RMI KA test=========================
:<<KA
nohup ./EUL_source.out RMI/RMI_K1A10 RMI/RMI_K1A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K1A20 RMI/RMI_K1A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K1A40 RMI/RMI_K1A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K1A80 RMI/RMI_K1A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K1A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A10 RMI/RMI_K2A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A20 RMI/RMI_K2A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A40 RMI/RMI_K2A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K2A80 RMI/RMI_K2A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K2A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A10 RMI/RMI_K4A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A20 RMI/RMI_K4A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A40 RMI/RMI_K4A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K4A80 RMI/RMI_K4A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K4A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A10 RMI/RMI_K8A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A20 RMI/RMI_K8A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A40 RMI/RMI_K8A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K8A80 RMI/RMI_K8A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K8A80.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A10 RMI/RMI_K16A10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A10.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A20 RMI/RMI_K16A20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A20.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A40 RMI/RMI_K16A40 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A40.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_K16A80 RMI/RMI_K16A80 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_K16A80.out 2>&1 &
KA
#=============================RMI KA test=========================
#./EUL_source.out RMI/RMI_161 RMI_bp/RMI_161 RMI_mesh GRP -1 0.4 second_order Two_species
#./EUL_source.out RMI/RMI_K1A10 RMI_bp/RMI_K1A10 RMI_mesh GRP -1 0.4 second_order Two_species
#=============================RMI lnKA more test==========================
#:<<lnKA
nohup ./EUL_source.out RMI/RMI_lnKA_-3 RMI/RMI_lnKA_-3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.9 RMI/RMI_lnKA_-2.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.8 RMI/RMI_lnKA_-2.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.7 RMI/RMI_lnKA_-2.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.6 RMI/RMI_lnKA_-2.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.5 RMI/RMI_lnKA_-2.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.4 RMI/RMI_lnKA_-2.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.3 RMI/RMI_lnKA_-2.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.2 RMI/RMI_lnKA_-2.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2.1 RMI/RMI_lnKA_-2.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-2 RMI/RMI_lnKA_-2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.9 RMI/RMI_lnKA_-1.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.8 RMI/RMI_lnKA_-1.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.7 RMI/RMI_lnKA_-1.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.6 RMI/RMI_lnKA_-1.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.5 RMI/RMI_lnKA_-1.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.4 RMI/RMI_lnKA_-1.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.3 RMI/RMI_lnKA_-1.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.2 RMI/RMI_lnKA_-1.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1.1 RMI/RMI_lnKA_-1.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-1 RMI/RMI_lnKA_-1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.9 RMI/RMI_lnKA_-0.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.8 RMI/RMI_lnKA_-0.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.7 RMI/RMI_lnKA_-0.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.6 RMI/RMI_lnKA_-0.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.5 RMI/RMI_lnKA_-0.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.4 RMI/RMI_lnKA_-0.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.3 RMI/RMI_lnKA_-0.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.2 RMI/RMI_lnKA_-0.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_-0.1 RMI/RMI_lnKA_-0.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_-0.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0 RMI/RMI_lnKA_0 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.1 RMI/RMI_lnKA_0.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.2 RMI/RMI_lnKA_0.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.3 RMI/RMI_lnKA_0.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.4 RMI/RMI_lnKA_0.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.5 RMI/RMI_lnKA_0.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.6 RMI/RMI_lnKA_0.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.7 RMI/RMI_lnKA_0.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.8 RMI/RMI_lnKA_0.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_0.9 RMI/RMI_lnKA_0.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_0.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1 RMI/RMI_lnKA_1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.1 RMI/RMI_lnKA_1.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.2 RMI/RMI_lnKA_1.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.3 RMI/RMI_lnKA_1.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.4 RMI/RMI_lnKA_1.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.5 RMI/RMI_lnKA_1.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.6 RMI/RMI_lnKA_1.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.7 RMI/RMI_lnKA_1.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.8 RMI/RMI_lnKA_1.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_1.9 RMI/RMI_lnKA_1.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_1.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2 RMI/RMI_lnKA_2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.1 RMI/RMI_lnKA_2.1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.1.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.2 RMI/RMI_lnKA_2.2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.2.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.3 RMI/RMI_lnKA_2.3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.3.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.4 RMI/RMI_lnKA_2.4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.4.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.5 RMI/RMI_lnKA_2.5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.5.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.6 RMI/RMI_lnKA_2.6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.6.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.7 RMI/RMI_lnKA_2.7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.7.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.8 RMI/RMI_lnKA_2.8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.8.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_2.9 RMI/RMI_lnKA_2.9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_2.9.out 2>&1 &
nohup ./EUL_source.out RMI/RMI_lnKA_3 RMI/RMI_lnKA_3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_record/RMI_lnKA_3.out 2>&1 &
#lnKA
#=============================RMI lnKA more test==========================

#=============================RMI At more test==========================
:<<At
nohup ./EUL_source.out RMI_At/RMI_At_1 RMI_At/RMI_At_1 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_1.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_2 RMI_At/RMI_At_2 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_2.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_3 RMI_At/RMI_At_3 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_3.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_4 RMI_At/RMI_At_4 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_4.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_5 RMI_At/RMI_At_5 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_5.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_6 RMI_At/RMI_At_6 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_6.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_7 RMI_At/RMI_At_7 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_7.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_8 RMI_At/RMI_At_8 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_8.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_9 RMI_At/RMI_At_9 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_9.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_10 RMI_At/RMI_At_10 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_10.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_11 RMI_At/RMI_At_11 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_11.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_12 RMI_At/RMI_At_12 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_12.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_13 RMI_At/RMI_At_13 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_13.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_14 RMI_At/RMI_At_14 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_14.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_15 RMI_At/RMI_At_15 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_15.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_16 RMI_At/RMI_At_16 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_16.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_17 RMI_At/RMI_At_17 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_17.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_18 RMI_At/RMI_At_18 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_18.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_19 RMI_At/RMI_At_19 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_19.out 2>&1 &
nohup ./EUL_source.out RMI_At/RMI_At_20 RMI_At/RMI_At_20 RMI_mesh GRP -1 0.4 second_order Two_species > ./RMI_At_record/RMI_At_20.out 2>&1 &
At
#=============================RMI At more test==========================



exit 0
