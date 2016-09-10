#!/bin/bash

ulimit -s  102400


INITIAL_PATH=$PWD

#Compile

cd ./tools
gcc -std=c99 -c ./math_algo.c
gcc -std=c99 -c ./sys_pro.c
ar crv tools.a math_algo.o sys_pro.o

cd ../file_io
gcc -std=c99 -c ./file_in.c -I ../include
gcc -std=c99 -c ./file_out.c -I ../include
gcc -std=c99 -c ./io_control.c -I ../include
ar crv file_io.a file_in.o io_control.o file_out.o

cd ../meshing
gcc -std=c99 -c ./mesh_init.c -I ../include
gcc -std=c99 -c ./msh_load.c -I ../include
gcc -std=c99 -c ./quad_mesh.c -I ../include
gcc -std=c99 -c ./line_mesh.c -I ../include
ar crv meshing.a mesh_init.o msh_load.o quad_mesh.o line_mesh.o

cd ../finite_volume
gcc -std=c99 -c ./Euler_scheme.c -I ../include
gcc -std=c99 -c ./cell_init.c -I ../include
ar crv finite_volume.a Euler_scheme.o cell_init.o

cd ../
gcc -std=c99 -c ./hydrocode.c -I ./include
gcc -std=c99 -o hydrocode.out ./hydrocode.o ./finite_volume/finite_volume.a ./meshing/meshing.a ./file_io/file_io.a ./tools/tools.a -lm

#./hydrocode.out Sod_test Sod_test/Sod_test_ROE 1 ROE free_1D_mesh
./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test_ROE 2 ROE Sod_mesh
#./hydrocode.out Sod_10_test2 Sod_10_test2/Sod_10_test_ROE 2 ROE t1

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
#./EUL_source.out odd_even odd_even/odd_even_vfix odd_even_mesh Riemann_exact_vfix -1 0.4


#./EUL_source.out odd_even_entropy_wave odd_even_EW/odd_even_EW_ROE odd_even_mesh ROE -1 0.4

#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream_ROE odd_even_mesh ROE -1 0.4
#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream_Roe_Goundov odd_even_mesh Roe_Goundov -1 0.4
#./EUL_source.out odd_even_EW_upstream odd_even_EW_upstream/odd_even_EW_upstream odd_even_mesh Riemann_exact -1 0.4

#./EUL_source.out drho_upstream drho_upstream/drho_upstream_ROE odd_even_EW_upstream_mesh ROE -1 0.4
#./EUL_source.out drho_upstream drho_upstream/drho_upstream_Roe_Goundov odd_even_EW_upstream_mesh Roe_Goundov -1 0.4

#./EUL_source.out du_upstream du_upstream/du_upstream_ROE odd_even_EW_upstream_mesh ROE -1 0.4
#./EUL_source.out du_upstream du_upstream/du_upstream_Roe_Goundov odd_even_EW_upstream_mesh Roe_Goundov -1 0.4

#./EUL_source.out stationary_shock_OE stationary_shock_OE/stationary_shock_OE_ROE odd_even_EW_upstream_mesh ROE -1 0.4
#./EUL_source.out stationary_shock_OE stationary_shock_OE/stationary_shock_OE_Roe_Goundov odd_even_EW_upstream_mesh Roe_Goundov -1 0.4

#./EUL_source.out EW_SI EW_SI/EW_SI_ROE odd_even_EW_upstream_mesh ROE -1 0.4
#./EUL_source.out EW_SI EW_SI/EW_SI_Roe_Goundov odd_even_EW_upstream_mesh Roe_Goundov -1 0.4

#./EUL_source.out VW_SI VW_SI/VW_SI_ROE odd_even_EW_upstream_mesh ROE -1 0.4
#./EUL_source.out VW_SI VW_SI/VW_SI_Roe_Goundov odd_even_EW_upstream_mesh Roe_Goundov -1 0.4

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

#./EUL_source.out Cylinder_fix Cylinder_fix/Cylinder_fix_ROE Cylinder_mesh ROE -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder_ROE Cylinder_mesh ROE -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder Cylinder_mesh Riemann_exact -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder_Roe_Goundov Cylinder_mesh Roe_Goundov -1 0.4
#./EUL_source.out Cylinder Cylinder/Cylinder_Roe_HLL Cylinder_mesh Roe_HLL -1 0.4

#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad_ROE Free_mesh ROE -1 0.4
#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad Free_mesh Riemann_exact -1 0.4
#./EUL_source.out Riemann_2D3_Quad Riemann_2D3_Quad/Riemann_2D3_Quad_vfix Free_mesh Riemann_exact_vfix -1 0.4

#./EUL_source.out shock-entropy_wave shock-entropy_wave/shock-entropy_wave_ROE Sod_mesh ROE -1 0.4

#./EUL_source.out steady_entropy_wave steady_entropy_wave/steady_entropy_wave_ROE RMI_mesh ROE -1 0.4

#./EUL_source.out NEW_TEST_Sod NEW_TEST_Sod/NEW_TEST_Sod Free_mesh Riemann_exact_vfix -1 0.4
#./EUL_source.out NEW_TEST_Sod NEW_TEST_Sod/NEW_TEST_Sod_R Free_mesh Riemann_exact -1 0.4

#./EUL_source.out NEW_TEST NEW_TEST/NEW_TEST Free_mesh Riemann_exact_vfix -1 0.4
#./EUL_source.out NEW_TEST NEW_TEST/NEW_TEST_R Free_mesh Riemann_exact -1 0.4

#./EUL_source.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE_R Free_mesh Riemann_exact -1 0.4
#./EUL_source.out NEW_TEST_OBLIQUE NEW_TEST_OBLIQUE/NEW_TEST_OBLIQUE Free_mesh Riemann_exact_vfix -1 0.4

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

#./RMI_test.sh


exit 0
