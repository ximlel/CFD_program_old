#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
//#include "../include/Riemann_solver.h"
#include "../include/finite_volume.h"


void Euler_scheme(struct flu_var * FV, struct mesh_var * mv, char * scheme)
{
	clock_t start_clock;
	double cpu_time = 0.0;	

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int el = isinf(config[8]);
	const int num_cell = mv->num_ghost + (int)config[3];

	struct cell_var cv = cell_mem_init(mv);

	cons_qty_init(&cv, FV);

	vol_comp(&cv, mv);

	if (dim == 2)
		cell_pt_clockwise(mv);

	cell_rel(&cv, mv);

	if (dim == 2 && order > 1)
		cell_centroid(&cv, mv);

	printf("Grid has been constructed.\n");
	
	for(int i = 0, stop_step = 0; i < (int)config[5] && stop_step == 0; ++i)
		{
			start_clock = clock();

			if (dim == 2 && order > 1 && el != 0 && i > 0)
				cell_centroid(&cv, mv);
			
			if (order > 1)				
				slope_limiter(&cv, mv, FV);
			


			
			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;

			DispPro(i*100.0/config[5], i);
			
			if (stop_step)
				break;
		}
	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
