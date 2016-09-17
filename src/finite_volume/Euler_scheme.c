#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "../include/var_struc.h"
#include "../include/tools.h"
//#include "../include/Riemann_solver.h"
#include "../include/finite_volume.h"


void Euler_scheme(struct flu_var *FV, const struct mesh_var mv, const char *scheme)
{
	clock_t start_clock;
	double cpu_time = 0.0;	

	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int el = isinf(config[8]);
	const int num_cell = (int)config[3];

	struct cell_var cv = cell_mem_init(mv);

	cons_qty_init(&cv, *FV);

	vol_comp(&cv, mv);

	cell_rel(&cv, mv);

	if (dim == 2 && order > 1)
		cell_centroid(&cv, mv);

	printf("Grid has been constructed.\n");


	
	double tau; // the length of the time step
	double t_all = 0.0;

	struct i_f_var ifv, ifv_R;

	double lambda_max, P;
	
	int k, j, ivi;
	int p_p, p_n;
	int ** cp = mv.cell_pt;
	for(int i = 0, stop_step = 0; i < (int)config[5] && stop_step == 0; ++i)
		{
			start_clock = clock();
			
			if (dim == 2 && order > 1 && el != 0 && i > 0)
				cell_centroid(&cv, mv);
			
			if (order > 1)				
				slope_limiter(&cv, mv, *FV);

			if (mv.bc != NULL)
				mv.bc(&cv, mv, t_all);
			
			if (order > 1 && dim == 2)
				tau_calc(cv, mv);

			for(k = 0; k < num_cell; k++)
				{
					if (stop_step == 1)
						break;
					for(j = 1; j <= cp[k][0]; j++)
						{
							ivi = interface_var_init(cv, mv, &ifv, &ifv_R, k, j);
																printf("%d\n", ivi);
							if (ivi == 0)
								break;
							else if(ivi == -1)
								{
									stop_step = 1;
									break;									
								}
						}
				}
			
			DispPro(i*100.0/config[5], i);
			
			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;
			
			if (stop_step == 1)
				{
					printf("NOW : %d", i);
					break;
				}
		}
	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
