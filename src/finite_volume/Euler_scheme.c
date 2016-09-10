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
	double cpu_time;	
	
	const int num_cell = mv->num_ghost + (int)config[3];

	struct cell_var cv = cell_mem_init(mv);

	cons_qty_init(&cv, FV);

	vol_comp(&cv, mv);

	//	cell_pt_clockwise(mv);

	/*		for(int i = 0; i < num_cell; i++)
			{
				printf("%d 1= %d\n", i, mv->cell_pt[i][1]);
				printf("%d 2= %d\n", i, mv->cell_pt[i][2]);
				printf("%d 3= %d\n", i, mv->cell_pt[i][3]);
				printf("%d 4= %d\n", i, mv->cell_pt[i][4]);
			}
	
	for(int i = 0; i <= mv->num_border[1]; i++)
		{
			printf("%d = %d\n", i, mv->border_pt[i]);
		}
	*/

	cell_rel(&cv, mv);

	printf("Grid has been constructed.\n");
	
	
	for(int i = 0, stop_step = 0; i < (int)config[5] && stop_step == 0; ++i)
		{
			start_clock = clock(); 		
			cpu_time += (clock() - start_clock) / (double)CLOCKS_PER_SEC;

			DispPro(i*100.0/config[5], i);
			
			if (stop_step)
				break;
		}
	printf("\nThe cost of CPU time for the Eulerian method is %g seconds.\n", cpu_time);
}
