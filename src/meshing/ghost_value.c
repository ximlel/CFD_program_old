#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"



void period_ghost(struct cell_var * cv, struct mesh_var mv, double t)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int num_cell_ghost = mv.num_ghost + (int)config[3];
	const int num_cell = (int)config[3];

	const int *per_cell = mv.peri_cell;	
	
	for(int i = num_cell; i < num_cell_ghost; i++)
		{
			cv->U_rho[i] = cv->U_rho[per_cell[i]];
			cv->U_e[i] = cv->U_e[per_cell[i]];
			cv->U_u[i] = cv->U_u[per_cell[i]];
			if (order > 1)
				{
					cv->gradx_rho[i] = cv->gradx_rho[per_cell[i]];
					cv->gradx_e[i] = cv->gradx_e[per_cell[i]];
					cv->gradx_u[i] = cv->gradx_u[per_cell[i]];
				}
			if (dim > 1)
				{
					cv->U_v[i] = cv->U_v[per_cell[i]];
					if (order > 1)
						{
							cv->grady_rho[i] = cv->grady_rho[per_cell[i]];
							cv->grady_e[i] = cv->grady_e[per_cell[i]];
							cv->grady_u[i] = cv->grady_u[per_cell[i]];
							cv->grady_v[i] = cv->grady_v[per_cell[i]];
							cv->gradx_v[i] = cv->gradx_v[per_cell[i]];
						}
				}
			if (dim > 2)
				{
					cv->U_w[i] = cv->U_w[per_cell[i]];
					if (order > 1)
						{
							cv->gradz_rho[i] = cv->gradz_rho[per_cell[i]];
							cv->gradz_e[i] = cv->gradz_e[per_cell[i]];
							cv->gradz_u[i] = cv->gradz_u[per_cell[i]];
							cv->gradz_v[i] = cv->gradz_v[per_cell[i]];
							cv->gradz_w[i] = cv->gradz_w[per_cell[i]];
							cv->grady_w[i] = cv->grady_w[per_cell[i]];
							cv->gradx_w[i] = cv->gradx_w[per_cell[i]];
						}
				}
			if ((int)config[2] == 2)
				{
					cv->U_phi[i] = cv->U_phi[per_cell[i]];
					if (order > 1)
						{
							cv->gradx_phi[i] = cv->gradx_phi[per_cell[i]];
							if (dim > 2)
								cv->grady_phi[i] = cv->grady_phi[per_cell[i]];
							if (dim > 3)
								cv->gradz_phi[i] = cv->gradz_phi[per_cell[i]];
						}
				}
			cv->gamma[i] = cv->gamma[per_cell[i]];
		}
}



void period_cell_modi(struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];

	int *pc = mv->peri_cell;	
	int per_num[num_cell], per_n = 0;
	for (int i = 0; i < num_cell; i++)
		{
			if (pc[i] >= 0)
				per_n++;

			per_num[i] == per_n;
		}

	for (int i = 0; i < num_cell; i++)
		if (pc[i] >= 0)
			pc[i] -= per_num[pc[i]];

	int *cc_tmp, pc_tmp; 
	for (int i = 1, j; i < num_cell; i++)
		{
			j = i;
			while (pc[j-1] < 0 || j < 1)
				{
					pc_tmp = pc[j-1];
					pc[j-1] = pc[j];
					pc[j] = pc_tmp;
					cc_tmp = mv->cell_pt[j];
					mv->cell_pt[j] = mv->cell_pt[j-1];
					mv->cell_pt[j-1] = cc_tmp;
				}	 				
		}
	mv->bc = period_ghost;
}
