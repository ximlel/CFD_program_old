#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"
#include "../include/tools.h"


#define MAX(a,b) ((a) > (b) ? (a) : (b))


#define cv_init_mem(v, n)												\
	do {																\
		cv.v = malloc((n) * sizeof(double));							\
		if(cv.v == NULL)												\
			{															\
				fprintf(stderr, "Not enough memory in cell var init!\n"); \
				goto return_NULL;										\
			}															\
	} while (0)															\

#define cp_init_mem(v, n)												\
	do {																\
		cv.v = malloc((n) * sizeof(double));							\
		if(cv.v == NULL)												\
			{															\
				fprintf(stderr, "Not enough memory in cell var init!\n"); \
				goto return_NULL;										\
			}															\
		init_mem(cv.v, n, mv->cell_pt);									\
	} while (0)															\

#define cp_init_mem_int(v, n)											\
	do {																\
		cv.v = malloc((n) * sizeof(int));								\
		if(cv.v == NULL)												\
			{															\
				fprintf(stderr, "Not enough memory in cell var init!\n"); \
				goto return_NULL;										\
			}															\
		init_mem_int(cv.v, n, mv->cell_pt);								\
	} while (0)															\

struct cell_var cell_mem_init(struct mesh_var * mv)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	const int num_cell = mv->num_ghost + (int)config[3];

	struct cell_var cv;

	cp_init_mem_int(cell_cell, num_cell);
	cv_init_mem(vol, num_cell);
	
	cp_init_mem(n_x, num_cell);
	cp_init_mem(F_u, (int)config[3]);
	cv_init_mem(U_u, num_cell);

	cp_init_mem(F_rho, (int)config[3]);
	cv_init_mem(U_rho, num_cell);

	cp_init_mem(F_e, (int)config[3]);
	cv_init_mem(U_e, num_cell);

	if (order > 1)
		{
			cv_init_mem(gradx_rho, (int)config[3]);
			cv_init_mem(gradx_e, (int)config[3]);
			cv_init_mem(gradx_u, (int)config[3]);			
		}
	
	if (dim > 1)
		{
			cp_init_mem(n_y, num_cell);
			cp_init_mem(F_v, (int)config[3]);
			cv_init_mem(U_v, num_cell);
			if (order > 1)
				{
					cv_init_mem(grady_rho, (int)config[3]);
					cv_init_mem(grady_e, (int)config[3]);
					cv_init_mem(grady_u, (int)config[3]);
					cv_init_mem(grady_v, (int)config[3]);
					cv_init_mem(gradx_v, (int)config[3]);
				}
		}
	if (dim > 2)
		{					
			cp_init_mem(n_z, num_cell);
			cp_init_mem(F_w, (int)config[3]);
			cv_init_mem(U_w, num_cell);		
			if (order > 1)
				{
					cv_init_mem(gradz_rho, (int)config[3]);
					cv_init_mem(gradz_e, (int)config[3]);
					cv_init_mem(gradz_u, (int)config[3]);
					cv_init_mem(gradz_v, (int)config[3]);
					cv_init_mem(gradz_w, (int)config[3]);
					cv_init_mem(grady_w, (int)config[3]);
					cv_init_mem(gradx_w, (int)config[3]);
				}
		}

	if ((int)config[2] == 2)
		{	
			cp_init_mem(F_phi, (int)config[3]);
			cv_init_mem(U_phi, num_cell);				
			if (order > 1)
				{
					cv_init_mem(gradx_phi, (int)config[3]);
					if (dim > 2)
						cv_init_mem(grady_phi, (int)config[3]);
					if (dim > 3)
						cv_init_mem(gradz_phi, (int)config[3]);	
				}
		}

	return cv;
	
 return_NULL:
	exit(5);
}



//Initialize conserved quantities.
void cons_qty_init(struct cell_var * cv, struct flu_var * FV)
{
	const int dim = (int)config[0];
	for(int k = 0; k < (int)config[3]; k++)
		{
			cv->U_rho[k] = FV->RHO[k];
			cv->U_e[k] = FV->P[k]/(FV->gamma[k]-1.0) + 0.5*FV->RHO[k]*FV->U[k]*FV->U[k];
			cv->U_u[k] = FV->RHO[k] * FV->U[k];
			
			if (dim > 1)
				{									
					cv->U_v[k] = FV->RHO[k] * FV->V[k];
					cv->U_e[k] += 0.5*FV->RHO[k]*FV->V[k]*FV->V[k];
				}
			if (dim > 2)
				{									
					cv->U_w[k] = FV->RHO[k] * FV->W[k];
					cv->U_e[k] += 0.5*FV->RHO[k]*FV->W[k]*FV->W[k];
				}

			if ((int)config[2] == 2)					
				cv->U_phi[k] = FV->RHO[k] * FV->PHI[k];
		}
}


//Calculate volume.
void vol_comp(struct cell_var * cv, struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int p_p, p_n;

	for(int k = 0; k < num_cell; k++)
		{			
			cv->vol[k] = 0.0;
			for(int j = 1; j <= mv->cell_pt[k][0]; j++)
				{
					if(j == mv->cell_pt[k][0]) 
						{
							p_p=mv->cell_pt[k][1];
							p_n=mv->cell_pt[k][j];
						}				  
					else
						{
							p_p=mv->cell_pt[k][j+1];
							p_n=mv->cell_pt[k][j];
						} 
					cv->vol[k] = cv->vol[k] + 0.5 * (mv->X[p_n]*mv->Y[p_p] - mv->Y[p_n]*mv->X[p_p]);
				}
		}
}


void cell_pt_clockwise(struct mesh_var * mv)
{
	const int num_cell = mv->num_ghost + (int)config[3];
	int **cp = mv->cell_pt;
	int p_p, p, p_n;
	
	int X_max, n_max;
	for(int k = 0; k < num_cell; k++)
		{
			n_max = 1;
			p = cp[k][n_max];
			X_max = mv->X[p];
			
			for(int j = 2; j <= cp[k][0]; j++)
				{
					n_max = mv->X[cp[k][j]] > X_max ? j : n_max;
					p = cp[k][n_max];
					X_max = mv->X[p];
				}

			if(n_max == cp[k][0]) 
				{
					p_p=cp[k][1];
					p_n=cp[k][n_max-1];
				}
			else if(n_max == 1)
				{
					p_p=cp[k][n_max+1];
					p_n=cp[k][cp[k][0]];
				}
			else
				{
					p_p=cp[k][n_max+1];
					p_n=cp[k][n_max-1];
				}

			if ((mv->X[p_p] - mv->X[p])*(mv->Y[p_n] - mv->Y[p]) - (mv->Y[p_p] - mv->Y[p])*(mv->X[p_n] - mv->X[p]) < 0.0)
				for(int j = 1, temp; j < cp[k][0]/2; j++)
					{
						temp = cp[k][j];
						cp[k][j] = cp[k][cp[k][0]+1-j];
						cp[k][cp[k][0]+1-j] = temp;
					}			
		}
}


//Determine the normal direction and relationship between cells.
void cell_rel(struct cell_var * cv, struct mesh_var * mv)
{	
	const int num_cell = mv->num_ghost + (int)config[3];
	
	int **cp = mv->cell_pt;
	int p_p, p_n, p2_p, p2_n;

	int cell_rec, n_border;
	int i, l, ts;
	double length;	
	for(int k = 0; k < num_cell; k++)
		{ 						
			for(int j = 1; j <= cp[k][0]; j++)
				{
					if(j == cp[k][0]) 
						{
							p_p=cp[k][1];
							p_n=cp[k][j];
						}				  
					else
						{
							p_p=cp[k][j+1];
							p_n=cp[k][j];
						}
					length =  sqrt((mv->Y[p_p] - mv->Y[p_n])*(mv->Y[p_p] - mv->Y[p_n])+(mv->X[p_n] - mv->X[p_p])*(mv->X[p_n] - mv->X[p_p]));
					cv->n_x[k][j] = (mv->Y[p_p] - mv->Y[p_n]) / length;
					cv->n_y[k][j] = (mv->X[p_n] - mv->X[p_p]) / length;
					
					cell_rec = 0;										
					
					ts = 1;
					while (ts <= MAX(num_cell-k-1, k))
						{
							// seek in two side
							i = k + ts;
							if (ts > 0)
								ts = -ts;
							else
								ts = -ts + 1;
							if (i < 0 || i >= num_cell)
								continue;
								
							for(l = 1; l <= cp[i][0]; l++)
								{
									if(l == cp[i][0]) 
										{
											p2_p=cp[i][1];
											p2_n=cp[i][l];
										}				  
									else
										{
											p2_p=cp[i][l+1];
											p2_n=cp[i][l];
										}											
									if((p_p == p2_n) && (p2_p == p_n))
										{
											cv->cell_cell[k][j] = i;
											cell_rec = 1;
											printf("%d %d %d %d\n",k, j, i, l);
											break;
										}
								}
							if (cell_rec)
								break;
						}
					//					printf("%d\n",k);
					if (cell_rec)
						continue;								
												printf("BC:%d\n",k);
					
					for(l = 1, i = -1; l <= mv->num_border[0]; l++)
						{
							i++;
							n_border = i + mv->num_border[l]; 
							for( ; i < n_border; i++)
								{				
									p2_p=mv->border_pt[i+1];
									p2_n=mv->border_pt[i];
									if((p_p == p2_p && p_n == p2_n) || (p_p == p2_n && p_n == p2_p))
										{
											cv->cell_cell[k][j] = mv->border_cond[i];
											cell_rec = 1;
											break;
										}
								}							
							if (cell_rec)
								break;
						}
			
					if(!cell_rec && k < (int)config[3])
						{
							fprintf(stderr, "Ther are some wrong cell relationships!\n");
							exit(2);
						}
				}							
		}
}
