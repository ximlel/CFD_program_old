#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "../include/var_struc.h"
#include "../include/tools.h"


static inline double mu_BJ(double x)
{
	return (x<1?x:1);
}

static inline double mu_Ven(double x)
{
	return ((x*x+2*x)/(x*x+x+2));
}


static void lsq_limiter
(struct cell_var * cv, struct mesh_var * mv, 
 double * grad_W_x, double * grad_W_y, double * W)
{
	const double eps = config[4];
	const int num_cell = (int)config[3];
	const int n_x = (int)config[13], n_y = (int)config[14];
	const int lim = isinf(config[40]) ? 1 : (int)config[40]; //limiter
	double (*mu[])(double) = { mu_BJ, mu_Ven };
	
	int **cp = mv->cell_pt;
	int **cc = cv->cell_cell;
	const double *X_c = cv->X_c;
	const double *Y_c = cv->Y_c;	
	const double *X = mv->X;
	const double *Y = mv->Y;

	int cell_R;
	double tmp_x, tmp_y;
	double M_c[2][2], M_c_i[2][2];
	double W_c_min, W_c_max, W_c_x_p;
	double fai_W;
	
	for(int k = 0; k < num_cell; ++k)
		{
			M_c[0][0] = 0.0;  M_c[0][1] = 0.0;
			M_c[1][0] = 0.0;  M_c[1][1] = 0.0;
			grad_W_x[k] = 0.0;
			grad_W_y[k] = 0.0;

			for(int j = 0; j < cc[k][0]; ++j)
				{
					if (cc[k][j] >= 0)
						cell_R = cc[k][j];
					else if (cc[k][j] == -1 || cc[k][j] == -2 || cc[k][j] == -3)
						continue;
					else
						{
							fprintf(stderr, "No suitable boundary!\n");
							exit(2);
						}
							
					M_c[0][0] += (X_c[cell_R] - X_c[k]) * (X_c[cell_R] - X_c[k]);
					M_c[0][1] += (X_c[cell_R] - X_c[k]) * (Y_c[cell_R] - Y_c[k]);
					M_c[1][0] += (Y_c[cell_R] - Y_c[k]) * (X_c[cell_R] - X_c[k]);
					M_c[1][1] += (Y_c[cell_R] - Y_c[k]) * (Y_c[cell_R] - Y_c[k]);
					grad_W_x[k] += (W[cell_R] - W[k]) * (X_c[cell_R] - X_c[k]);
					grad_W_y[k] += (W[cell_R] - W[k]) * (Y_c[cell_R] - Y_c[k]);
				}
								
								M_c_i[0][0] =   M_c[1][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
								M_c_i[0][1] = - M_c[0][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
								M_c_i[1][0] = - M_c[1][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
								M_c_i[1][1] =   M_c[0][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
			/*
			M_c_i[0][0] = M_c[0][0];  M_c_i[0][1] = M_c[0][1];
			M_c_i[1][0] = M_c[1][0];  M_c_i[1][1] = M_c[1][1];
			if(rinv(M_c_i[0], 2) == 0)
				exit(3);
			*/

			tmp_x = M_c_i[0][0] * grad_W_x[k] + M_c_i[0][1] * grad_W_y[k];
			tmp_y = M_c_i[1][0] * grad_W_x[k] + M_c_i[1][1] * grad_W_y[k];
			grad_W_x[k] = tmp_x;
			grad_W_y[k] = tmp_y;						
		}

	for(int k = 0; k < num_cell; ++k)
		{
			W_c_min = W[k];
			W_c_max = W[k];
			for(int j = 1; j <= cc[k][0]; ++j)
				{
					if (cc[k][j] >= 0)						
						cell_R = cc[k][j];
					else if (cc[k][j] == -1 || cc[k][j] == -2 || cc[k][j] == -3)
						continue;
					else
						{
							printf("No suitable boundary!\n");
							exit(2);
						}	
							
					if(W[cell_R] < W_c_min)
						W_c_min = W[cell_R];
					else if(W[cell_R] > W_c_max)
						W_c_max = W[cell_R];
				}
			fai_W = 1.0;			
			for(int j = 1; j <= cp[k][0]; ++j)
				{							
					W_c_x_p = W[k] + grad_W_x[k] * (X[cp[k][j]] - X_c[k]) + grad_W_y[k] * (Y[cp[k][j]] - Y_c[k]);
					if (fabs(W_c_x_p - W[k]) < eps)
						;	
					else if((W_c_x_p - W[k]) > 0.0)
						fai_W = fmin(fai_W, mu[lim]((W_c_max - W[k])/(W_c_x_p - W[k])));
					else
						fai_W = fmin(fai_W, mu[lim]((W_c_min - W[k])/(W_c_x_p - W[k])));
				}
			grad_W_x[k] = grad_W_x[k] * fai_W;
			grad_W_y[k] = grad_W_y[k] * fai_W;					   	
		}
}

void slope_limiter(struct cell_var * cv, struct mesh_var * mv, struct flu_var * FV)
{
	const int dim = (int)config[0];

	if (dim > 1)
		{
			lsq_limiter(cv, mv, cv->gradx_rho, cv->grady_rho, cv->U_rho);
			lsq_limiter(cv, mv, cv->gradx_e, cv->grady_e, cv->U_e);
			lsq_limiter(cv, mv, cv->gradx_u, cv->grady_u, cv->U_u);
			lsq_limiter(cv, mv, cv->gradx_v, cv->grady_v, cv->U_v);		
			if ((int)config[2] == 2)
				lsq_limiter(cv, mv, cv->gradx_phi, cv->grady_phi, cv->U_phi);	
		}
}
