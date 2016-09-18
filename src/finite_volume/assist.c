#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"


static int cons2prim(struct i_f_var * ifv)
{
	const int dim = (int)config[0];
	const double eps = (int)config[4];
	
	ifv->RHO = ifv->U_rho;
	ifv->U   = ifv->U_u/ifv->U_rho;
	ifv->P   = (ifv->U_e - 0.5*(ifv->U_u*ifv->U_u)/ifv->U_rho) * (ifv->gamma-1.0);
	if (dim > 1)
		{
			ifv->V  = ifv->U_v/ifv->U_rho;
			ifv->P -= (0.5*(ifv->U_v*ifv->U_v)/ifv->U_rho) * (ifv->gamma-1.0);
			if (isnan(ifv->V) || isinf(ifv->V))									
				return 0;
		}
	if (dim > 2)
		{			
			ifv->W  = ifv->U_w/ifv->U_rho;
			ifv->P -= (0.5*(ifv->U_w*ifv->U_w)/ifv->U_rho) * (ifv->gamma-1.0);
			if (isnan(ifv->W) || isinf(ifv->W))
				return 0;
		}
	if ((int)config[2] == 2)
		{
			ifv->PHI = ifv->U_phi/ifv->U_rho;
			if (isnan(ifv->PHI) || ifv->PHI < -10*eps || ifv->PHI > 1.0 + 10*eps)
				return 0;
		}

	if (isnan(ifv->RHO + ifv->U + ifv->P) || isinf(ifv->RHO + ifv->U + ifv->P) || ifv->RHO < -10*eps || ifv->P < -10*eps)
		return 0;

	return 1;
}


void fluid_var_update(struct flu_var *FV, struct cell_var cv)
{
	const int dim = (int)config[0];		
	const int num_cell = (int)config[3];
	struct i_f_var ifv;
	
	for(int k = 0; k < num_cell; k++)
		{					
			ifv.U_rho = cv.U_rho[k];
			ifv.U_e   = cv.U_e[k];
			ifv.U_u   = cv.U_u[k];
			if (dim > 1)
				ifv.U_v = cv.U_v[k];
			if (dim > 2)
				ifv.U_w = cv.U_w[k];
			if ((int)config[2] == 2)
				ifv.U_phi = cv.U_phi[k];
			ifv.gamma   = cv.gamma[k];

			if(cons2prim(&ifv) == 0)
				printf("WROUNG!\n");


			FV->RHO[k] = ifv.RHO;
			FV->P[k]   = ifv.P;
			FV->U[k]   = ifv.U;
			if (dim > 1)
				FV->V[k] = ifv.V;
			if (dim > 2)
				FV->W[k] = ifv.W;
			if ((int)config[2] == 2)
				FV->PHI[k] = ifv.PHI;
			FV->gamma[k] = ifv.gamma;
		}
}


static void order2_i_f_var_init(const struct cell_var cv, struct i_f_var * ifv, const int k)
{
	const int dim = (int)config[0];
	
	const double n_x = ifv->n_x, n_y = ifv->n_y, n_z = ifv->n_z;
	const double delta_x = ifv->delta_x, delta_y = ifv->delta_y, delta_z = ifv->delta_z;
												
	ifv->d_rho  = cv.gradx_rho[k]*n_x;
	ifv->U_rho += cv.gradx_rho[k]*delta_x;
	ifv->d_e    = cv.gradx_e[k]  *n_x;
	ifv->U_e   += cv.gradx_e[k]  *delta_x;
	ifv->d_u    = cv.gradx_u[k]  *n_x;
	ifv->U_u   += cv.gradx_u[k]  *delta_x;
	if ((int)config[2] == 2)
		{									
			ifv->d_phi  = cv.gradx_phi[k]*n_x;
			ifv->U_phi += cv.gradx_phi[k]*delta_x;
		}

	if (dim > 1)
		{
			ifv->d_rho += cv.grady_rho[k]*n_y;
			ifv->U_rho += cv.grady_rho[k]*delta_y;
			ifv->d_e   += cv.grady_e[k]  *n_y;
			ifv->U_e   += cv.grady_e[k]  *delta_y;
			ifv->d_u   += cv.grady_u[k]  *n_y;
			ifv->U_u   += cv.grady_u[k]  *delta_y;
			ifv->d_v    = cv.gradx_v[k]  *n_x     + cv.grady_v[k]*n_y;
			ifv->U_v   += cv.gradx_v[k]  *delta_x + cv.grady_v[k]*delta_y;
			if ((int)config[2] == 2)
				{
					ifv->d_phi += cv.grady_phi[k]*n_y;
					ifv->U_phi += cv.grady_phi[k]*delta_y;
				}
		}
	if (dim > 2)
		{
			ifv->d_rho += cv.gradz_rho[k]*n_z;
			ifv->U_rho += cv.gradz_rho[k]*delta_z;
			ifv->d_e   += cv.gradz_e[k]  *n_z;
			ifv->U_e   += cv.gradz_e[k]  *delta_z;
			ifv->d_u   += cv.gradz_u[k]  *n_z;
			ifv->U_u   += cv.gradz_u[k]  *delta_z;
			ifv->d_v   += cv.gradz_v[k]  *n_z;
			ifv->U_v   += cv.gradz_v[k]  *delta_z;
			ifv->d_w    = cv.gradx_w[k]  *n_x     + cv.grady_w[k]*n_y     + cv.gradz_w[k]*n_z;
			ifv->U_w   += cv.gradx_w[k]  *delta_x + cv.grady_w[k]*delta_y + cv.gradz_w[k]*delta_z;
			if ((int)config[2] == 2)
				{													
					ifv->d_phi += cv.gradz_phi[k]*n_z;
					ifv->U_phi += cv.gradz_phi[k]*delta_z;				
				}
		}
}

static void order2_i_f_var0(struct i_f_var * ifv)
{
	const int dim = (int)config[0];
		
	ifv->d_rho= 0.0;
	ifv->d_e  = 0.0;
	ifv->d_u  = 0.0;
	if (dim > 1)						
		ifv->d_v = 0.0;
	if (dim > 2)
		ifv->d_w  = 0.0;				
	if ((int)config[2] == 2)
		ifv->d_phi = 0.0;
}

int interface_var_init
(const struct cell_var cv, const struct mesh_var mv, struct i_f_var * ifv,
 struct i_f_var * ifv_R, const int k, const int j)
{
	const int dim = (int)config[0];
	const int order = (int)config[9];
	int **cc = cv.cell_cell;
	int **cp = mv.cell_pt;	

	int p_p, p_n;
	if (dim == 2)
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
			ifv->length = sqrt((mv.X[p_p] - mv.X[p_n])*(mv.X[p_p] - mv.X[p_n]) + (mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p] - mv.Y[p_n]));
		}
	ifv->n_x = cv.n_x[k][j];
	ifv->U_rho = cv.U_rho[k];
	ifv->U_e = cv.U_e[k];
	ifv->U_u = cv.U_u[k];
	if (dim > 1)
		{			
			ifv->n_y = cv.n_y[k][j];
 			ifv->U_v = cv.U_v[k];
		}
	if (dim > 2)
		{
			ifv->n_z = cv.n_z[k][j];
			ifv->U_w = cv.U_w[k];		
		}
	if ((int)config[2] == 2)
		ifv->U_phi = cv.U_phi[k];
	ifv->gamma = cv.gamma[k];
	
	
	if (order == 2)
		{
			if (dim == 1)
				ifv->delta_x = mv.X[cp[k][j]] - cv.X_c[k];
			else if (dim == 2)
				{					

					ifv->delta_x = 0.5*(mv.X[p_p] + mv.X[p_n]) - cv.X_c[k];
					ifv->delta_y = 0.5*(mv.Y[p_p] + mv.Y[p_n]) - cv.Y_c[k];
				}
			
			order2_i_f_var_init(cv, ifv, k);			
		}
	
	int cR; //cell_right
	
	ifv_R->n_x = ifv->n_x;
	if (dim > 1)			
		ifv_R->n_y = ifv->n_y;
	if (dim > 2)
		ifv_R->n_z = ifv->n_z;
	
	if (cc[k][j] >= 0)
		{
			cR = cc[k][j];
			ifv_R->U_rho = cv.U_rho[cR];
			ifv_R->U_e = cv.U_e[cR];
			ifv_R->U_u = cv.U_u[cR];
			if (dim > 1)
				ifv_R->U_v = cv.U_v[cR];
			if (dim > 2)
				ifv_R->U_w = cv.U_w[cR];
			if ((int)config[2] == 2)
				ifv_R->U_phi = cv.U_phi[cR];
			ifv_R->gamma = cv.gamma[cR];

			if (order == 2)
				{
					if (dim == 1)
						ifv_R->delta_x = mv.X[cp[k][j]] - cv.X_c[cR];
					else if (dim == 2)
						{
							if(j == cp[k][0]) 
								{
									p_p = cp[k][1];
									p_n = cp[k][j];
								}				  
							else
								{
									p_p = cp[k][j+1];
									p_n = cp[k][j];
								}
							ifv_R->delta_x = 0.5*(mv.X[p_p] + mv.X[p_n]) - cv.X_c[cR];
							ifv_R->delta_y = 0.5*(mv.Y[p_p] + mv.Y[p_n]) - cv.Y_c[cR];
						}
					order2_i_f_var_init(cv, ifv_R, cR);						
				}
		}
	else if (cc[k][j] == -1)//initial boundary condition.		
		{
			ifv_R->U_rho = cv.U0_rho[k];
			ifv_R->U_e = cv.U0_e[k];
			ifv_R->U_u = cv.U0_u[k];
			if (dim > 1)
				ifv_R->U_v = cv.U0_v[k];
			if (dim > 2)
				ifv_R->U_w = cv.U0_w[k];
			if ((int)config[2] == 2)
				ifv_R->U_phi = cv.U0_phi[k];
			ifv_R->gamma = cv.gamma[k];

			if (order == 2)
				order2_i_f_var0(ifv_R);
		}
	else if (cc[k][j] == -2)//reflecting boundary condition.
		{
			if(cons2prim(ifv) == 0)
				{
					fprintf(stderr, "Error happens on primitive variable!\n");
					return -1;
				}
			ifv->F_rho = 0.0;
			ifv->F_u = ifv->P*ifv->n_x;
			ifv->F_v = ifv->P*ifv->n_y;
			ifv->F_e = 0.0;
			return 0;
		}
	else if (cc[k][j] == -3)//prescribed boundary condition.
		{
			ifv_R->U_rho = cv.U_rho[k];
			ifv_R->U_e = cv.U_e[k];
			ifv_R->U_u = cv.U_u[k];
			if (dim > 1)
				ifv_R->U_v = cv.U_v[k];
			if (dim > 2)
				ifv_R->U_w = cv.U_w[k];
			if ((int)config[2] == 2)
				ifv_R->U_phi = cv.U_phi[k];
			ifv_R->gamma = cv.gamma[k];

			if (order == 2)
				order2_i_f_var0(ifv_R);
		}		
	else
		{
			printf("No suitable boundary!\n");
			return -1;
		}

	if(cons2prim(ifv) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return -1;
		}
	if(cons2prim(ifv_R) == 0)
		{
			fprintf(stderr, "Error happens on primitive variable!\n");
			return -1;
		}

	return 1;
}

double tau_calc(const struct cell_var cv, const struct mesh_var mv)
{
	const double CFL = config[7];
	if (CFL < 0.0)
		return -CFL;
	const int dim = (int)config[0];
		
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	double tau = config[1];
	struct i_f_var ifv, ifv_R;
	double cum, lambda_max;
	int ivi;
	
	double qn, qn_R;
	double c, c_R;	
	
	for(int k = 0; k < num_cell; ++k)
		{
			cum = 0.0;
			
			for(int j = 1; j <= cp[k][0]; ++j)
				{
					ivi = interface_var_init(cv, mv, &ifv, &ifv_R, k, j);
					if (ivi == 0)
						;
					else if(ivi == -1)
						return -1.0;
					else if (dim == 2)
						{
							qn = ifv.U*ifv.n_x + ifv.V*ifv.n_y; 
							qn_R = ifv_R.U*ifv_R.n_x + ifv_R.V*ifv_R.n_y;
							c = sqrt(ifv.gamma * ifv.P / ifv.RHO);
							c_R = sqrt(ifv_R.gamma * ifv_R.P / ifv_R.RHO);
							lambda_max = fmax(c+fabs(qn), c_R+fabs(qn_R));					
							cum +=  lambda_max * ifv.length;
						}
				}		
			tau = fmin(tau, 2.0*cv.vol[k]/cum * CFL);
		}	//To decide tau.
	return tau;
}


void cons_qty_update(struct cell_var * cv, const struct mesh_var mv, const double tau)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	int ** cp = mv.cell_pt;
	
	int p_p, p_n;
	double length;
	for(int k = 0; k < num_cell; ++k)
		{
			for(int j = 1; j <= cp[k][0]; ++j)
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
					length = sqrt((mv.X[p_p] - mv.X[p_n])*(mv.X[p_p]-mv.X[p_n]) + (mv.Y[p_p] - mv.Y[p_n])*(mv.Y[p_p]-mv.Y[p_n]));
					
					cv->U_rho[k] += - tau*cv->F_rho[k][j] * length / cv->vol[k];
					cv->U_e[k]   += - tau*cv->F_e[k][j]   * length / cv->vol[k];
					cv->U_u[k]   += - tau*cv->F_u[k][j]   * length / cv->vol[k];
					if (dim > 1)
						cv->U_v[k] += - tau*cv->F_v[k][j] * length / cv->vol[k];
					if (dim > 2)
						cv->U_w[k] += - tau*cv->F_w[k][j] * length / cv->vol[k];	
					if ((int)config[2] == 2)
						cv->U_phi[k] += - tau*cv->F_phi[k][j] * length / cv->vol[k];
				}
		}
}
