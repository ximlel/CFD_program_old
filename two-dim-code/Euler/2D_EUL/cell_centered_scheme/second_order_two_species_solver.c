#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#define min(x,y)  ( x<y?x:y )
#define max(x,y)  ( x>y?x:y )

#include <time.h>


#include "../../../lib/Riemann_solver.h"
#include "../../../lib/custom.h"
#include "../cell_centered_scheme.h"


/* This function use second order scheme to solve 2-D
 * Euler equations for two species by Eulerian method.
 *
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *NUM_CELL is the number of the grids.
 */


void second_order_two_species_solver
(int * STEP, double * E_M_wave, double * h, double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[], double * Z[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL/* the CFL number */)
{	

	int i, j, k, l; 

	clock_t tic, toc;
	double sum = 0.0;

	int const N = (int)(config[5]);  // the number of time steps
	double const eps = config[4];    // the largest value could be seen as zero
	double delta =0.2; // the coefficient of entropy fix
	double delta_God =0.2;
	double lambda_max; // the maximum eigenvalue
	
	double tau; // the length of the time step
	double t_all = 0.0;
	int stop_step = 0;

	double wave_speed[2], dire[4], mid[4];
	double rho_L, rho_R, u_L, u_R, qn_L, qn_R, v_L, v_R, qt_L, qt_R, p_L, p_R, c_L, c_R, z_L, z_R;
	double d_rho_L, d_rho_R, d_u_L, d_u_R, d_qn_L, d_qn_R, d_v_L, d_v_R, d_qt_L, d_qt_R, d_p_L, d_p_R, d_z_L, d_z_R;
	double u_mid, v_mid, rho_mid, p_mid, z_mid; // the Riemann solutions


	int CELL_RIGHT, STEP_RIGHT; //For boundary condition
	double cum;
	
	double u_con[NUM_CELL],  v_con[NUM_CELL],  e_con[NUM_CELL], z_con[NUM_CELL];

	for(k = 0; k < NUM_CELL; ++k)
		{			
			u_con[k] = RHO[0][k]*U[0][k];
			v_con[k] = RHO[0][k]*V[0][k];
			e_con[k] = P[0][k]/(gamma[k]-1.0) + 0.5*(U[0][k]*U[0][k]+V[0][k]*V[0][k])*RHO[0][k];
			z_con[k] = RHO[0][k]*Z[0][k];
		}	//Initialize conserved variables.

	
	int p_p,p_n;
  
	double VOLUME[NUM_CELL];
	for(k = 0; k < NUM_CELL; ++k)
		{
			VOLUME[k] = 0;
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					if(j == CELL_POINT[k][0]-1) 
					  {
						  p_p=CELL_POINT[k][1];
						  p_n=CELL_POINT[k][j+1];
					  }				  
					else
						{
							p_p=CELL_POINT[k][j+2];
							p_n=CELL_POINT[k][j+1];
						} 
					VOLUME[k] = VOLUME[k] + 0.5 * (X[p_n] * Y[p_p] - Y[p_n] * X[p_p]);
				}
		} //Calculate volume.
	
	
	double * n_x[NUM_CELL];
	initialize_memory(n_x,NUM_CELL,CELL_POINT);
	double * n_y[NUM_CELL];
	initialize_memory(n_y,NUM_CELL,CELL_POINT);
	
	
	int p_p_2,p_n_2;
	
	int * CELL_CELL[NUM_CELL];
	initialize_memory_int(CELL_CELL,NUM_CELL,CELL_POINT);
	int CELL_record;

  	for(k = 0; k < NUM_CELL; ++k)
		{ 						
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					if(j == CELL_POINT[k][0]-1) 
						{
							p_p=CELL_POINT[k][1];
							p_n=CELL_POINT[k][j+1];
						}				  
					else
						{
							p_p=CELL_POINT[k][j+2];
							p_n=CELL_POINT[k][j+1];
						} 
					
					CELL_record = 0;					
					
					for(i = 0; i < NUM_CELL; ++i)
						{
							if(i!=k)
								{
									for(l = 0; l < CELL_POINT[i][0]; ++l)
										{
											if(l == CELL_POINT[i][0]-1) 
												{
													p_p_2=CELL_POINT[i][1];
													p_n_2=CELL_POINT[i][l+1];
									         	}				  
											else
												{
													p_p_2=CELL_POINT[i][l+2];
													p_n_2=CELL_POINT[i][l+1];
												}
											
											if((p_p==p_n_2)&&(p_p_2==p_n))
												{
													CELL_CELL[k][j] = i;
													n_x[k][j] = (Y[p_p] - Y[p_n]) / sqrt((Y[p_p] - Y[p_n])*(Y[p_p] - Y[p_n])+(X[p_n] - X[p_p])*(X[p_n] - X[p_p]));
													n_y[k][j] = (X[p_n] - X[p_p]) / sqrt((Y[p_p] - Y[p_n])*(Y[p_p] - Y[p_n])+(X[p_n] - X[p_p])*(X[p_n] - X[p_p]));
													CELL_record = 1;
													break;
												}
										}
									
								}
							if(CELL_record)	break;								
						}
					
					if(CELL_record==0)
						{

							for(i = 0; i < NUM_BOUNDARY; ++i)
								{			
									if(i==(NUM_BOUNDARY-1))
										{
											p_p_2=BOUNDARY_POINT[0][0];
											p_n_2=BOUNDARY_POINT[0][i];
										} 
									else 
										{
											p_p_2=BOUNDARY_POINT[0][i+1];
											p_n_2=BOUNDARY_POINT[0][i];
										}

									if((p_p==p_p_2)&&(p_n==p_n_2))
										{
											CELL_CELL[k][j] = BOUNDARY_POINT[1][i];
											n_x[k][j] = (Y[p_p] - Y[p_n]) / sqrt((Y[p_p] - Y[p_n])*(Y[p_p] - Y[p_n])+(X[p_n] - X[p_p])*(X[p_n] - X[p_p]));
											n_y[k][j] = (X[p_n] - X[p_p]) / sqrt((Y[p_p] - Y[p_n])*(Y[p_p] - Y[p_n])+(X[p_n] - X[p_p])*(X[p_n] - X[p_p]));
											CELL_record = 1;
											break;
										}
								}							
						}
					if(CELL_record==0)
						{
							printf("Some errors happen when meshing.\n");
							exit(6);
						}
					
				}
		} //Determine the normal direction and relationship between cells.


	printf("Grid has been constructed.\n");
	
	double F_mk[5];

	double * F_mk_1[NUM_CELL];
	initialize_memory(F_mk_1,NUM_CELL,CELL_POINT);
	double * F_mk_2[NUM_CELL];
	initialize_memory(F_mk_2,NUM_CELL,CELL_POINT);
	double * F_mk_3[NUM_CELL];
	initialize_memory(F_mk_3,NUM_CELL,CELL_POINT);
	double * F_mk_4[NUM_CELL];
	initialize_memory(F_mk_4,NUM_CELL,CELL_POINT);
	double * F_mk_5[NUM_CELL];
	initialize_memory(F_mk_5,NUM_CELL,CELL_POINT);



	double S, S_tri;
	double X_c[NUM_CELL], Y_c[NUM_CELL];
	double delta_X, delta_Y;

	double grad_RHO_x[NUM_CELL], grad_RHO_y[NUM_CELL];
	double grad_P_x[NUM_CELL], grad_P_y[NUM_CELL];
	double grad_U_x[NUM_CELL], grad_U_y[NUM_CELL];
	double grad_V_x[NUM_CELL], grad_V_y[NUM_CELL];
	double grad_Z_x[NUM_CELL], grad_Z_y[NUM_CELL];




//------------THE MAIN LOOP-------------

	for(i = 0; i < * STEP; ++i)
		{	

			tic = clock();		
			tau = DBL_MAX; 

			for(k = 0; k < NUM_CELL; ++k)
				{
					S = 0.0;
					X_c[k] = 0.0;
					Y_c[k] = 0.0;
					for(j = 2; j < CELL_POINT[k][0]; ++j)
						{
							S_tri =X[CELL_POINT[k][1]]*Y[CELL_POINT[k][j]] + X[CELL_POINT[k][j+1]]*Y[CELL_POINT[k][1]] + X[CELL_POINT[k][j]]*Y[CELL_POINT[k][j+1]] - X[CELL_POINT[k][j+1]]*Y[CELL_POINT[k][j]] - X[CELL_POINT[k][1]]*Y[CELL_POINT[k][j+1]] - X[CELL_POINT[k][j]]*Y[CELL_POINT[k][1]];
							X_c[k] = X_c[k] + (X[CELL_POINT[k][1]] + X[CELL_POINT[k][j]] + X[CELL_POINT[k][j+1]]) * S_tri;
							Y_c[k] = Y_c[k] + (Y[CELL_POINT[k][1]] + Y[CELL_POINT[k][j]] + Y[CELL_POINT[k][j+1]]) * S_tri;
							S = S + S_tri;
						}			 
					X_c[k] = X_c[k]/S/3;
					Y_c[k] = Y_c[k]/S/3;
				}

			slope_limiter_Ven(X_c, Y_c, X, Y, grad_RHO_x, grad_RHO_y, RHO, NUM_CELL, config, CELL_CELL, CELL_POINT, m, n);
			slope_limiter_Ven(X_c, Y_c, X, Y, grad_U_x, grad_U_y, U, NUM_CELL, config, CELL_CELL, CELL_POINT, m ,n);
			slope_limiter_Ven(X_c, Y_c, X, Y, grad_V_x, grad_V_y,  V, NUM_CELL, config, CELL_CELL, CELL_POINT, m, n);
			slope_limiter_Ven(X_c, Y_c, X, Y, grad_P_x, grad_P_y, P, NUM_CELL, config, CELL_CELL, CELL_POINT, m, n);
			slope_limiter_Ven(X_c, Y_c, X, Y, grad_Z_x, grad_Z_y, Z, NUM_CELL, config, CELL_CELL, CELL_POINT, m, n);

			
			for(k = 0; k < NUM_CELL; ++k)
				{
					cum = 0.0;
					
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							if(j == CELL_POINT[k][0]-1) 
								{
									p_p=CELL_POINT[k][1];
									p_n=CELL_POINT[k][j+1];
								}				  
							else
								{
									p_p=CELL_POINT[k][j+2];
									p_n=CELL_POINT[k][j+1];
								}
															
							if (CELL_CELL[k][j]==-2)//reflecting boundary condition.
								{
									lambda_max = 0.0;
								}
							else
								{
									if (CELL_CELL[k][j]>=0)
										{											
											STEP_RIGHT = 1;							
											CELL_RIGHT = CELL_CELL[k][j];
										}
									else if (CELL_CELL[k][j]==-1)//initial boundary condition.
										{
											STEP_RIGHT = 0;
											CELL_RIGHT = k;
										}
									else if (CELL_CELL[k][j]==-3)//prescribed boundary condition.
										{
											STEP_RIGHT = 1;
											CELL_RIGHT = k;
										}
									else if (CELL_CELL[k][j]==-4)//periodic boundary condition in x-direction.
										{
											STEP_RIGHT = 1;
											if(!(k%n))
												CELL_RIGHT = k+n-1;
											else if(k%n==n-1)
												CELL_RIGHT = k-n+1;
											else
												{																									
													printf("Something wrong as we construct periodic boundary condition in x-direction.\n");
													exit(8);
												}
										}
									else
										{
											printf("No suitable boundary!\n");
											exit(7);
										}						 
							   																									
											qn_L = U[1][k]*n_x[k][j] + V[1][k]*n_y[k][j]; 
											qn_R = U[STEP_RIGHT][CELL_RIGHT]*n_x[k][j] + V[STEP_RIGHT][CELL_RIGHT]*n_y[k][j];
											c_L = sqrt(gamma[k] * P[1][k] / RHO[1][k]);
											c_R = sqrt(gamma[k] * P[STEP_RIGHT][CELL_RIGHT] / RHO[STEP_RIGHT][CELL_RIGHT]);
											lambda_max = max(c_L+fabs(qn_L),c_R+fabs(qn_R));					
								}								
							
							cum +=  lambda_max*sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n]));						
						}
					
					tau = min(tau, 2.0*VOLUME[k]/cum*CFL);
				}	//To decide tau.
			
			if (CFL<0.0)
				{
					tau=-CFL;
				}
			t_all += tau;
			
			if(tau < config[1]/N/5)
				{
					printf("The length of the time step is so small at step %d, t_all=%lf, tau=%lf.\n",i,t_all,tau);
					stop_step=1;
				}
			
			if(t_all > config[1])
				{
					printf("The time is enough at step %d.\n",i);
					tau = tau - (t_all-config[1]);
					stop_step=1;
				} // Time



			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							if(j == CELL_POINT[k][0]-1) 
								{
									p_p=CELL_POINT[k][1];
									p_n=CELL_POINT[k][j+1];
								}				  
							else
								{
									p_p=CELL_POINT[k][j+2];
									p_n=CELL_POINT[k][j+1];
								}
											delta_X = 0.5*(X[p_p]+X[p_n]) - X_c[k];
											delta_Y = 0.5*(Y[p_p]+Y[p_n]) - Y_c[k];
											
							if (CELL_CELL[k][j]==-2)//reflecting boundary condition.
								{
									F_mk[0] = 0.0;
									F_mk[1] = P[1][k]*n_x[k][j];
									F_mk[2] = P[1][k]*n_y[k][j];
									F_mk[3] = 0.0;
									F_mk[4] = 0.0;
								}
							else
								{
									if (CELL_CELL[k][j]>=0)
										{											
											STEP_RIGHT = 1;							
											CELL_RIGHT = CELL_CELL[k][j];
											
											d_rho_R = grad_RHO_x[CELL_RIGHT]*n_x[k][j] + grad_RHO_y[CELL_RIGHT]*n_y[k][j];
											d_u_R = grad_U_x[CELL_RIGHT]*n_x[k][j] + grad_U_y[CELL_RIGHT]*n_y[k][j];
											d_v_R = grad_V_x[CELL_RIGHT]*n_x[k][j] + grad_V_y[CELL_RIGHT]*n_y[k][j];
											d_p_R = grad_P_x[CELL_RIGHT]*n_x[k][j] + grad_P_y[CELL_RIGHT]*n_y[k][j];
											d_z_R = grad_Z_x[CELL_RIGHT]*n_x[k][j] + grad_Z_y[CELL_RIGHT]*n_y[k][j];

											rho_R = RHO[STEP_RIGHT][CELL_RIGHT] - grad_RHO_x[CELL_RIGHT] * delta_X - grad_RHO_y[CELL_RIGHT] * delta_Y;
											u_R = U[STEP_RIGHT][CELL_RIGHT] - grad_U_x[CELL_RIGHT] * delta_X - grad_U_y[CELL_RIGHT] * delta_Y;
											v_R = V[STEP_RIGHT][CELL_RIGHT] - grad_V_x[CELL_RIGHT] * delta_X - grad_V_y[CELL_RIGHT] * delta_Y;
											p_R = P[STEP_RIGHT][CELL_RIGHT] - grad_P_x[CELL_RIGHT] * delta_X - grad_P_y[CELL_RIGHT] * delta_Y;
											z_R = Z[STEP_RIGHT][CELL_RIGHT] - grad_Z_x[CELL_RIGHT] * delta_X - grad_Z_y[CELL_RIGHT] * delta_Y;
										}
									else if (CELL_CELL[k][j]==-1)//initial boundary condition.
										{
											STEP_RIGHT = 0;
											CELL_RIGHT = k;
											
											d_rho_R = 0.0;
											d_u_R = 0.0;
											d_v_R = 0.0;
											d_p_R = 0.0;
											d_z_R = 0.0;

											rho_R = RHO[STEP_RIGHT][CELL_RIGHT];
											u_R = U[STEP_RIGHT][CELL_RIGHT];
											v_R = V[STEP_RIGHT][CELL_RIGHT];
											p_R = P[STEP_RIGHT][CELL_RIGHT];
											z_R = Z[STEP_RIGHT][CELL_RIGHT];
										}
									else if (CELL_CELL[k][j]==-3)//prescribed boundary condition.
										{
											STEP_RIGHT = 1;
											CELL_RIGHT = k;

											d_rho_R = grad_RHO_x[CELL_RIGHT]*n_x[k][j] + grad_RHO_y[CELL_RIGHT]*n_y[k][j];
											d_u_R = grad_U_x[CELL_RIGHT]*n_x[k][j] + grad_U_y[CELL_RIGHT]*n_y[k][j];
											d_v_R = grad_V_x[CELL_RIGHT]*n_x[k][j] + grad_V_y[CELL_RIGHT]*n_y[k][j];
											d_p_R = grad_P_x[CELL_RIGHT]*n_x[k][j] + grad_P_y[CELL_RIGHT]*n_y[k][j];
											d_z_R = grad_Z_x[CELL_RIGHT]*n_x[k][j] + grad_Z_y[CELL_RIGHT]*n_y[k][j];
											
											rho_R = RHO[STEP_RIGHT][CELL_RIGHT];
											u_R = U[STEP_RIGHT][CELL_RIGHT];
											v_R = V[STEP_RIGHT][CELL_RIGHT];
											p_R = P[STEP_RIGHT][CELL_RIGHT];
											z_R = Z[STEP_RIGHT][CELL_RIGHT];
										}
									else if (CELL_CELL[k][j]==-4)//periodic boundary condition in x-direction.
										{
											STEP_RIGHT = 1;
											if(!(k%n))
												CELL_RIGHT = k+n-1;
											else if(k%n==n-1)
												CELL_RIGHT = k-n+1;
											else
												{																									
													printf("Something wrong as we construct periodic boundary condition in x-direction.\n");
													exit(8);
												}

											d_rho_R = grad_RHO_x[CELL_RIGHT]*n_x[k][j] + grad_RHO_y[CELL_RIGHT]*n_y[k][j];
											d_u_R = grad_U_x[CELL_RIGHT]*n_x[k][j] + grad_U_y[CELL_RIGHT]*n_y[k][j];
											d_v_R = grad_V_x[CELL_RIGHT]*n_x[k][j] + grad_V_y[CELL_RIGHT]*n_y[k][j];
											d_p_R = grad_P_x[CELL_RIGHT]*n_x[k][j] + grad_P_y[CELL_RIGHT]*n_y[k][j];
											d_z_R = grad_Z_x[CELL_RIGHT]*n_x[k][j] + grad_Z_y[CELL_RIGHT]*n_y[k][j];

											rho_R = RHO[STEP_RIGHT][CELL_RIGHT] - grad_RHO_x[CELL_RIGHT] * delta_X - grad_RHO_y[CELL_RIGHT] * delta_Y;
											u_R = U[STEP_RIGHT][CELL_RIGHT] - grad_U_x[CELL_RIGHT] * delta_X - grad_U_y[CELL_RIGHT] * delta_Y;
											v_R = V[STEP_RIGHT][CELL_RIGHT] - grad_V_x[CELL_RIGHT] * delta_X - grad_V_y[CELL_RIGHT] * delta_Y;
											p_R = P[STEP_RIGHT][CELL_RIGHT] - grad_P_x[CELL_RIGHT] * delta_X - grad_P_y[CELL_RIGHT] * delta_Y;
											z_R = Z[STEP_RIGHT][CELL_RIGHT] - grad_Z_x[CELL_RIGHT] * delta_X - grad_Z_y[CELL_RIGHT] * delta_Y;	
										}
									else
										{
											printf("No suitable boundary!\n");
											exit(7);
										}						 								
								
									//============================================the GRP scheme with exact Riemann solver==========================================

									if(strcmp(scheme,"GRP")==0)
										{
											d_rho_L =  grad_RHO_x[k]*n_x[k][j] + grad_RHO_y[k]*n_y[k][j];
											d_u_L =  grad_U_x[k]*n_x[k][j] + grad_U_y[k]*n_y[k][j];
											d_v_L =  grad_V_x[k]*n_x[k][j] + grad_V_y[k]*n_y[k][j];
											d_p_L =  grad_P_x[k]*n_x[k][j] + grad_P_y[k]*n_y[k][j];
											d_z_L =  grad_Z_x[k]*n_x[k][j] + grad_Z_y[k]*n_y[k][j];

											rho_L = RHO[1][k] + grad_RHO_x[k] * delta_X + grad_RHO_y[k] * delta_Y;
											u_L = U[1][k] + grad_U_x[k] * delta_X + grad_U_y[k] * delta_Y;
											v_L = V[1][k] + grad_V_x[k] * delta_X + grad_V_y[k] * delta_Y;
											p_L = P[1][k] + grad_P_x[k] * delta_X + grad_P_y[k] * delta_Y;
											z_L = Z[1][k] + grad_Z_x[k] * delta_X + grad_Z_y[k] * delta_Y;

											qn_L = u_L*n_x[k][j] + v_L*n_y[k][j];
											qt_L = -u_L*n_y[k][j] + v_L*n_x[k][j];
											d_qn_L = d_u_L*n_x[k][j] + d_v_L*n_y[k][j];
											d_qt_L = -d_u_L*n_y[k][j] + d_v_L*n_x[k][j];
											qn_R = u_R*n_x[k][j] + v_R*n_y[k][j];
											qt_R = -u_R*n_y[k][j] + v_R*n_x[k][j];
											d_qn_R = d_u_R*n_x[k][j] + d_v_R*n_y[k][j];
											d_qt_R = -d_u_R*n_y[k][j] + d_v_R*n_x[k][j];				   						

											linear_GRP_solver_Edir_2D(wave_speed, dire, mid, 0.0, rho_L,rho_R, d_rho_L, d_rho_R, qn_L, qn_R, d_qn_L, d_qn_R, qt_L, qt_R, d_qt_L, d_qt_R, p_L, p_R, d_p_L, d_p_R, gamma[k], eps);

											rho_mid = mid[0] + 0.5*tau*dire[0];
											p_mid = mid[3] + 0.5*tau*dire[3];
											u_mid = (mid[1] + 0.5*tau*dire[1])*n_x[k][j] - (mid[2] + 0.5*tau*dire[2])*n_y[k][j];
											v_mid = (mid[1] + 0.5*tau*dire[1])*n_y[k][j] + (mid[2] + 0.5*tau*dire[2])*n_x[k][j];

											F_mk[0] = rho_mid*u_mid*n_x[k][j] + rho_mid*v_mid*n_y[k][j];
											F_mk[1] = F_mk[0]*u_mid + p_mid*n_x[k][j];
											F_mk[2] = F_mk[0]*v_mid + p_mid*n_y[k][j];
											F_mk[3] = (gamma[k]/(gamma[k]-1.0))*p_mid/rho_mid + 0.5*u_mid*u_mid + 0.5*v_mid*v_mid;
											F_mk[3] = F_mk[0]*F_mk[3];

											linear_GRP_solver_Edir_2D(wave_speed, dire, mid, 0.0, rho_L,rho_R, d_rho_L, d_rho_R, qn_L, qn_R, d_qn_L, d_qn_R, z_L, z_R, d_z_L, d_z_R, p_L, p_R, d_p_L, d_p_R, gamma[k], eps);

											z_mid = mid[2] + 0.5*tau*dire[2];
											
											F_mk[4] = F_mk[0]*z_mid;											
										}								
									else
										{
											printf("No Riemann solver!\n");
											exit(7);
										}
								}							
							F_mk_1[k][j] = F_mk[0];
							F_mk_2[k][j] = F_mk[1];
							F_mk_3[k][j] = F_mk[2];
							F_mk_4[k][j] = F_mk[3];
							F_mk_5[k][j] = F_mk[4];
						}					
				}	//Solver



	for(k = 0; k < NUM_CELL; ++k)
		{
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					if(j == CELL_POINT[k][0]-1) 
						{
							p_p=CELL_POINT[k][1];
							p_n=CELL_POINT[k][j+1];
						}				  
					else
						{
							p_p=CELL_POINT[k][j+2];
							p_n=CELL_POINT[k][j+1];
						}
												
					RHO[1][k] += - tau*F_mk_1[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					u_con[k] += - tau*F_mk_2[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					v_con[k] += - tau*F_mk_3[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					e_con[k] += - tau*F_mk_4[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					z_con[k] += - tau*F_mk_5[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
				}	
	
			U[1][k] = u_con[k]/RHO[1][k];
			V[1][k] = v_con[k]/RHO[1][k];
			P[1][k] = (e_con[k] - 0.5*(U[1][k]*U[1][k]+V[1][k]*V[1][k])*RHO[1][k])*(gamma[k]-1.0);
			Z[1][k] = z_con[k]/RHO[1][k];			
			if((RHO[1][k] < eps) || (P[1][k] < eps)||/* (Z[1][k] < -1.0*eps)|| (Z[1][k] > 1.0+eps) ||*/ isnan(RHO[1][k])||isnan(U[1][k])||isnan(V[1][k])||isnan(P[1][k])||isnan(Z[1][k]))
				{
					if(!stop_step)
						printf("Error firstly happens at t_all=%lf, step=%d, on cell=%d", t_all, i, k);
					else
						printf (",%d",k);
					stop_step=2;
				}
		}
	if(stop_step==2)
		printf(".\n");	


	
//==================================================

    toc = clock();
    cpu_time[i] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += cpu_time[i];

	if(stop_step)
		break;
		
		}



//=======================kinetic energy of mixing region=====================
	
	double phi[n];
	double const phi_0 = 1;
	double const phi_m = 0.25;
	
	for(k = 0; k < m; ++k)
		{
			phi[k] = 0;
			for(j = 0; j < n ;++j)
				{									
					phi[k] += Z[1][m*j + k];
				}
			phi[k] = phi[k]/n;
		}
	
	* h = 0.0;
	for(k = 0; k < m; ++k)
		* h += phi[k]*(phi_0-phi[k])/phi_m/phi_m*config[3];

	
	double rho_bar, U_bar, V_bar, U_wave, V_wave, VOL_MIX;

	for(k = 0; k < m; ++k)
		{
			if(phi[k]>0.05&&phi[k]<0.95) 							   
				for(j = 0; j < n ;++j)
					{
						VOL_MIX += VOLUME[m*j+k];
						rho_bar += RHO[1][m*j+k]*VOLUME[m*j+k];
						U_bar += U[1][m*j+k]*VOLUME[m*j+k];
						V_bar += V[1][m*j+k]*VOLUME[m*j+k];
						U_wave += RHO[1][m*j+k]*U[1][m*j+k]*VOLUME[m*j+k];
						V_wave += RHO[1][m*j+k]*V[1][m*j+k]*VOLUME[m*j+k];
					}
		}
	U_bar = U_bar/VOL_MIX;
	V_bar = V_bar/VOL_MIX;
	U_wave = U_wave/rho_bar;
	V_wave = V_wave/rho_bar;


	double U_M, V_M, E_M;
	
	for(k = 0; k < m; ++k)
		{
			if(phi[k]>0.05&&phi[k]<0.95) 							   
				for(j = 0; j < n ;++j)
					{
						U_M = U[1][m*j+k];
						V_M = V[1][m*j+k];
						E_M += 0.5*RHO[1][m*j+k]*(U_M*U_M+V_M*V_M);
					}
		}
	E_M_wave[0] = E_M/rho_bar;

	for(k = 0; k < m; ++k)
		{
			if(phi[k]>0.05&&phi[k]<0.95) 							   
				for(j = 0; j < n ;++j)
					{
						U_M = U[1][m*j+k] - U_bar;
						V_M = V[1][m*j+k] - V_bar;
						E_M += 0.5*RHO[1][m*j+k]*(U_M*U_M+V_M*V_M);
					}
		}
	E_M_wave[1] = E_M/rho_bar;

	for(k = 0; k < m; ++k)
		{
			if(phi[k]>0.05&&phi[k]<0.95) 							   
				for(j = 0; j < n ;++j)
					{
						U_M = U[1][m*j+k] - U_wave;
						V_M = V[1][m*j+k] - V_wave;
						E_M += 0.5*RHO[1][m*j+k]*(U_M*U_M+V_M*V_M);
					}
		}
	E_M_wave[2] = E_M/rho_bar;

//========================================================	


	
	printf("The cost of CPU time for 2D equations of motion by Eulerian method is %g seconds.\n", sum);

	if(!stop_step)
		printf("The maximum number of time steps is not enough for this calculation, t_all=%lf.\n",t_all);

//------------END OF THE MAIN LOOP-------------

  for(k = 0; k < NUM_CELL; ++k)
	  { 	
		  free(CELL_CELL[k]);
		  free(n_x[k]);
		  free(n_y[k]);
		  free(F_mk_1[k]);
		  free(F_mk_2[k]);
		  free(F_mk_3[k]);
		  free(F_mk_4[k]);
		  free(F_mk_5[k]);
		  CELL_CELL[k] = NULL;
		  n_x[k] = NULL;
		  n_y[k] = NULL;
		  F_mk_1[k] = NULL;
		  F_mk_2[k] = NULL;
		  F_mk_3[k] = NULL;
		  F_mk_4[k] = NULL;
		  F_mk_5[k] = NULL;
	  }

  * STEP = i;  
}
