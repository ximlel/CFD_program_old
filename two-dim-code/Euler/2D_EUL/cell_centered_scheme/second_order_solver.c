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
 * equations of motion by Eulerian method.
 *
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *NUM_CELL is the number of the grids.
 */


int second_order_solver
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[],
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

	double dire[3], mid[3];
	double mid_qt;
	double u_L, u_R, c_L, c_R;
	double u_mid, v_mid, rho_mid, p_mid; // the Riemann solutions

	int CELL_RIGHT, STEP_RIGHT; //For boundary condition
	double cum;
	
	double u_con[NUM_CELL],  v_con[NUM_CELL],  e_con[NUM_CELL];
	for(k = 0; k < NUM_CELL; ++k)
		{			
			u_con[k] = RHO[0][k]*U[0][k];
			v_con[k] = RHO[0][k]*V[0][k];
			e_con[k] = P[0][k]/(gamma[k]-1.0) + 0.5*(U[0][k]*U[0][k]+V[0][k]*V[0][k])*RHO[0][k];				
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
	
	double F_mk[4];

	double * F_mk_1[NUM_CELL];
	initialize_memory(F_mk_1,NUM_CELL,CELL_POINT);
	double * F_mk_2[NUM_CELL];
	initialize_memory(F_mk_2,NUM_CELL,CELL_POINT);
	double * F_mk_3[NUM_CELL];
	initialize_memory(F_mk_3,NUM_CELL,CELL_POINT);
	double * F_mk_4[NUM_CELL];
	initialize_memory(F_mk_4,NUM_CELL,CELL_POINT);



	double S, S_tri;
	double X_c[NUM_CELL], Y_c[NUM_CELL];

	double grad_RHO_x[NUM_CELL], grad_RHO_y[NUM_CELL];
	double grad_P_x[NUM_CELL], grad_P_y[NUM_CELL];
	double grad_U_x[NUM_CELL], grad_U_y[NUM_CELL];
	double grad_V_x[NUM_CELL], grad_V_y[NUM_CELL];
	double grad_CC_x[NUM_CELL], grad_CC_y[NUM_CELL];


//------------THE MAIN LOOP-------------

	for(i = 0; i < N; ++i)
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

			slop_limiter_Ven(X_c, Y_c, X, Y, grad_RHO_x, grad_RHO_y, RHO, NUM_CELL, config, CELL_CELL, CELL_POINT, i, m, n);
			slop_limiter_Ven(X_c, Y_c, X, Y, grad_U_x, grad_U_y, U, NUM_CELL, config, CELL_CELL, CELL_POINT, i, m ,n);
			slop_limiter_Ven(X_c, Y_c, X, Y, grad_V_x, grad_V_y,  V, NUM_CELL, config, CELL_CELL, CELL_POINT, i, m, n);
			slop_limiter_Ven(X_c, Y_c, X, Y, grad_P_x, grad_P_y, P, NUM_CELL, config, CELL_CELL, CELL_POINT, i, m, n);
			
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
									F_mk[0] = 0.0;
									F_mk[1] = P[i][k]*n_x[k][j];
									F_mk[2] = P[i][k]*n_y[k][j];
									F_mk[3] = 0.0;
									lambda_max = 0.0;
								}
							else
								{
									if (CELL_CELL[k][j]>=0)
										{											
											STEP_RIGHT = i;							
											CELL_RIGHT = CELL_CELL[k][j];
										}
									else if (CELL_CELL[k][j]==-1)//initial boundary condition.
										{
											STEP_RIGHT = 0;
											CELL_RIGHT = k;
										}
									else if (CELL_CELL[k][j]==-3)//prescribed boundary condition.
										{
											STEP_RIGHT = i;
											CELL_RIGHT = k;
										}
									else if (CELL_CELL[k][j]==-4)//periodic boundary condition in x-direction.
										{
											STEP_RIGHT = i;
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
								
								
//============================================the Goundov scheme of exact Riemann solver==========================================

									if(strcmp(scheme,"GRP")==0)
										{																				
											u_L = U[i][k]*n_x[k][j] + V[i][k]*n_y[k][j]; 
											u_R = U[STEP_RIGHT][CELL_RIGHT]*n_x[k][j] + V[STEP_RIGHT][CELL_RIGHT]*n_y[k][j];
											c_L = sqrt(gamma[k] * P[i][k] / RHO[i][k]);
											c_R = sqrt(gamma[k] * P[STEP_RIGHT][CELL_RIGHT] / RHO[STEP_RIGHT][CELL_RIGHT]);
											linear_GRP_solver_Edir(dire, mid, RHO[i][k],RHO[STEP_RIGHT][CELL_RIGHT], 0.0, 0.0, u_L, u_R, 0.0, 0.0, P[i][k], P[STEP_RIGHT][CELL_RIGHT], 0.0, 0.0, gamma[k], eps);
											rho_mid = mid[0];
											p_mid = mid[2];

											/*if(fabs(mid[1])<delta_God)
											  mid_qt = 0.5*(-U[i][k]*n_y[k][j] + V[i][k]*n_x[k][j])*(mid[1]/delta_God+1) + 0.5*(-U[STEP_RIGHT][CELL_RIGHT]*n_y[k][j] + V[STEP_RIGHT][CELL_RIGHT]*n_x[k][j])*(-mid[1]/delta_God+1);
											  else*/ if(mid[1]>0)
												mid_qt = -U[i][k]*n_y[k][j] + V[i][k]*n_x[k][j];
											else
												mid_qt = -U[STEP_RIGHT][CELL_RIGHT]*n_y[k][j] + V[STEP_RIGHT][CELL_RIGHT]*n_x[k][j];
											u_mid = mid[1]*n_x[k][j] - mid_qt*n_y[k][j];
											v_mid = mid[1]*n_y[k][j] + mid_qt*n_x[k][j];
											
											F_mk[0] = rho_mid*u_mid*n_x[k][j] + rho_mid*v_mid*n_y[k][j];
											F_mk[1] = F_mk[0]*u_mid + p_mid*n_x[k][j];
											F_mk[2] = F_mk[0]*v_mid + p_mid*n_y[k][j];
											F_mk[3] = (gamma[k]/(gamma[k]-1.0))*p_mid/rho_mid + 0.5*u_mid*u_mid + 0.5*v_mid*v_mid;
											F_mk[3] = F_mk[0]*F_mk[3];
											lambda_max = max(c_L+fabs(u_L),c_R+fabs(u_R));					
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
							
							cum +=  lambda_max*sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n]));
							
						}
					
					tau = min(tau, 2.0*VOLUME[k]/cum*CFL);
				}	//Solver
			
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
			RHO[i+1][k] = RHO[i][k];
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
												
					RHO[i+1][k] += - tau*F_mk_1[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					u_con[k] += - tau*F_mk_2[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					v_con[k] += - tau*F_mk_3[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
					e_con[k] += - tau*F_mk_4[k][j] * sqrt((X[p_p]-X[p_n])*(X[p_p]-X[p_n])+(Y[p_p]-Y[p_n])*(Y[p_p]-Y[p_n])) / VOLUME[k];
				}
			
			U[i+1][k] = u_con[k]/RHO[i+1][k];
			V[i+1][k] = v_con[k]/RHO[i+1][k];
			P[i+1][k] = (e_con[k] - 0.5*(U[i+1][k]*U[i+1][k]+V[i+1][k]*V[i+1][k])*RHO[i+1][k])*(gamma[k]-1.0);
			if(P[i+1][k] < eps)
				{
					printf ("P is smaller than 0, error firstly happens in cell %d and step %d, t_all=%lf.\n",k,i,t_all);
					stop_step=1;
				}
		}
	


	
//==================================================

    toc = clock();
    cpu_time[i] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += cpu_time[i];

	if(stop_step)
		break;
		
				}

	printf("The cost of CPU time for 2D equations of motion by Eulerian method is %g seconds.\n", sum);

	if(!stop_step)
		printf("The maximum number of time steps is not enough for this calculation, t_tall=%lf.\n",t_all);

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
		  CELL_CELL[k] = NULL;
		  n_x[k] = NULL;
		  n_y[k] = NULL;
		  F_mk_1[k] = NULL;
		  F_mk_2[k] = NULL;
		  F_mk_3[k] = NULL;
		  F_mk_4[k] = NULL;
	  }

  return i+1;

  
}
