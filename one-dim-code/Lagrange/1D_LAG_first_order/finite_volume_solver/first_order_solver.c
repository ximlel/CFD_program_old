#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <time.h>

#include "../finite_volume_solver.h"
#include "../../../lib/Riemann_solver.h"

#include <float.h>
#define min(x,y)  ( x<y?x:y )

/* This function use finite volume scheme to solve 1-D
 * equations of motion On Lagrange coordinate.
 *
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *m      is the number of the grids.
 */

//
int first_order_solver
(double * config, int m, 
 double * RHO[], double * U[], double * P[], double * E[], double * X[], double * cpu_time, char * scheme, double CFL)
{
  int j, k;  /* j is a frequently used index for
		      * spatial variables. k is a frequ-
		      * ently used index for the time
		      * step.
		      */

  clock_t tic, toc;
  double sum = 0.0;

  int const N = (int)(config[4]);  // the number of time steps
  double const eps = config[3];    // the largest value could be
                                   // seen as zero
  double tau;
  tau = DBL_MAX;    // the length of the time step
  double t_all = 0.0;
 
  double const gamma = config[0];      // the constant of the perfect gas
  int stop_step = 0;

  double const zeta = (gamma-1.0)/(gamma+1.0);

  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double c_L, c_R;               // the speed of sound
  double c[m];
  double u_star, p_star;
  double u_mid[m+1], p_mid[m+1]; // the Riemann solutions
  int CRW[2];


  double MASS[m];

  for(j = 0; j < m; ++j)
		  MASS[j] = (X[0][j+1]-X[0][j])*RHO[0][j];                                               
  for(j = 0; j < m; ++j)
		  E[0][j] = 0.5*U[0][j]*U[0][j] + P[0][j]/(gamma - 1.0)/RHO[0][j]; /* initialize the values of mass and energy.
*/

  double UL[N], PL[N], RHOL[N];
  double UR[N], PR[N], RHOR[N];


  
  for(k = 0; k < N; ++k)
  {
    UL[k] = U[0][0];
    UR[k] = U[0][m-1];
    PL[k] = P[0][0];
    PR[k] = P[0][m-1];
    RHOL[k] = RHO[0][0];
    RHOR[k] = RHO[0][m-1];
  }


//------------THE MAIN LOOP-------------
  for(k = 1; k <= N; ++k)
  {
    
    tic = clock();

    tau = 1.01*tau;

    for(j = 0; j < m; ++j)
    {
	c[j] = sqrt(gamma * P[k-1][j] / RHO[k-1][j]);
	tau = min(tau,CFL*(X[k-1][j+1]-X[k-1][j])/(c[j]+fabs(U[k-1][j])));	
    }
	
	if (CFL<0.0)
		{
			tau=-CFL;
		}
	t_all += tau;

	if(tau < config[1]/N/5)
		{
			printf("The length of the time step is so small at step %d, t_all=%lf, tau=%lf.\n",k,t_all,tau);
			stop_step=1;
		}
	
	if(t_all > config[1])
		{
			printf("The time is enough at step %d.\n",k);
			tau = tau - (t_all-config[1]);
			stop_step=1;
		}


    for(j = 0; j <= m; ++j)
    { /*
       *  j-1          j          j+1
       * j-1/2  j-1  j+1/2   j   j+3/2  j+1
       *   o-----X-----o-----X-----o-----X--...
       */

      if(j)
      { rho_L = RHO[k-1][j-1];
	u_L   =   U[k-1][j-1];
	p_L   =   P[k-1][j-1];
      }
      else
      { rho_L = RHOL[k-1];
	u_L   = UL[k-1];
	p_L   = PL[k-1];
      }

      if(j < m)
      { rho_R = RHO[k-1][j];
	u_R   =   U[k-1][j];
	p_R   =   P[k-1][j];
      }
      else
      { rho_R = RHOR[k-1];
	u_R = UR[k-1];
	p_R = PR[k-1];
      } /* initialize the values */

      c_L = sqrt(gamma * p_L / rho_L);
      c_R = sqrt(gamma * p_R / rho_R);

      if((p_L < eps) || (p_R < eps) || (rho_L < eps) || (rho_R < eps)||isnan(p_L)||isnan(p_R)||isnan(u_L)||isnan(u_R)||isnan(rho_L)||isnan(rho_R))
	{
		printf("Error firstly happens on step=%d, grid=%d.\n", k, j);
		stop_step=1;
		break;
	}
//===============Solve Riemann Problem==============
	
	if(strcmp(scheme,"Riemann_exact")==0)
		{
			Riemann_solver_exact(&u_star, &p_star, gamma, u_L, u_R, p_L, p_R, c_L, c_R, CRW, eps, 500);
			
			u_mid[j] = u_star;
			p_mid[j] = p_star;		
		}
	else
		{
			printf("No Riemann solver!");
			exit(7);
		}
	
//==================================================

	  X[k][j] = X[k-1][j] + tau*u_mid[j];//motion along the contact discontinuity
  
    }
	
//===============THE CORE ITERATION=================(On Lagrange Coordinate)

    for(j = 0; j < m; ++j)
    { /*
       *  j-1          j          j+1
       * j-1/2  j-1  j+1/2   j   j+3/2  j+1
       *   o-----X-----o-----X-----o-----X--...
       */
		RHO[k][j] = 1 / ( 1/RHO[k-1][j] + tau/MASS[j]*(u_mid[j+1] - u_mid[j]) );
		U[k][j] = U[k-1][j] - tau/MASS[j]*(p_mid[j+1] - p_mid[j]);
		E[k][j] = E[k-1][j] - tau/MASS[j]*(p_mid[j+1]*u_mid[j+1] - p_mid[j]*u_mid[j]);
		P[k][j] = (E[k][j] - 0.5*U[k][j]*U[k][j]) * (gamma - 1.0) * RHO[k][j];
    } //forward Euler


//==================================================

    toc = clock();
    cpu_time[k] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += cpu_time[k];


	if(stop_step)
		break;
  }

  printf("The cost of CPU time for this problem in 1D is %g seconds.\n", sum);

  if(!stop_step)
	printf("The maximum number of time steps is not enough for this calculation, t_tall=%lf.\n",t_all);
//------------END OF THE MAIN LOOP-------------

	return k;

}
