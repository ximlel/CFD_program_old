#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <time.h>

#include "../finite_volume_solver.h"
#include "../../../lib/Riemann_solver.h"

#include <float.h>
#define min(x,y)  ( x<y?x:y )

/* This function use GRP scheme to solve 1-D
 * equations of motion On Lagrange coordinate.
 *
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *m      is the number of the grids.
 */

//
int second_order_solver
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
  double const h = config[2];      // the length of the initial spatial grids
  double tau;
  tau = DBL_MAX;    // the length of the time step
  double t_all = 0.0;
 
  double const gamma = config[0];      // the constant of the perfect gas
  int stop_step = 0;

  double const zeta = (gamma-1.0)/(gamma+1.0);


  double s_rho[m], s_u[m], s_p[m];
  double s_L, s_R;
  double u_L, p_L, rho_L;
  double u_R, p_R, rho_R;
  double t_u_L, t_p_L, t_rho_L;
  double t_u_R, t_p_R, t_rho_R;
  double dire[4], mid[4];  //RHO_L_t,U_t,P_t,RHO_R_t.
  double s_u_L, s_p_L, s_rho_L;
  double s_u_R, s_p_R, s_rho_R;
  double slope_temp;
  double X_mass[m];
  double c[m]; // the speed of sound

  double U_next[m+1], P_next[m+1], RHO_next_L[m+1], RHO_next_R[m+1];
  
  double U_F[m+1], P_F[m+1];// the numerical flux at t_{n+1/2}

  double s_time = 1.9; // the paramater in slope limiters



  double MASS[m];

  for(j = 0; j < m; ++j)
		  MASS[j] = (X[0][j+1]-X[0][j])*RHO[0][j];   
  for(j = 0; j < m; ++j)
		  E[0][j] = 0.5*U[0][j]*U[0][j] + P[0][j]/(gamma - 1.0)/RHO[0][j]; /* initialize the values of mass and energy.
*/

  double UL[N], PL[N], RHOL[N];
  double UR[N], PR[N], RHOR[N];
  double SUL[N], SPL[N], SRHOL[N];
  double SUR[N], SPR[N], SRHOR[N];

  
  for(k = 0; k < N; ++k)
  {
    UL[k] = U[0][0];
    UR[k] = U[0][m-1];
    PL[k] = P[0][0];
    PR[k] = P[0][m-1];
    RHOL[k] = RHO[0][0];
    RHOR[k] = RHO[0][m-1];
    SUL[k] = 0.0;
    SUR[k] = 0.0;
    SPL[k] = 0.0;
    SPR[k] = 0.0;
    SRHOL[k] = 0.0;
    SRHOR[k] = 0.0;
  }


  for(j = 0; j < m; ++j)
  {
    if(j)
    {
      s_u_L = (U[0][j] - U[0][j-1]) / h;
      s_p_L = (P[0][j] - P[0][j-1]) / h;
      s_rho_L = (RHO[0][j] - RHO[0][j-1]) / h;
    }
    else
    {
      s_u_L = (U[0][j] - UL[0]) / h;
      s_p_L = (P[0][j] - PL[0]) / h;
      s_rho_L = (RHO[0][j] - RHOL[0]) / h;
    }
    if(j < m-1)
    {
      s_u_R = (U[0][j+1] - U[0][j]) / h;
      s_p_R = (P[0][j+1] - P[0][j]) / h;
      s_rho_R = (RHO[0][j+1] - RHO[0][j]) / h;
    }
    else
    {
      s_u_R = (UR[0] - U[0][j]) / h;
      s_p_R = (PR[0] - P[0][j]) / h;
      s_rho_R = (RHOR[0] - RHO[0][j]) / h;
    }

    if(s_u_L * s_u_R < 0.0)
      s_u[j] = 0.0;
    else if(fabs(s_u_R) < fabs(s_u_L))
      s_u[j] = s_u_R;
    else
      s_u[j] = s_u_L;

    if(s_p_L * s_p_R < 0.0)
      s_p[j] = 0.0;
    else if(fabs(s_p_R) < fabs(s_p_L))
      s_p[j] = s_p_R;
    else
      s_p[j] = s_p_L;

    if(s_rho_L * s_rho_R < 0.0)
      s_rho[j] = 0.0;
    else if(fabs(s_rho_R) < fabs(s_rho_L))
      s_rho[j] = s_rho_R;
    else
      s_rho[j] = s_rho_L;
  }

//===============================================

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
      { rho_L = RHO[k-1][j-1] + 0.5*(X[k-1][j]-X[k-1][j-1])*s_rho[j-1];
	u_L   =   U[k-1][j-1] + 0.5*(X[k-1][j]-X[k-1][j-1])*s_u[j-1];
	p_L   =   P[k-1][j-1] + 0.5*(X[k-1][j]-X[k-1][j-1])*s_p[j-1];
      }
      else
      { rho_L = RHOL[k-1];
	u_L = UL[k-1];
	p_L = PL[k-1];
      }

      if(j < m)
      { rho_R = RHO[k-1][j] - 0.5*(X[k-1][j+1]-X[k-1][j])*s_rho[j];
	u_R   =   U[k-1][j] - 0.5*(X[k-1][j+1]-X[k-1][j])*s_u[j];
	p_R   =   P[k-1][j] - 0.5*(X[k-1][j+1]-X[k-1][j])*s_p[j];
      }
      else
      { rho_R = RHOR[k-1];
	u_R = UR[k-1];
	p_R = PR[k-1];
      }


      if(j)
      { t_u_L   = s_u[j-1]/rho_L;
	t_p_L   = s_p[j-1]/rho_L;
	t_rho_L = s_rho[j-1]/rho_L;
      }
      else
      { t_rho_L = SRHOL[k-1]/rho_L;
	t_u_L = SUL[k-1]/rho_L;
	t_p_L = SPL[k-1]/rho_L;
      }

      if(j < m)
      { t_u_R   = s_u[j]/rho_R;
	t_p_R   = s_p[j]/rho_R;
	t_rho_R = s_rho[j]/rho_R;
      }
      else
      { t_rho_R = SRHOR[k-1]/rho_R;
	t_u_R = SUR[k-1]/rho_R;
	t_p_R = SPR[k-1]/rho_R;
      }//calculate the material derivative

      if((p_L < eps) || (p_R < eps) || (rho_L < eps) || (rho_R < eps)||isnan(p_L)||isnan(p_R)||isnan(u_L)||isnan(u_R)||isnan(rho_L)||isnan(rho_R))
	{
		printf("Error firstly happens on step=%d, grid=%d.\n", k, j);
		stop_step=1;
		break;
	}
//===============GRP scheme======================
	if(strcmp(scheme,"GRP")==0)
		{
			linear_GRP_solver_LAG(dire, mid, rho_L, rho_R, t_rho_L, t_rho_R, u_L, u_R, t_u_L, t_u_R, p_L, p_R, t_p_L, t_p_R, gamma, eps);
		}
	else
		{
			printf("No Riemann solver!\n");
			exit(7);
		}
	

      U_F[j] = mid[1] + 0.5*tau*dire[1];
      P_F[j] = mid[2] + 0.5*tau*dire[2];         

      X[k][j] = X[k-1][j] + tau*U_F[j];

      RHO_next_L[j] = mid[0] + tau*dire[0];
      RHO_next_R[j] = mid[3] + tau*dire[3];
      U_next[j] = mid[1] + tau*dire[1];
      P_next[j] = mid[2] + tau*dire[2];

//===============================================
 
   if(j) X_mass[j-1] = 0.5*(X[k][j-1]+X[k][j]);
   }
	
//===============THE CORE ITERATION=================(On Lagrange Coordinate)

    for(j = 0; j < m; ++j)
    { /*
       *  j-1          j          j+1
       * j-1/2  j-1  j+1/2   j   j+3/2  j+1
       *   o-----X-----o-----X-----o-----X--...
       */
	RHO[k][j] = 1 / ( 1/RHO[k-1][j] + tau/MASS[j]*(U_F[j+1] - U_F[j]) );
	U[k][j] = U[k-1][j] - tau/MASS[j]*(P_F[j+1] - P_F[j]);
	E[k][j] = E[k-1][j] - tau/MASS[j]*(P_F[j+1]*U_F[j+1] - P_F[j]*U_F[j]);
	P[k][j] = (E[k][j] - 0.5*U[k][j]*U[k][j]) * (gamma - 1.0) * RHO[k][j];
	/* forward Euler */

//-----------------------determind the slope-----------------------------

        s_u[j] = (  U_next[j+1] -   U_next[j])/(X[k][j+1]-X[k][j]);
        s_p[j] = (  P_next[j+1] -   P_next[j])/(X[k][j+1]-X[k][j]);
        s_rho[j] = (RHO_next_L[j+1] - RHO_next_R[j])/(X[k][j+1]-X[k][j]);
   }

    for(j = 0; j < m; ++j)
    {
      if(j)
      {
	s_u_L = s_time*(U[k][j] - U[k][j-1]) / (X_mass[j]-X_mass[j-1]);
	s_p_L = s_time*(P[k][j] - P[k][j-1]) / (X_mass[j]-X_mass[j-1]);
	s_rho_L = s_time*(RHO[k][j] - RHO[k][j-1]) / (X_mass[j]-X_mass[j-1]);
      }
      else
      {
	s_u_L = s_time*(U[k][j] - UL[k]) / (X_mass[j]-X[k][j])/2;
	s_p_L = s_time*(P[k][j] - PL[k]) / (X_mass[j]-X[k][j])/2;
	s_rho_L = s_time*(RHO[k][j] - RHOL[k]) / (X_mass[j]-X[k][j])/2;
      }
      if(j < m-1)
      {
	s_u_R = s_time*(U[k][j+1] - U[k][j]) / (X_mass[j+1]-X_mass[j]);
	s_p_R = s_time*(P[k][j+1] - P[k][j]) / (X_mass[j+1]-X_mass[j]);
	s_rho_R = s_time*(RHO[k][j+1] - RHO[k][j]) / (X_mass[j+1]-X_mass[j]);
      }
      else
      {
	s_u_R = s_time*(UR[k] - U[k][j]) / (X[k][j+1]-X_mass[j])/2;
	s_p_R = s_time*(PR[k] - P[k][j]) / (X[k][j+1]-X_mass[j])/2;
	s_rho_R = s_time*(RHOR[k] - RHO[k][j]) / (X[k][j+1]-X_mass[j])/2;
      }

      if(s_u_L * s_u_R < 0.0)
	s_u[j] = 0.0;
      else if((s_u_L * s_u[j]) < 0.0)
        s_u[j] = 0.0;
      else
      {
	slope_temp = ( (fabs(s_u_L) < fabs(s_u_R)) ? s_u_L : s_u_R );
	if( fabs(slope_temp) < fabs(s_u[j]) )
	    s_u[j] = slope_temp;
      }

      if(s_p_L * s_p_R < 0.0)
	s_p[j] = 0.0;
      else if((s_p_L * s_p[j]) < 0.0)
        s_p[j] = 0.0;
      else
      {
	slope_temp = ((fabs(s_p_L) < fabs(s_p_R)) ? s_p_L : s_p_R);
	if( fabs(slope_temp) < fabs(s_p[j]) )
	    s_p[j] = slope_temp;
      }

      if(s_rho_L * s_rho_R < 0.0)
	s_rho[j] = 0.0;
      else if((s_rho_L * s_rho[j]) < 0.0)
        s_rho[j] = 0.0;
      else
      {
	slope_temp = ((fabs(s_rho_L) < fabs(s_rho_R)) ? s_rho_L : s_rho_R);
	if( fabs(slope_temp) < fabs(s_rho[j]) )
	    s_rho[j] = slope_temp;
      }
    }

//---------------------------------------------------------------------------

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
