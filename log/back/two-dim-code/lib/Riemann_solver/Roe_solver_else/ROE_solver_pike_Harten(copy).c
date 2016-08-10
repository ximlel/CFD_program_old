#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void ROE_solver(double *F, double gamma, double pL, double rhoL, double uL, double vL,double n_x, double n_y, double pR, double rhoR, double uR,double vR,double *lambda_max, double delta)
{  


	// Local variables
	double m_x, m_y;                  // Tangent vector: mx*nx+my*ny = 0
	double unL, unR, umL, umR;      // Normal and tangent velocities
	double aL, aR, HL, HR;          // Speeds of sound.
	double RT,rho,u,v,H,a,un, um;   // Roe-averages
	double drho,dun,dum,dp,LdU[5];  // Wave strenghs
	double ws[5], Rv[5][5];          // Wave speeds and right-eigevectors
	double fL[5], fR[5], diss[5];   // Fluxes ad dissipation term
	double dws[5] ;                 // User-specified width for entropy fix
	int i, j;
	
	m_x = -n_y;
	m_y =  n_x;

	//Primitive and other variables.
	//  Left state

	unL = uL*n_x+vL*n_y;
	umL = uL*m_x+vL*m_y;
	aL = sqrt(gamma*pL/rhoL);
	HL = aL*aL/(gamma-1.0) + 0.5*(uL*uL+vL*vL);
	  // Right state
	
	unR = uR*n_x+vR*n_y;
	umR = uR*m_x+vR*m_y;
	aR = sqrt(gamma*pR/rhoR);
	HR = aR*aR/(gamma-1) + 0.5*(uR*uR+vR*vR);

	//First compute the Roe Averages
    	RT = sqrt(rhoR/rhoL);
	rho = RT*rhoL;
	u = (uL+RT*uR)/(1.0+RT);
	v = (vL+RT*vR)/(1.0+RT);
	H = (HL+RT* HR)/(1.0+RT);
	a = sqrt( (gamma-1.0)*(H-0.5*(u*u+v*v)) );
    	un = u*n_x+v*n_y;
  	um = u*m_x+v*m_y;

	//Wave Strengths
	drho = rhoR - rhoL ;
	dp =   pR - pL;
  	dun =  unR - unL;
    	dum =  umR - umL;

	LdU[1] = (dp - rho*a*dun )/(2.0*a*a);
	LdU[2] = rho*dum;
	LdU[3] =  drho - dp/(a*a);
	LdU[4] = (dp + rho*a*dun )/(2.0*a*a);

	//Wave Speed
	ws[1] = fabs(un-a);
	ws[2] = fabs(un);
	ws[3] = fabs(un);
	ws[4] = fabs(un+a);

	//Harten's Entropy Fix JCP(1983), 49, pp357-393:
	//only for the nonlinear fields.
	dws[1] = delta;
	if ( ws[1] < dws[1] )
		ws[1] = 0.5 * (ws[1]*ws[1]/dws[1]+dws[1] );
	dws[4] = delta;
	if ( ws[4] < dws[4] )
		ws[4] = 0.5 * ( ws[4]*ws[4]/dws[4]+dws[4] );

	//Right Eigenvectors
	Rv[1][1] = 1.0   ;
	Rv[2][1] = u - a*n_x;
	Rv[3][1] = v - a*n_y;
	Rv[4][1] = H - un*a;

	  Rv[1][2] = 1.0;
	  Rv[2][2] = m_x;
	  Rv[3][2] = m_y;
	  Rv[4][2] = um;

	  Rv[1][3] = 1.0;
	  Rv[2][3] = u;
	  Rv[3][3] = v ;
	  Rv[4][3] = 0.5*(u*u+v*v);

	  Rv[1][4] = 1.0;
	  Rv[2][4] = u + a*n_x;
	  Rv[3][4] = v + a*n_y;
	  Rv[4][4] = H + un*a;

	//Dissipation Term
	  for(i=1;i<=4;i++)
		  {
			  diss[i] = 0.0;		  		  
			  for(j=1;j<=4;j++)
				  diss[i] = diss[i] + ws[j]*LdU[j]*Rv[i][j];
		  }

	//Compute the flux.
	  fL[1] = rhoL*unL;
	  fL[2] = rhoL*unL * uL + pL*n_x;
	  fL[3] = rhoL*unL * vL + pL*n_y;
	  fL[4] = rhoL*unL * HL;

	  fR[1] = rhoR*unR;
	  fR[2] = rhoR*unR * uR + pR*n_x;
	  fR[3] = rhoR*unR * vR + pR*n_y;
	  fR[4] = rhoR*unR * HR;

	  for(i = 0; i < 4; i++)
					F[i] = 0.5 * (fL[i+1] + fR[i+1] - diss[i+1]);				
					
	  * lambda_max = fabs(un) + a;  //Normal max wave speed times half


}


