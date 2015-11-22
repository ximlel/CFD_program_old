#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define max(x,y)  ( x>y?x:y )



void ROE_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta)
{	
	double H_L, H_R;
	H_L = gamma/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L);
	H_R = gamma/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R);

	F[0] = 0.5*(RHO_L*U_L+RHO_R*U_R);
	F[1] = 0.5*(RHO_L*U_L*U_L+P_L+RHO_R*U_R*U_R+P_R);
	F[2] = 0.5*(RHO_L*U_L*H_L+RHO_R*U_R*H_R);

	double RHO_S, U_S, H_S, C_S;
	RHO_S = sqrt(RHO_L*RHO_R);
	U_S = (U_L*sqrt(RHO_L)+U_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	H_S = (H_L*sqrt(RHO_L)+H_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	C_S = sqrt((gamma-1.0)*(H_S-0.5*U_S*U_S));

	
	double R[3][3];
	double lambda[3], W[3];
	R[0][0] = 1.0;
	R[0][1] = 1.0;
	R[0][2] = 1.0;
	R[1][0] = U_S - C_S;
	R[1][1] = U_S;
	R[1][2] = U_S + C_S;
	R[2][0] = H_S - U_S*C_S;
	R[2][1] = 0.5*(U_S*U_S);
	R[2][2] = H_S + U_S*C_S;

	int i, j;
	
	W[0] = 0.5*((P_R-P_L)-RHO_S*C_S*(U_R-U_L))/(C_S*C_S);
	W[1] = (RHO_R-RHO_L)-(P_R-P_L)/(C_S*C_S);
	W[2] = 0.5*((P_R-P_L)+RHO_S*C_S*(U_R-U_L))/(C_S*C_S);
	
		
	lambda[0] = fabs(U_S - C_S);
	lambda[1] = fabs(U_S);
 	lambda[2] = fabs(U_S + C_S);
	
	if(lambda[0]<delta)
			lambda[0] = 0.5/delta*(lambda[0]*lambda[0] + delta*delta);	
//	if(lambda[1]<delta)
//			lambda[1] = 0.5/delta*(lambda[1]*lambda[1] + delta*delta);		   
	if(lambda[2]<delta)
			lambda[2] = 0.5/delta*(lambda[2]*lambda[2] + delta*delta);		   

	*lambda_max = 0;
	for(i = 0; i < 3; i++)
		{
			*lambda_max = max(*lambda_max, lambda[i]);
			for(j = 0; j < 3 ; j++)
				{
					F[i] += -0.5*lambda[j]*W[j]*R[i][j];				
				}
		}
//	* lambda_max = 0.5*(fabs(U_S)+C_S);	  
}


