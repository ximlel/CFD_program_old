#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define max(x,y)  ( x>y?x:y )

#include "../custom.h"


void ROE_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max, double delta)
{	
	double H_L, H_R;
	H_L = gamma/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L+V_L*V_L);
	H_R = gamma/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R+V_R*V_R);

	F[0] = 0.5*(RHO_L*U_L+RHO_R*U_R)*n_x + 0.5*(RHO_L*V_L+RHO_R*V_R)*n_y;
	F[1] = 0.5*(RHO_L*U_L*U_L+P_L+RHO_R*U_R*U_R+P_R)*n_x + 0.5*(RHO_L*U_L*V_L+RHO_R*U_R*V_R)*n_y;
	F[2] = 0.5*(RHO_L*U_L*V_L+RHO_R*U_R*V_R)*n_x + 0.5*(RHO_L*V_L*V_L+P_L+RHO_R*V_R*V_R+P_R)*n_y;
	F[3] = 0.5*(RHO_L*U_L*H_L+RHO_R*U_R*H_R)*n_x+0.5*(RHO_L*V_L*H_L+RHO_R*V_R*H_R)*n_y;

	double RHO_S, U_S, V_S, H_S, C_S;
	RHO_S = sqrt(RHO_L*RHO_R);
	U_S = (U_L*sqrt(RHO_L)+U_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	V_S = (V_L*sqrt(RHO_L)+V_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	H_S = (H_L*sqrt(RHO_L)+H_R*sqrt(RHO_R)) / (sqrt(RHO_L)+sqrt(RHO_R));
	C_S = sqrt((gamma-1.0)*(H_S-0.5*U_S*U_S-0.5*V_S*V_S));

	double qn_S, qt_S;
	qn_S = U_S*n_x + V_S*n_y;
	qt_S = -U_S*n_y + V_S*n_x;
	
	double R[4][4];
	double L[4][4];
	double lambda[4], W[4];
	R[0][0] = 1.0;
	R[0][1] = 1.0;
	R[0][2] = 1.0;
	R[0][3] = 0.0;
	R[1][0] = U_S - C_S*n_x;
	R[1][1] = U_S;
	R[1][2] = U_S + C_S*n_x;
	R[1][3] = -n_y;
	R[2][0] = V_S - C_S*n_y;
	R[2][1] = V_S;
	R[2][2] = V_S + C_S*n_y;
	R[2][3] = n_x;
	R[3][0] = H_S - qn_S*C_S;
	R[3][1] = 0.5*(U_S*U_S+V_S*V_S);
	R[3][2] = H_S + qn_S*C_S;
	R[3][3] = qt_S;

	int i, j, k;

	for(i = 0; i < 4; i++)
		for(j = 0; j<4; j++)
			{
				L[i][j]=R[i][j];
			}	

	
	k = rinv(L[0],4);
	if(!k)
		{
			printf("ERROR: matrix R is not reversible.");
			exit(8);
		}

	
	double E_L, E_R;
	E_L = 1.0/(gamma-1.0)*P_L/RHO_L + 0.5*(U_L*U_L+V_L*V_L);
	E_R = 1.0/(gamma-1.0)*P_R/RHO_R + 0.5*(U_R*U_R+V_R*V_R);
	
	double delta_U[4];
	delta_U[0] = RHO_R - RHO_L;
	delta_U[1] = RHO_R*U_R - RHO_L*U_L;
	delta_U[2] = RHO_R*V_R - RHO_L*V_L;
	delta_U[3] = RHO_R*E_R - RHO_L*E_L;
	
	for(i = 0; i < 4; i++)
		{
			W[i] = 0.0;
			for(j = 0; j<4; j++)
				{
					W[i] += L[i][j]*delta_U[j];				
				}
		}
	
	lambda[0] = fabs(qn_S - C_S);
	lambda[1] = fabs(qn_S);
 	lambda[2] = fabs(qn_S + C_S);
	lambda[3] = fabs(qn_S);

	
	double delta_1=0.1;
	double delta_2=0.18;
	
//	if(lambda[0]<delta)
//			lambda[0] = 0.5/delta*(lambda[0]*lambda[0] + delta*delta);	
//	if(lambda[1]<delta_1)
//			lambda[1] = 0.5/delta_1*(lambda[1]*lambda[1] + delta_1*delta_1);		   
//	if(lambda[2]<delta)
//			lambda[2] = 0.5/delta*(lambda[2]*lambda[2] + delta*delta);
//	if(lambda[3]<delta_2)
//			lambda[3] = 0.5/delta_2*(lambda[3]*lambda[3] + delta_2*delta_2);		   

	*lambda_max = 0;
	for(i = 0; i < 4; i++)
		{
			*lambda_max = max(*lambda_max, lambda[i]);
			for(j = 0; j<4; j++)
				{
					F[i] += -0.5*lambda[j]*W[j]*R[i][j];				
				}
		}
}


