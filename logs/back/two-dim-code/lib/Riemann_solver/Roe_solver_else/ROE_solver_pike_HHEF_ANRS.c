#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define max(x,y)  ( x>y?x:y )
#define min(x,y)  (x<y?x:y)



void ROE_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max, double delta)
{
	double const  Q_user = 2.0;
		
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
	double qn_R, qt_R;
	qn_R = U_R*n_x + V_R*n_y;
	qt_R = -U_R*n_y + V_R*n_x;
	double qn_L, qt_L;
	qn_L = U_L*n_x + V_L*n_y;
	qt_L = -U_L*n_y + V_L*n_x;


	double C_L, C_R;
	C_L = sqrt(gamma*P_L/RHO_L);
	C_R = sqrt(gamma*P_R/RHO_R);
	double z = 0.5 *  (gamma-1.0) / gamma;

	double Q, P_pvrs, P_max, P_min, RHO_bar, C_bar;
	P_min = min(P_L,P_R);
	P_max = max(P_L,P_R);
	Q = P_max/P_min;
	RHO_bar = 0.5*(RHO_L+RHO_R);
	C_bar = 0.5*(C_L+C_R);
	P_pvrs = 0.5*(P_L+P_R)+0.5*(qn_L-qn_R)*RHO_bar*C_bar;

	double A_L,A_R,B_L,B_R;
	A_L = 2.0/(gamma+1.0)/RHO_L;
	A_R = 2.0/(gamma+1.0)/RHO_R;
	B_L = (gamma-1)/(gamma+1)*P_L;
	B_R = (gamma-1)/(gamma+1)*P_R;

	double P_star, U_star, U_star_L, U_star_R, RHO_star_L, RHO_star_R, C_star_L, C_star_R, P_0, g_L_0, g_R_0, lambda_L_1, lambda_R_1, lambda_L_4, lambda_R_4;

	if(Q<Q_user&&P_min<P_pvrs&&P_pvrs<P_max) //PVRS
		{
			P_star = max(0,P_pvrs);
			U_star = 0.5*(qn_L+qn_R)+0.5*(P_L-P_R)/(RHO_bar*C_bar);
			RHO_star_L = RHO_L + (qn_L-U_star)*RHO_bar/C_bar;
			RHO_star_R = RHO_R + (U_star - qn_R)*RHO_bar/C_bar;
			C_star_L = sqrt(gamma*P_star/RHO_star_L);
			C_star_R = sqrt(gamma*P_star/RHO_star_R);
			U_star_L = U_star;
			U_star_R = U_star;
		}
	else if(P_pvrs<P_min) //TRRS
		{	   	
			P_star = pow((C_L + C_R - (gamma-1.0)/2.0*(qn_R - qn_L))/(C_L/pow(P_L,z) + C_R/pow(P_R,z)), 1.0/z);	
			C_star_L = C_L * pow(P_star/P_L, z);
			U_star_L = qn_L + 2.0/(gamma - 1.0)*(C_L - C_star_L);
			C_star_R = C_R * pow(P_star/P_R, z);
			U_star_R = qn_R + 2.0/(gamma - 1.0)*(C_star_R - C_R);
		}
	else //TSRS
		{
			P_0 = max(0,P_pvrs);
			g_L_0 = sqrt(A_L/(P_0+B_L));
			g_R_0 = sqrt(A_R/(P_0+B_R));
			P_star = (g_L_0*P_L+g_R_0*P_R-(qn_R-qn_L))/(g_L_0+g_R_0);
			U_star = 0.5*(qn_R+qn_L)+0.5*((P_star-P_R)*g_R_0-(P_star-P_L)*g_L_0);
			RHO_star_L = RHO_L*(P_star/P_L+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)*P_star/(gamma+1.0)/P_L+1.0);
			RHO_star_R = RHO_R*(P_star/P_R+(gamma-1.0)/(gamma+1.0))/((gamma-1.0)*P_star/(gamma+1.0)/P_R+1.0);
			C_star_L = sqrt(gamma*P_star/RHO_star_L);
			C_star_R = sqrt(gamma*P_star/RHO_star_R);
			U_star_L = U_star;
			U_star_R = U_star;
		}
	
	
	lambda_L_1 = qn_L - C_L;
	lambda_R_1 = U_star_L - C_star_L;
	lambda_L_4 = U_star_R + C_star_R;
	lambda_R_4 = qn_R + C_R;

	
	double R[4][4];
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

	int i, j;
	
	W[0] = 0.5*((P_R-P_L)-RHO_S*C_S*(qn_R-qn_L))/(C_S*C_S);
	W[1] = (RHO_R-RHO_L)-(P_R-P_L)/(C_S*C_S);
	W[2] = 0.5*((P_R-P_L)+RHO_S*C_S*(qn_R-qn_L))/(C_S*C_S);
	W[3] = RHO_S*(qt_R-qt_L);
	
		
	lambda[0] = fabs(qn_S - C_S);
	lambda[1] = fabs(qn_S);
 	lambda[2] = fabs(qn_S + C_S);
	lambda[3] = fabs(qn_S);
	

	*lambda_max = fabs(U_S) + C_S;
	
	if(lambda_L_1<0&&lambda_R_1>0)
		{
			F[0] = RHO_L*U_L*n_x+RHO_L*V_L*n_y;
			F[1] = (RHO_L*U_L*U_L+P_L)*n_x+RHO_L*U_L*V_L*n_y;
			F[2] = RHO_L*U_L*V_L*n_x+(RHO_L*V_L*V_L+P_L)*n_y;
 			F[3] = RHO_L*U_L*H_L*n_x+RHO_L*V_L*H_L*n_y;
			lambda[0] = lambda_L_1*(lambda_R_1-(qn_S-C_S))/(lambda_R_1-lambda_L_1);
			for(i = 0; i < 4; i++)
				{					
					F[i] += lambda[0]*W[0]*R[i][0];				
				}
		}
	else if(lambda_L_4<0&&lambda_R_4>0)
		{
			F[0] = RHO_R*U_R*n_x+RHO_R*V_R*n_y;
			F[1] = (RHO_R*U_R*U_R+P_R)*n_x+RHO_R*U_R*V_R*n_y;
			F[2] = RHO_R*U_R*V_R*n_x+(RHO_R*V_R*V_R+P_R)*n_y;
 			F[3] = RHO_R*U_R*H_R*n_x+RHO_R*V_R*H_R*n_y;
			lambda[2] = lambda_R_4*((qn_S+C_S)-lambda_L_4)/(lambda_R_4-lambda_L_4);
			for(i = 0; i < 4; i++)
				{
					F[i] += -lambda[2]*W[2]*R[i][2];				
				}
		}
	else
		{	
			double delta_1=0.003;
			double delta_2=0.2;
		
//			if(lambda[1]<delta_1)
//				lambda[1] = 0.5/delta_1*(lambda[1]*lambda[1] + delta_1*delta_1);		   
//			if(lambda[3]<delta_2)
//				lambda[3] = 0.5/delta_2*(lambda[3]*lambda[3] + delta_2*delta_2);

			for(i = 0; i < 4; i++)
				{
					for(j = 0; j<4; j++)
						{
							F[i] += -0.5*lambda[j]*W[j]*R[i][j];				
						}
				}
		}

}


