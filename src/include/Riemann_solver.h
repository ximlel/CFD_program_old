#ifndef Riemann_solver_H
#define Riemann_solver_H

void ROE_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max, double delta);

void HLL_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max);

double Riemann_solver_exact(double * U_star, double * P_star, double gamma, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, int * CRW, double tol, int N);

void linear_GRP_solver_Edir
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);

void Roe_Goundov_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta);

void Roe_HLL_solver(double *V_mk, double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L, double P_R, double RHO_R, double U_R, double V_R, double *lambda_max, double delta);


#endif
