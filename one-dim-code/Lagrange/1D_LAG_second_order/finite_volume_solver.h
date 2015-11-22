#ifndef U_MIN_BURGERS
#define U_MIN_BURGERS 0.0
#endif /* U_MIN_BURGERS */


int second_order_solver
(double * config, int m, 
 double * RHO[], double * U[], double * P[], double * E[], double * X[], double * cpu_time, char * scheme, double CFL);

void linear_GRP_solver_LAG
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);
