#ifndef U_MIN_BURGERS
#define U_MIN_BURGERS 0.0
#endif /* U_MIN_BURGERS */


int first_order_solver
(double * config, int m, 
 double * RHO[], double * U[], double * P[], double * E[], double * X[], double * cpu_time, char * scheme, double CFL);
