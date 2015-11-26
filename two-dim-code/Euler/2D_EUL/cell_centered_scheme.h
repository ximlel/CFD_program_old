int first_order_solver
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], double * RHO[], double * U[], double * V[], double * P[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL);

int first_order_two_species_solver
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[], double * Z[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL);

void linear_GRP_solver_Edir
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);