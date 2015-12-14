void first_order_solver
(int * STEP, double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL);

void first_order_two_species_solver
(int * STEP, double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[], double * Z[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL);

void linear_GRP_solver_Edir
(double * direvative, double * mid,
 double rho_L, double rho_R, double s_rho_L, double s_rho_R,
 double   u_L, double   u_R, double   s_u_L, double   s_u_R,
 double   p_L, double   p_R, double   s_p_L, double   s_p_R,
 double gamma, double eps);

void linear_GRP_solver_Edir_2D
(double *wave_speed, double *direvative, double *source, double lambda,
 double rho_L, double rho_R, double d_rho_L, double d_rho_R,
 double   u_L, double   u_R, double   d_u_L, double   d_u_R,
 double   v_L, double   v_R, double   d_v_L, double   d_v_R,
 double   p_L, double   p_R, double   d_p_L, double   d_p_R,
 double gamma, double eps);

void slope_limiter_Ven
(double * X_c, double * Y_c, double  * X, double * Y,
 double * grad_W_x, double * grad_W_y, 
 double * W[],  int NUM_CELL, double * config,
 int * CELL_CELL[], int * CELL_POINT[],
 int m, int n);


void second_order_solver
(int * STEP, double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL);

void second_order_two_species_solver
(int * STEP, double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[], double * Z[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL);
