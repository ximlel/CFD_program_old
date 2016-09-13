struct cell_var cell_mem_init(struct mesh_var * mv);

void cons_qty_init(struct cell_var * cv, struct flu_var * FV);

void vol_comp(struct cell_var * cv, struct mesh_var * mv);

void cell_pt_clockwise(struct mesh_var * mv);

void cell_rel(struct cell_var * cv, struct mesh_var * mv);

void cell_centroid(struct cell_var * cv, struct mesh_var * mv);



void slope_limiter(struct cell_var * cv, struct mesh_var * mv, struct flu_var * FV);



void Euler_scheme(struct flu_var * FV, struct mesh_var * mv, char * scheme);


/*
void first_order_solver
(int * STEP, double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[],
 int * BOUNDARY_POINT[], int m, int n, double * RHO[], double * U[], double * V[], double * P[],
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL, char * example);

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
 double * X, double * Y, double * gamma, double * cpu_time, char * scheme, double CFL, char * example);
*/
