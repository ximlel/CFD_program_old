int second_order_solver
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[], int * BOUNDARY_POINT[],
 double * RHO[], double * U[], double * V[], double * P[], double * X[], double * Y[],
 double * NORMAL_VELOCITY[], double * gamma, double * cpu_time, char * limiter, double CFL);

int second_order_solver_GLACE
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[], int * BOUNDARY_POINT[],
 double * RHO[], double * U[], double * V[], double * P[], double * X[], double * Y[],
 double * NORMAL_VELOCITY[], double * gamma, double * cpu_time, char * limiter, double CFL);

int second_order_iteration_solver
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[], int * BOUNDARY_POINT[],
 double * RHO[], double * U[], double * V[], double * P[], double * X[], double * Y[],
 double * NORMAL_VELOCITY[], double * gamma, double * cpu_time, char * limiter, double CFL);
