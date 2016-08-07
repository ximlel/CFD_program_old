int Sedov_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Sod_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Saltzman_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Shock_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Blunts_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);
