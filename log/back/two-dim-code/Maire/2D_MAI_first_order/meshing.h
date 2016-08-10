int Sedov_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Sod_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);
int Sod_2material_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Saltzman_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Riemann_2D3_Tria_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);
int Riemann_2D3_Quad_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);

int Blunt_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n);
