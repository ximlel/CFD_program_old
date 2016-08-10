int Sod_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int Free_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int odd_even_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int odd_even_EW_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int odd_even_EW_upstream_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int Shear_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int Cylinder_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int RMI_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
