#ifndef VARSTRUC_H
#define VARSTRUC_H

double * RHO = NULL;
double * U = NULL;
double * P = NULL;
double * V = NULL;
double * PHI = NULL;

struct flu_var {
	double *RHO, *U, *V, *W, *P, *PHI, *gamma;
};

struct mesh_var {
	int **cell_point, *border[2];
	double *X, *Y, *Z;
};

extern double config[];

#endif
