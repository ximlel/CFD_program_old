#ifndef VARSTRUC_H
#define VARSTRUC_H

struct flu_var {
	double *RHO, *U, *V, *W, *P, *PHI, *gamma;
};

struct mesh_var {
	int **cell_point, *border[2];
	double *X, *Y, *Z;
};

extern double config[];

#define N_CONF 400

#define CONF_ERR(n)													\
	do {																\
		fprintf(stderr, "Error in the %d-th value of the configuration!", n); \
		exit(2);														\
	} while (0)

#define CONF_INI(i,j) printf("The %d-th value of configuration is initialized to %g .\n", i, j)

#endif
