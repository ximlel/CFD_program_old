void initialize_memory(double * p[],int N,int * CELL_POINT[]);

void initialize_memory_int(int * p[],int N,int * CELL_POINT[]);

int CreateDir(const char* pPath);

int rinv(double a[], int n);

double rnd( double *r);

double miu_BJ(double x);

double miu_Ven(double x);

void Gauss_elimination(int n, double (*a)[n+1], double *x);
