#ifndef TOOLS_H
#define TOOLS_H

/*!\file tools.h
 * \brief  Some mathematical algorithm functions.
 * \author Lei Xin
 */

void initialize_memory(double * p[],int N,int * CELL_POINT[]);

void initialize_memory_int(int * p[],int N,int * CELL_POINT[]);

int CreateDir(const char* pPath);

/*!\brief A function to caculate the inverse of the input square matrix.
 * \param a The pointer of the input square matrix.
 * \param[in] n The order of the input square matrix.
 */
int rinv(double a[], int n);

double rnd( double *r);

/*!\brief  \f$\mu\f$ of Barth Jesperse limiter.
 * \param[in] x Variable \f$x\f$ in \f$\mu(x)\f$
 * \return Value of \f$\mu(x)\f$
 */
double miu_BJ(double x);

/*!\brief  \f$\mu\f$ of Venkatakrishnan limiter.
 * \param[in] x Variable \f$x\f$ in \f$\mu(x)\f$
 * \return Value of \f$\mu(x)\f$
 */
double miu_Ven(double x);

void Gauss_elimination(int n, double (*a)[n+1], double *x);

#endif
