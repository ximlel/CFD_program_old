#include <math.h>
#include <stdio.h>



void GLACE(double (* M_pc)[2], double gamma, double P, double RHO, double l_pc, double n_pc_x, double n_pc_y)
{
	double z;
	z = RHO * sqrt(gamma * P / RHO);
	M_pc[0][0] = z * l_pc * n_pc_x * n_pc_x;
	M_pc[0][1] = z * l_pc * n_pc_x * n_pc_y;
	M_pc[1][0] = z * l_pc * n_pc_y * n_pc_x;
	M_pc[1][1] = z * l_pc * n_pc_y * n_pc_y;
}
