#include <math.h>
#include <stdio.h>



void GLACE_dt(double (* M_pc)[2], double gamma, double P, double RHO, double l_pc, double n_pc_x, double n_pc_y, double dt_l_n_pc_x, double dt_l_n_pc_y, double dt_n_pc_x, double dt_n_pc_y)
{
	double z;
	z = RHO * sqrt(gamma * P / RHO);
	M_pc[0][0] = z * (dt_l_n_pc_x * n_pc_x + l_pc * n_pc_x * dt_n_pc_x);
	M_pc[0][1] = z * (dt_l_n_pc_x * n_pc_y + l_pc * n_pc_x * dt_n_pc_y);
	M_pc[1][0] = z * (dt_l_n_pc_y * n_pc_x + l_pc * n_pc_y * dt_n_pc_x);
	M_pc[1][1] = z * (dt_l_n_pc_y * n_pc_y + l_pc * n_pc_y * dt_n_pc_y);
}
