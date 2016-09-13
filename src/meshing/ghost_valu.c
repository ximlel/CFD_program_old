#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"

void period_ghost(struct mesh_var * mv)
{
	mv->num_ghost = 0;
	int n_border;
	for(int l = 1, i = -1; l <= mv->num_border[0]; l++)
		{
			i++;
			n_border = i + mv->num_border[l]; 
			for( ; i < n_border; i++)
				{
					if (mv->border_cond[i] >= 0)
						mv->num_ghost++;
					p2_p = mv->border_pt[i+1];
					p2_n = mv->border_pt[i];
					if((p_p == p2_p && p_n == p2_n) || (p_p == p2_n && p_n == p2_p))
						{
							cv->cell_cell[k][j] = mv->border_cond[i];
							cell_rec = 1;
							break;
						}
				}							
			if (cell_rec)
				break;
		}
}
