struct cell_var cell_mem_init(const struct mesh_var mv);
void cons_qty_init(struct cell_var * cv, const struct flu_var FV);
void vol_comp(struct cell_var * cv, const struct mesh_var mv);
void cell_pt_clockwise(struct mesh_var * mv);
void cell_rel(struct cell_var * cv, const struct mesh_var mv);
void cell_centroid(struct cell_var * cv, const struct mesh_var mv);


void slope_limiter(struct cell_var * cv, const struct mesh_var mv, const struct flu_var FV);


int interface_var_init
(const struct cell_var cv, const struct mesh_var mv, struct i_f_var * ifv,
 struct i_f_var * ifv_R, const int k, const int j);
double tau_calc(const struct cell_var cv, const struct mesh_var mv);
void cons_qty_update(struct cell_var * cv, const struct mesh_var mv, const double tau);
void fluid_var_update(struct flu_var *FV, struct cell_var cv);


void Roe_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void HLL_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void Riemann_exact_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);


void Euler_scheme(struct flu_var *FV, const struct mesh_var mv, const char *scheme);
