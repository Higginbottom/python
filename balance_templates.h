/* balance_abso.c */
int absolute_saha(PlasmaPtr www, PhotPtr p, double nh, double t_r, double t_e, double weight, int freq_sampling);
double absolute_total_emission(WindPtr one, double t_e, double f1, double f2);
double absolute_total_line_emission(PlasmaPtr ww, double t_e, double f1, double f2);
double absolute_line_heating(PlasmaPtr w, PhotPtr p, double ds);
double absolute_line_nsigma(struct lines *line_ptr, PlasmaPtr w);
double absolute_two_level_atom(struct lines *line_ptr, PlasmaPtr www, double *d1, double *d2);
double fb_cooling(double t_e, double f1, double f2, double *fb_h, double *fb_he1, double *fb_he2);
double fb_cooling_d(double f);
/* balance_bb.c */
int xbb(PhotPtr p, double t, double weight, double f1, double f2, int freq_sampling);
int xmake_bb(PhotPtr p, double t_r, double freqmin, double freqmax, double weight, int iphot_start, int nphot);
/* balance.c */
int main(int argc, char *argv[]);
int multicycle(PlasmaPtr www, PhotPtr p, int mode, int freq_sampling);
int print_levels(PlasmaPtr w);
int multi_concen(PlasmaPtr www);
int check_field(PlasmaPtr one);
int do_detail(PlasmaPtr one);
/* balance_gen.c */
int summary(PlasmaPtr one);
double line_heating(PlasmaPtr w, PhotPtr p, double ds);
double sobolev_line_heating(PlasmaPtr w, PhotPtr p, double ds);
/* balance_sub.c */
int xsaha(WindPtr one, PhotPtr p, double nh, double t_r, double t_e, double weight, int ioniz_mode, int freq_sampling);
int cycle(PlasmaPtr www, PhotPtr p, double nh, double t_r, double t_e, double weight, int mode, int freq_sampling);
int dumb_step(WindPtr one, double *te);
int find_te(WindPtr one);
/* bal_photon2d.c */
int translate(WindPtr w, PhotPtr pp, double tau_scat, double *tau, int *nres);
int translate_in_space(PhotPtr pp);
double ds_to_wind(PhotPtr pp);
int translate_in_wind(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres);
int walls(PhotPtr p, PhotPtr pold);
/* bal_trans_phot.c */
int trans_phot(WindPtr w, PhotPtr p, int iextract);
/* plane.c */
int pl_wind_define(WindPtr w);
double pl_rho(double x[]);
double pl_velocity(double x[], double v[]);
int define_wind_grid(WindPtr w, double xmin, double xmax, double zmin, double zmax);
int pl_copy_conditions(PlasmaPtr win, PlasmaPtr wout);
