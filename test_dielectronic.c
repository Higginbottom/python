 /***********************************************************
                                       Southampton university

 Synopsis:
	test_dielectronic is a stand alone code to produce data relating to trying to understand dielectronic recombination in all its forms!
 
Arguments:		


Returns:
 
Description:	
	
Notes:

History:
	 
	0512 - nsh wrote it
 	


**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
struct topbase_phot *fb_xtop;
struct photoionization *fb_xver;
double g1, g2, xne, xip;
#define SAHA 4.82907e15
double xip, xxxne, qromb_temp, dnu;
double fbt;			// Temperature at which thee emissivity is calculated
int fbfr;
#include "python.h"
int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;
  PhotPtr p;
  int n, i, j, nion, ion_lower;
  int ntemp, nne, temp_start, tempdiv, ntmin, nvmin;
  double ne[11], temp1, nh, xtemp, weight, xsaha, t_r, t_e, temp, zeta,
    alpha_all, alpha_gs, comp_zeta;
  double b, *partition, recomb_fudge, pi_fudge, gs_fudge, test, integral;
  double fthresh, fmax;
  double temp_func ();
  int compute_dr_coeffs ();
  FILE *fp_c_zeta;
  FILE *fp_c_die;
  FILE *fp_c_rec_all;
  FILE *fp_c_rec_gs, *fp_he_rec_gs;
  FILE *fp_c_comp_zeta;
  FILE *fig_data;
  double qromb ();
  double fb_topbase_partial ();

  PlasmaPtr xplasma;
  partition = xplasma->partition;	/* Set the partition function array to that held for the cell */
  temp_start = 299;
  ntemp = 401;
  tempdiv = 100.0;
  nne = 1;
  nh = 1e10;
  weight = 0.01;
  xne = nh;

//  fp_c_zeta=fopen("carbon_zeta.out","w");
//  fp_c_die=fopen("carbon_die.out","w");
//  fp_c_rec_all=fopen("carbon_rec_all.out","w");
  fp_c_rec_gs = fopen ("carbon_rec_gs.out", "w");
//  fp_he_rec_gs=fopen("helium_rec_gs.out","w");
//  fp_c_comp_zeta=fopen("carbon_comp_zeta.out","w");
//  fig_data=fopen("badnell_fig.out","w");

  strcpy (geo.atomic_filename, "atomic/standard73");
  printf ("atomic_filename=%s\n", geo.atomic_filename);
  get_atomic_data (geo.atomic_filename);


//     init_freebound (1.e2, 1.e9, 0, 1e50);
  for (i = 200; i < 801; i++)
    {
      temp = pow (10, i / 100.);
//              compute_dr_coeffs (temp);
//              fprintf(fp_c_zeta," %6.3e",temp);
//              fprintf(fp_c_die," %6.3e",temp);
//              fprintf(fp_c_rec_all," %6.3e",temp);
      fprintf (fp_c_rec_gs, " %6.3e", temp);
//              fprintf(fp_he_rec_gs," %6.3e",temp);
//              fprintf(fp_c_comp_zeta," %6.3e",temp);
//              fprintf(fig_data," %6.3e",temp);
      fbt = temp;		/* Externally transmitted variable */
      fbfr = 2;
      printf ("Temperature =%f\n", temp);
//              j=41;
//              alpha_all = xinteg_fb (temp, 1e14, 1e19, j-1, 2);
//              fprintf(fig_data," %6.3e",alpha_all);
//              fprintf(fig_data," %6.3e",dr_coeffs[j]);
//              j=54;
//              alpha_all = xinteg_fb (temp, 1e14, 1e19, j-1, 2);
//              fprintf(fig_data," %6.3e",alpha_all);
//              fprintf(fig_data," %6.3e",dr_coeffs[j]);
//              fprintf(fig_data,"\n");
      for (j = 6; j < 12; j++)
	{
//                      printf ("Ion %i, is z=%i, istate=%i\n",j,ion[j].z,ion[j].istate);               
//                      zeta=compute_zeta (temp, j, 1e14, 1e19, 1);
//                      alpha_all = xinteg_fb (temp, 1e14, 1e19, j-1, 2);
	  if (ion[j - 1].phot_info == 1)	//topbase
	    {
	      ntmin = ion[j - 1].ntop_ground;
//                              printf ("We are looking at topbase data for ion %i \n",j-1);
	      fb_xtop = &phot_top[ntmin];
	      fthresh = fb_xtop->freq[0];
	      fmax = fb_xtop->freq[fb_xtop->np - 1];
	      dnu = 100.0 * (fbt / H_OVER_K);
	      if (fthresh + dnu < fmax)
		{
		  fmax = fthresh + dnu;
		}
	      alpha_gs = qromb (fb_topbase_partial, fthresh, fmax, 1e-5);
	    }
	  else if (ion[j - 1].phot_info == 0)	// verner
	    {
	      nvmin = j - 1;
	      n = nvmin;
	      fb_xver = &xphot[ion[n].nxphot];
//                              printf ("We are looking at verner data for ion %i \n",j);
	      fthresh = fb_xver->freq_t;
	      fmax = fb_xver->freq_max;
	      dnu = 100.0 * (fbt / H_OVER_K);
	      if (fthresh + dnu < fmax)
		{
		  fmax = fthresh + dnu;
		}
	      alpha_gs = qromb (fb_verner_partial, fthresh, fmax, 1e-6);
	    }
//                      printf ("Integrating cross section %i from %e to %e and getting %e\n",ntmin,fthresh,fmax,alpha_gs);     
//                      comp_zeta=alpha_gs/alpha_all;
//                      fprintf(fp_c_zeta," %6.3e",zeta);
//                      fprintf(fp_c_die," %6.3e",dr_coeffs[j]);
//                      fprintf(fp_c_rec_all," %6.3e",alpha_all);
	  fprintf (fp_c_rec_gs, " %6.3e", alpha_gs);
//                      fprintf(fp_c_comp_zeta," %6.3e",comp_zeta);
	}
//                      ntmin = ion[2].ntop_ground;
//                      fb_xtop = &phot_top[ntmin];
//              fthresh = fb_xtop->freq[0];
//              fmax = fb_xtop->freq[fb_xtop->np - 1];
//              fprintf(fp_he_rec_gs," %6.3e\n",qromb (fb_topbase_partial, fthresh,fmax, 1e-4));
//              fprintf(fp_c_rec_all,"\n");
//              fprintf(fp_c_zeta,"\n");
//              fprintf(fp_c_die,"\n");
      fprintf (fp_c_rec_gs, "\n");
//              fprintf(fp_c_comp_zeta,"\n");
    }





  return EXIT_SUCCESS;
}
