
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "balance_templates.h"

//OLD #define LINELENGTH 132
#define NP    10000

/* xsaha determines the abundances using the electron temperature and the Saha equation
   and then calculates the heating and cooling 

  02jun	ksl	Moved initializations of heading etc. after ionization calculation in
 		an attempt to get xsaha to behave more like python since some of ionization
	        calculations use the current value of the the total heating.	
	06nov	ksl	58b - Began work to get this to work again
 */

int
xsaha (one, p, nh, t_r, t_e, weight, ioniz_mode, freq_sampling)
     WindPtr one;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int ioniz_mode;
     int freq_sampling;
{
  double total_fb ();
  double total_free ();
  double total_line_emission ();
  double total_emission ();
  double s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
  double fb_cooling ();
  double vol;
//  double freqmin, freqmax;
  int n;
  PlasmaPtr xplasma;


  xplasma = &plasmamain[one->nplasma];

/* Now determine the abundances, concentrations, and partition functions based
on the conditions as currently defined.

ion_abundances is a python routine.  Currently (Sept 01) legal values are

0   on-the-spot approximation using existing t_e
1   LTE (e.g. concentrations)
2   Fixed concentrations
3   On the spot, with one_shot at updating t_e before calculating densities

It is far from clear that all of these will work at present with balance.  

01sept ksl

*/

  if (ioniz_mode < 0)
    {
      Log
	("xsaha: Calculating heating and cooling without changing abundances\n");
    }
  else
    {
      if (ioniz_mode == 0)
	{
	  Log
	    ("xsaha: Calculating abundances from On_the_spot or nebular approx\n");
	}
      else if (ioniz_mode == 1)
	{
	  Log ("xsaha: Calculating abundances from LTE\n");
	}
      else if (ioniz_mode == 2)
	{
	  Log ("xsaha: Fixed concentrations\n");
	}
      else if (ioniz_mode == 3)
	{
	  Log ("xsaha: Abundances using one shot (one_shot at t update\n");
	}
      else if (ioniz_mode == 4)
	{
	  Log ("xsaha: Test of detailed balance\n");
	}

      else
	{
	  Error ("xsaha: ioniz_mode unknown. Good luck!\n");
	}

      ion_abundances (xplasma, ioniz_mode);
    }


/* initialize some parts of the wind ptr (see wind_rad_init) */

  xplasma->j = xplasma->ave_freq = xplasma->lum = xplasma->heat_tot =
    xplasma->ntot = 0;
  xplasma->heat_ff = xplasma->heat_photo = xplasma->heat_lines = 0.0;
  xplasma->heat_z = 0;
  xplasma->lum_z = 0.0;
  xplasma->ntot = xplasma->nioniz = 0;
  for (n = 0; n < nions; n++)
    {
      xplasma->ioniz[n] = 0;
      xplasma->recomb[n] = 0;
      xplasma->heat_ion[n] = 0;
      xplasma->lum_ion[n] = 0;
    }

  if (xplasma->t_e != t_e)
    Error ("xsaha xplasma->t_e %f != t_e %f. Why?\n", xplasma->t_e, t_e);
  if (xplasma->t_r != t_r)
    Error ("xsaha xplasma->t_r %f != t_r %f. Why?\n", xplasma->t_r, t_r);
  if (xplasma->w != weight)
    Error ("xsaha xplasma->w %f != weight %f. Why?\n", xplasma->w, weight);

  vol = wmain[xplasma->nwind].vol;
  s = vol;
//  s = xplasma->vol;

/* Now calculate the heating  based on these abundances, etc.

There appear to be 3 types of heating, one in which the line is assumed
to have a simple flat topped profile, one in which the heating is calculated
using Python (i.e. there is no absorption until the photon scatters in a MC
sense), and one with a simple sobolev_line assumption (in which case one calculates
tau to define how much got absorbed.)

There appearted to be an error in the routine below which I fixed up, just so that
the choices actually matched what was said.  I fixed this.  01sept ksl

*/

/* Create photons which comprise a weighted BB distribution */

/* 0810 - ksl - with next call xband->nbands, is set to 1, and the 
 * frequency limites depned on t_r 
 */

  init_bands (t_r, 0.0, 1e50, 0, &xband);

//???? New  make_bb (p, t_r, weight);
  xbb (p, t_r, weight, xband.f1[0], xband.f2[0], 1);

  if (geo.rt_mode == 2)
    {
      for (n = 0; n < NPHOT; n++)
	{
	  sobolev_line_heating (xplasma, &p[n], s);
	  radiation (&p[n], s);
	}
    }
  else if (geo.rt_mode == 1)
    {
      trans_phot (one, p, 0);
      xplasma[0].heat_tot += xplasma[0].heat_lines;
    }
  else
    {
      for (n = 0; n < NPHOT; n++)
	{
	  line_heating (xplasma, &p[n], s);
	  radiation (&p[n], s);
	}
    }

  luminosity = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

  return (0);
}


/***********************************************************
               Space Telescope Science Institute

Synopsis: 
	The routine cycle is intended to mimic what happens in Python in
	a very simple way, updating the ion abundances, etc. in a single
	cell, calculating the heating of the cell based on the updated
	abundances and then calling one shot to try to adjust the temperature.

   
Arguments:		

Returns:
 
 
Description:	

		
Notes:
	If t_r or t_e are changed then the routine will initialize everything
	to LTE, but otherwise it will go straight to the heating, and updating
	of abundances based on the new value of the heating.

	NB: t_r, t_e, weight are only used for generating the photon field. 
	The ionization calculation is based solely on what is contained in
	xplasma.


History:
	01dec	ksl	Updated to reflect new ways of calling routines
		    	like concentrations.	
	0810	ksl	67 - Relook and slight cleanup
 
**************************************************************/

double cy_tr = 0.0;
double cy_nh = 0.0;

int
cycle (xplasma, p, nh, t_r, t_e, weight, mode, freq_sampling)
     PlasmaPtr xplasma;
     PhotPtr p;
     double nh;
     double t_r;
     double t_e;
     double weight;
     int mode;
     int freq_sampling;
{
  double total_fb ();
  double total_free ();
  double total_line_emission ();
  double total_emission ();
  double ne, s;
  double alpha_recomb ();
  double planck (), kappa_ff ();
  double luminosity;
  double vol;
  int n;
  WindPtr one;

  one = &wmain[xplasma->nwind];

/* initialize some parts of the wind ptr (see wind_rad_init) */

  xplasma->j = xplasma->ave_freq = xplasma->lum = xplasma->heat_tot =
    xplasma->ntot = 0;
  xplasma->heat_ff = xplasma->heat_photo = xplasma->heat_lines = 0.0;
  xplasma->heat_z = 0;
  xplasma->lum_z = 0.0;
  xplasma->ntot = xplasma->nioniz = 0;
  for (n = 0; n < nions; n++)
    {
      xplasma->ioniz[n] = 0;
      xplasma->recomb[n] = 0;
      xplasma->heat_ion[n] = 0;
      xplasma->lum_ion[n] = 0;
    }

  if (xplasma->t_e != t_e)
    Error ("cycle: xsaha xplasma->t_e %f != t_e %f. Why?\n", xplasma->t_e,
	   t_e);
  if (xplasma->t_r != t_r)
    Error ("cycle: xsaha xplasma->t_r %f != t_r %f. Why?\n", xplasma->t_r,
	   t_r);
  if (xplasma->w != weight)
    Error ("cycle: xsaha xplasma->w %f != weight %f. Why?\n", xplasma->w,
	   weight);

  vol = wmain[xplasma->nwind].vol;
  s = vol;


/* 
Now set up the initial concentrations.  This is to parallel what is
done in python where we initially begin with an LTE like assumption.
On successive calls to cycle this section will be skipped as a result
of the if statement
 */

  if (t_r != cy_tr || nh != cy_nh)
    {
      Log
	("Cycle: Fixing abundances to Saha values since new t_r %g or nh %g\n",
	 t_r, nh);
      xplasma->t_e = 0.9 * t_e;	// Lucy guess
      nebular_concentrations (xplasma, 2);
      Log ("Cycle: On the spot estimate  t_r of %g gives ne %g\n", t_r, ne);
      cy_tr = t_r;
      cy_nh = nh;
      xplasma->dt_e = 0.0;	// Allow restart of oneshot limits
    }


/* Now calculate the heating */

/* First make some photons which comprise a BB distribution */

/* Set the frequency limits and then generate the photons */

  init_bands (t_r, 0.0, 1e50, 0, &xband);
  xbb (p, t_r, weight, xband.f1[0], xband.f2[0], freq_sampling);


/* Shift values to old */
  xplasma->heat_tot_old = xplasma->heat_tot;
  xplasma->dt_e_old = xplasma->dt_e;

/* Next calculate the heating for this distribution */
  for (n = 0; n < NPHOT; n++)
    {
      line_heating (xplasma, &p[n], s);
      radiation (&p[n], s);
    }

/* Now calculate the luminosity for these conditions */

  luminosity = total_emission (one, xband.f1[0], xband.f2[0]);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

/* Shift values to old */
  xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
  xplasma->t_e_old = xplasma->t_e;
  xplasma->t_r_old = xplasma->t_r;
  xplasma->lum_rad_old = xplasma->lum_rad;

  one_shot (xplasma, mode);


  Log ("Cycle: one_shot old t_r t_e %8.2g %8.2g, new t_r t_e %8.2g %8.2g\n",
       t_r, t_e, xplasma->t_r, xplasma->t_e);

//  Error ?? -- Next statement is probably superfluous, since calc_te, called by oneshot calculated total emission
//  luminosity = total_emission (xplasma, freqmin, freqmax);
  luminosity = total_emission (one, xband.f1[0], xband.f2[0]);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

/* Convergence check */
  convergence (xplasma);

  return (0);

}



/* 

The routine dumb_step increments te upward or downward depending on whether the total luminosity at
this temperature is less than or greater than the current value of the heating. It does not
recalculate the heating at any point.   The resulting te is incorporated into the windptr
	01dec	ksl	Updated to reflect new calls to routines like total_emission.  The routine
			seemed to work following these changes.
*/

int
dumb_step (one, te)
     WindPtr one;
     double *te;
{
  double lum_current, lum_hotter, lum_cooler;
  double tstep, hardness;
  double total_emission ();
  PlasmaPtr xplasma;

  xplasma = &plasmamain[one->nplasma];

  if (*te > 20000.)
    tstep = 0.1 * (*te);
  else
    tstep = 0.05 * (*te);

  xplasma->t_e = *te;
  lum_current = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  Log
    ("Current values of luminosities compared to previously calculated heating\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);

  xplasma->t_e = *te + tstep;
  lum_hotter = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  Log ("Results of raising t_e\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);

  xplasma->t_e = *te - tstep;
  lum_cooler = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  Log ("Results of lowering t_e\n");
  Log
    ("t_e %8.2e lum_tot  %8.2e lum_lines  %8.2e lum_ff  %8.2e lum_fb     %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_e, xplasma->lum_rad, xplasma->lum_lines, xplasma->lum_ff,
     xplasma->lum_fb, xplasma->lum_ion[0], xplasma->lum_ion[2],
     xplasma->lum_ion[3], xplasma->lum_z);
  Log
    ("t_r %8.2e heat_tot %8.2e heat_lines %8.2e heat_ff %8.2e heat_photo %8.2e %8.2e %8.2e %8.2e %8.2e\n",
     xplasma->t_r, xplasma->heat_tot, xplasma->heat_lines, xplasma->heat_ff,
     xplasma->heat_photo, xplasma->heat_ion[0], xplasma->heat_ion[2],
     xplasma->heat_ion[3], xplasma->heat_z);

  Log ("te %8.2g lum  %8.2e\n", *te, lum_current);
  Log ("te %8.2g lum  %8.2e\n", *te + tstep, lum_hotter);
  Log ("te %8.2g lum  %8.2e\n", *te - tstep, lum_cooler);

  Log ("tr %8.2g heat %8.2e\n", xplasma->t_r, xplasma->heat_tot);

  hardness =
    (lum_current - xplasma->heat_tot) / (lum_current + xplasma->heat_tot);

  if (fabs (hardness) < 0.1)
    {
      Log ("Convergence is achieved!\n");
    }
  else if (hardness > 0.0)
    {				// luminosity is greater the heating 

      if (lum_cooler < lum_hotter)
	{
	  *te -= tstep;
	  Log ("Dropping t_e to %8.2g\n", *te);
	}
      else
	{
	  *te += tstep;
	  Log ("Raising t_e to % 8.2g\n", *te);
	}
    }
  else
    {				//luminosity is less than heating

      if (lum_cooler > lum_hotter)
	{
	  *te -= tstep;
	  Log ("Dropping t_e to %8.2g\n", *te);
	}
      else
	{
	  *te += tstep;
	  Log ("Raising t_e to % 8.2g\n", *te);
	}
    }

  xplasma->t_e = *te;
  return (0);
}

int
find_te (one)
     WindPtr one;
{
  double calc_te ();
  double luminosity, total_emission ();
  PlasmaPtr xplasma;

  xplasma = &plasmamain[one->nplasma];

  luminosity = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

  xplasma->t_e = calc_te (xplasma, TMIN, 1.2 * xplasma->t_r);

  luminosity = total_emission (one, 0.0, VERY_BIG);
  num_recomb (&xplasma[0], xplasma->t_e);

  summary (xplasma);

  return (0);
}
