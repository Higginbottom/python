


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111001	ksl	Nearly all of this was coded by Nick in 2011, but ksl 
		added header this data, and reformated with indent
		this date so that the code was easier to read and so
		he could debug problems associated with indent.
111126	ksl	Renamed this routine power_sub.c in partial attempt
		to spread the blame for this particular approach
		to calculating the ionization equilibrium
120127  nsh	Changed the storage arrays fudge_store, num_store
		denom_store and ion_store to have dimension NIONS 
		these arrays were overflowing. Only picked up in
		Jan 12, when in python71 they started overflowing 
		and zeroing geo.wcycle. 

 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

struct photoionization *sim_xver;	//Verner & Ferland description of a photoionization x-section
struct topbase_phot *sim_xtop;	//Topbase description of a photoionization x-section
double sim_te;			//This is a storage variable for the current electron temperature so it is available for qromb calls

// These are variables need by tb_pow
PlasmaPtr xxxplasma;
double weight;			//This is a storage variable for the current PL geometric weight so it can be used in qromb
double xsim_alpha;		//This is a storage variable so the value of alpha (PL spectral index) can be used in qromb
double xsim_w;			//Storage variable for the sim W factor (either computed from photons going thruogh a cell, or before photon generation it is the power law constant x radiative weight / 4pi */

// 120127 NSH changed the dimension of following arrays from 300 to NIONS to prevent overflow
double fudge_store[NIONS];	//this is a store so we can see how the sim fudge factor varies
double num_store[NIONS];		// store for the numbnerator of the sim factor (PL)
double denom_store[NIONS];	// store for the denominator of the sim factor (BB)
double ion_store[NIONS];

//#define SAHA 4.82907e15		/* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
//#define MAXITERATIONS	200
//#define FRACTIONAL_ERROR 0.03
//#define THETAMAX	 1e4
//#define MIN_TEMP         100.  Put into python.h 

#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111217	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/

int
sim_driver (xplasma)
     PlasmaPtr xplasma;
{
  int nelem, nion, niterate;
  double xne, xnew;
  double newden[NIONS];
  double t_r, nh;
  double t_e, www;
  double f1, f2;
  int idiag;

  xxxplasma = xplasma;

  t_r = xplasma->t_r;
  t_e = xplasma->t_e;
  www = weight = xplasma->w;

  printf ("DDDDDDD We are just going to calculate DR rates\n");
  compute_dr_coeffs (t_e);	//Populate the dielectronic recombination coefficient array for the electron temperature for this cell. If none have been read in, then they will all be zero, and nothing will happen different.
  printf ("DDDDDDD and we have returned\n");



  /* 71 - 111229 - ksl - Small modifications to reflect moving nxfreq, xfreq into the geo structure */
  f1 = geo.xfreq[0];		/*NSH 1108 Use the lower bound of the power law estimator bands */
  f2 = geo.xfreq[geo.nxfreq];	/*NSH 1108 and the upper bound of the power law estimator bands */

  Log ("sim_driver: t_e=%f, t_r=%f, f1=%e, f2=%e\n", t_e, t_r, f1, f2);


  /* Initally assume electron density from the LTE densities */

  xne = xplasma->ne;
  if (xne < DENSITY_MIN)
    {
      Error ("sim_driver: Very low ionization: ne initially %8.2e\n", xne);
      xne = DENSITY_MIN;
    }

  nh = xplasma->rho * rho2nh;	//LTE -- Not clear needed at this level


//  printf ("Starting off with ne=%e and nh=%e\n", xne, nh);

  /* ksl - This is the standard loop over elements. Except for the call to sim_pl it is identical I expect
   * to the corresponding routing for saha.  And this suggests that one should consolodate this */

  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
      for (nelem = 0; nelem < nelements; nelem++)
	{
	  /* Do a one cycle for a single elementt */
	  idiag = sim_pl (nh, t_r, t_e, www, nelem, xplasma->ne,
			  xplasma->density, xne, newden, f1, f2);
	  if (idiag)
	    {
	      Error
		("sim_driver: sim_pl reported problem for element %d in wind cell %d of type %d\n",
		 nelem, xplasma->nwind, wmain[xplasma->nwind].inwind);
	      exit (0);
	    }
	  /* Re solve for the macro atom populations with the current guess for ne */
	  if (geo.macro_ioniz_mode == 1)
	    {
	      macro_pops (xplasma, xne);
	    }
	}

      /* OK now look at the macro atom information and if this is a macro atom calculate the abundances
       * from that
       */
      for (nion = 0; nion < nions; nion++)
	{
	  /* if the ion is being treated by macro_pops then use the populations just computed */
	  if ((ion[nion].macro_info == 1) && (geo.macro_simple == 0)
	      && (geo.macro_ioniz_mode == 1))
	    {
	      newden[nion] = xplasma->density[nion];
	    }

	  /*Set some floor so future divisions are sensible */
	  if (newden[nion] < DENSITY_MIN)
	    newden[nion] = DENSITY_MIN;
	}

      /* OK get the new ne and see whether we are close enough */
      xnew = get_ne (newden);
//      printf ("current estimate of ne after loop %i=%e\n",niterate,xnew);
      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;

      /* Check to see whether the search for xne has converged and if so exit loop */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      /* else start another iteration of the main loop */
      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }
  /* End of main iteration loop */

  if (niterate == MAXITERATIONS)
    {
      Error
	("sim_driver: failed to converge:nh %8.2e www %8.2e t_e %8.2e  t_r %8.2e \n",
	 nh, www, t_e, t_r);
      return (-1);
    }

/* Finally transfer the calculated densities to the real density array */

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nion] = newden[nion];
//      printf("Density of ion %i following sim correction is %e\n",nion,newden[nion]);

	}
    }
  return (0);
}




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
  This is where we solve for the ionization of a single element given many variables
  inlduding the electon density.  

  Arguments:  		


  Returns:

  Results are obtained in the array newden

  Notes:

  This routine parallels lucy_mazzali1 and lucy in saha.c.  They really ought to be combiend into
  a single routine as the only thing that is different is the way the fudge factors are calculated.


 

  History:
111217	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/

int
sim_pl (nh, t_r, t_e, www, nelem, ne, density, xne, newden, f1, f2)
     double nh, t_r, t_e, www;
     int nelem;
     double ne, density[], xne, newden[], f1, f2;
{
  double fudge, dummy, interpfrac;
  double fudge2, q;
  double sum, a;
  int ilow, ihi;
  int first, last, nion;
  double numerator, denominator;
  double max_ratio, current_ratio;
  FILE *fopen (), *fudgefile, *numfile, *denomfile, *ionfile;



  weight = www;




  if (t_e > MIN_TEMP)
    {
      fudge = 1.;

      /* now get the right place in the ground_frac tables  CK */
      dummy = t_e / TMIN - 1.;
      ilow = dummy;		/* have now truncated to integer below */
      ihi = ilow + 1;		/*these are the indeces bracketing the true value */
      interpfrac = (dummy - ilow);	/*this is the interpolation fraction */
      if (ilow < 0)
	{
	  ilow = 0;
	  ihi = 0;
	  interpfrac = 1.;
	}
      if (ihi > 19)
	{
	  ilow = 19;
	  ihi = 19;
	  interpfrac = 1.;
	}

    }

  else
    {
      Error_silent ("sim_pl: t_r too low www %8.2e t_e %8.2e  t_r %8.2e \n",
		    www, t_e, t_r);
      fudge = 0.0;
      interpfrac = 0.0;
      ihi = ilow = 0;
    }
  if (fudge < MIN_FUDGE || MAX_FUDGE < fudge)
    {
      Error_silent
	("sim_pl: fudge>10 www %8.2e t_e %8.2e  t_r %8.2e \n", www, t_e, t_r);
    }

  /* Initialization of fudges complete.  -- ksl - this comment about the fudge was copied from lucy_mazaali1 where there were two parts to
   * the fudge factor calcuation.   Not clear what you mean here 
   * open lots of files to log the sim factor */

  if (diag_on_off)
    {
      fudgefile = fopen ("fudge_summary.out", "a");
      numfile = fopen ("num_summary.out", "a");
      denomfile = fopen ("denom_summary.out", "a");
      ionfile = fopen ("ion_summary.out", "a");
      fprintf (fudgefile, "Te= %e, Element %s", t_e, ele[nelem].name);
      fprintf (numfile, "Te= %e, Element %s", t_e, ele[nelem].name);
      fprintf (denomfile, "Te= %e, Element %s", t_e, ele[nelem].name);
      fprintf (ionfile, "Te= %e, Element %s", t_e, ele[nelem].name);
    }



  first = ele[nelem].firstion;	/*identify the position of the first and last ion in the array */
  last = first + (ele[nelem].nions);	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */


/* 111219 - ksl - I'm nost sure whay the next for loop is necessary, we are not going to change densit */

  for (nion = ele[nelem].firstion;
       nion < (ele[nelem].firstion + ele[nelem].nions); nion++)
    {
      ion_store[nion] = density[nion];
    }

/* 111001 - ksl - As these next two loops were written, there was nothing keeping them from 
 * runnng over the ends of the arrays.  This is obviously a problem.  It is also not obvious
 * why the factor of 1.1 comes in.  The densities are being set to DENSITY_MIN (which is defined
 * in atomic.h  ????  The code is similar to that in saha.c, which seems similarly silly.  The
 * reason I ran into an error here is that for some reason here all the densities are zero and so
 * it does not know how to stop.  This has to do with the torus implementation.  
 */

  while (density[first] < 1.1 * DENSITY_MIN)
    {
      newden[first] = DENSITY_MIN;
      fudge_store[first] = -1.;
      num_store[first] = -1.;
      denom_store[first] = -1.;
      first++;
      if (first == last)
	{
	  Error
	    ("sim_pl: All of the densities for element %d were below the minimim density\n",
	     nelem);
	  return (1);
	}

    }

  while (density[last - 1] < 1.1 * DENSITY_MIN)
    {
      newden[last - 1] = DENSITY_MIN;
      fudge_store[last - 1] = -1;
      num_store[last - 1] = -1;
      denom_store[last - 1] = -1;
      last--;
    }

  /* At this point we have located the first and last ion we think is important enough to worry about */

  max_ratio = ((ele[nelem].abun * nh) / DENSITY_MIN) / 10.0;	//This is the maximum ratio between two ions of the same element
  sum = newden[first] = 1.;
  fudge_store[first] = 0.0;

  for (nion = first + 1; nion < last; nion++)
    {
//      printf ("Working on ion %i, with density %e\n",nion,density[nion]); 
      current_ratio = (density[nion] / density[nion - 1]);

      /* Calculate the correction factor for this particular ion for the power law approximatin */
      /* THIS IS THE MAIN DIFFERENCE BETWEEN THIS AND LUCY MAZZALI */
      fudge = (xinteg_sim (t_e, f1, f2, nion - 1, max_ratio, current_ratio));

      fudge_store[nion] = fudge;

      /* ksl - not sure what is different about this than the value in the for loop above; density has not been changed I hope */
      ion_store[nion] = density[nion];

      numerator = newden[nion - 1] * fudge * (ne) * density[nion];
      denominator = density[nion - 1] * xne;
      q = numerator / denominator;

      /* find fraction of recombinations going directly to
         the ground state for this temperature and ion */
//old this is the standard calculation of fudge2 (zeta)     fudge2 =
//      ground_frac[nion - 1].frac[ilow] +
//      interpfrac * (ground_frac[nion - 1].frac[ihi] -
//                    ground_frac[nion - 1].frac[ilow]); 

/* nsh 110826 This is the new call to compute_zeta, which is a function which incorporates 
the calculation of the recombination to ground state plus the dielectronic recombination. 
Mode 1 is just recomb to gs by the old method. Mode 2 is recomb to gs plus DR correction. */

//!!! ksl Eliminate dielectronic combination for now      fudge2 = compute_zeta (t_e, nion, ilow, ihi, interpfrac, f1, f2, 2);
      fudge2 = 1.;


/*  From lucymazzili : Since nion-1 points at state i-1 (saha: n_i/n_i-1) we want ground_frac[nion-1].
         Note that we NEVER access the last ion for any element in that way
         which is just as well since you can't have recombinations INTO
         a fully ionized state -- The corresponding lines in recombination.out 
         are just dummies to make reading the data easier */

      newden[nion] = fudge2 * q;	//we multiply by the recombination fraction going directly to ground state
      sum += newden[nion];
//   From lucymazzili : This is the equation being calculated
//              sum+=newden[nion]=newden[nion-1]*fudge*(*ne)*density[nion]/density[nion-1]/xne;
    }

  /* Print out the accumulated diagnostic information */
  if (diag_on_off)
    {

      for (nion = ele[nelem].firstion;
	   nion < (ele[nelem].firstion + ele[nelem].nions); nion++)
	{
	  fprintf (fudgefile, ",%e", fudge_store[nion]);
	  fprintf (numfile, ",%e", num_store[nion]);
	  fprintf (denomfile, ",%e", denom_store[nion]);
	  fprintf (ionfile, ",%e", ion_store[nion]);
	}
      fprintf (fudgefile, "\n");
      fprintf (numfile, "\n");
      fprintf (denomfile, "\n");
      fprintf (ionfile, "\n");

      fclose (fudgefile);
      fclose (numfile);
      fclose (denomfile);
      fclose (ionfile);
    }
  /* End of the diagnostic information */


  /* Renormalize the abundances so that we do not end up with too much of the element */
  a = nh * ele[nelem].abun / sum;
  for (nion = first; nion < last; nion++)
    newden[nion] *= a;


  return (0);
}





/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  This is where the fudge factor in the power law case is calculated

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111217	ksl	Added standard documentation header to encourage more systematic documentation
111229	ksl	Minor modifications reflecting the fact that nxbands, xfreq has been moved into
		goe structure.  Added a counter to turn off the error message about missing
		photoionization data.

 ************************************************************************/

int xinteg_sim_err=0;

double
xinteg_sim (t, f1, f2, nion, max_ratio, current_ratio)
     double t;			// The temperature at which to integrate the cross sections
     double f1, f2, max_ratio, current_ratio;	// The frequencies overwhich to integrate the cross sections
     int nion;			// The ion for which the integrated cross section is calculated
{
  int n, j;
  double sim_num, sim_denom, sim_frac;	// The upper and lower part of the fraction used to recalculate the ion fractions
  double fthresh, fmax;
//  double den_config ();
//  double sigma_phot (), sigma_phot_topbase ();
  int ntmin, ntmax;		// These are the Topbase photo-ionization levels that are used
  int nvmin, nvmax;		// These are the limits on the Verland x-sections
  double qromb ();
//  double fb_planck (), verner_planck ();


  if (-1 < nion && nion < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[nion].ntop_first;
      ntmax = ntmin + ion[nion].ntop;
      nvmin = nion;
      nvmax = nvmin + 1;
//;     printf ("For ion %i, the first topbase level is %i, and the first verner level is %i\n",nion,ion[nion].ntop_first,nvmin);
    }
  else				// Get the total emissivity
    {
      Error ("xinteg_sim: %d is unacceptable value of nion\n", nion);
      mytrap ();
      exit (0);
      return (0);
    }

// Put information where it can be used by the integrating function
  sim_te = t;

// Place over all limits on the integration interval if they are very large
/* We need to limit the frequency range to one that is reasonable
if we are going to integrate */
/* Because either f2 corresponded to something redward of 1000 A or f1 
was blueward of 10 Angstroms 

Note - ksl - It is not clear why this is necessary, we have already defined the frequency limits
*/

  if (f1 < 3e14)
    f1 = 3e14;
  if (f2 > 3e19)
    f2 = 3e19;
  if (f2 < f1)
    return (0);

  sim_num = 0.0;
  sim_denom = 0.0;

/* First work on ions with topbase data */
  if (ion[nion].phot_info == 1)
    {
      n = ntmin;		// we only want to use the ground state topbase ionisation cross section
      sim_xtop = &phot_top[n];
      /* Adding an if statement here so that photoionization that's part of a macro atom is 
         not included here (these will be dealt with elsewhere). (SS, Apr04) */
      if (sim_xtop->macro_info == 0 || geo.macro_simple == 1 || geo.rt_mode == 1)	//Macro atom check. (SS)
	{
	  fthresh = sim_xtop->freq[0];
	  fmax = sim_xtop->freq[sim_xtop->np - 1];	// Argues that this should be part of structure
	  if (f1 > fthresh)
	    fthresh = f1;
	  if (f2 < fmax)
	    fmax = f2;


	  if (fmax > fthresh)
	    {
//      printf ("We are going to topbase qromb for ion %i with fthresh=%e and fmax=%e\n",n,fthresh,fmax);
// old -we are now going to split this up                               sim_num = qromb (tb_pow, fthresh, fmax, 1.e-4);
/* for the power law, we need to carry out several seperate integrations since there is 
 * no stricture that the power law description of the spectrum is continuous over the 
 * frequency boundaries. This seems to cause serious problems in qromb which hangs the code. 
 * We must loop over the frequency intervals 
 */

/* There are five cases we need to consider.

1: If the current band contains the threshold frequency and the max frequency, we set 
alpha and w to the values for this band, and integreate from fthresh to fmax.

2: If the threshold frequency is in the band, but the max frequency isn't, we need to 
integrate from fthresh to the maximum frequency in the band

3: If the threshold frequency isnt in the band, but the max frequency is, we need to 
integrate from the minimum frequency in the band to fmax

4: If fthresh is less than the min, and fmax is greater than the max for a band, we need to 
integreate over the whole band. 

5: If the maximum frequency for this band is less than the threshold frequency, or the minuum 
frequency for the band is greater than the maximum frequency, we dont need to do anything

*/

	      /*71 - 111229 - ksl - Small changes reflecting the move of nxfreq and xfreq into geo structure */

	      for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
		{
		  if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//Case 1
		    {
		      xsim_alpha = xxxplasma->pl_alpha[j];
		      xsim_w = xxxplasma->pl_w[j];
		      sim_num += qromb (tb_pow, fthresh, fmax, 1.e-4);
		    }
		  else if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j + 1] < fmax)	//case 2 
		    {
		      xsim_alpha = xxxplasma->pl_alpha[j];
		      xsim_w = xxxplasma->pl_w[j];
		      sim_num +=
			qromb (tb_pow, fthresh, geo.xfreq[j + 1], 1.e-4);
		    }
		  else if (geo.xfreq[j] > fthresh && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//case 3
		    {
		      xsim_alpha = xxxplasma->pl_alpha[j];
		      xsim_w = xxxplasma->pl_w[j];
		      sim_num += qromb (tb_pow, geo.xfreq[j], fmax, 1.e-4);
		    }
		  else if (geo.xfreq[j] > fthresh && geo.xfreq[j + 1] < fmax)	// case 4
		    {
		      xsim_alpha = xxxplasma->pl_alpha[j];
		      xsim_w = xxxplasma->pl_w[j];
		      sim_num +=
			qromb (tb_pow, geo.xfreq[j], geo.xfreq[j + 1], 1.e-4);
		    }
		  else		//case 5 - should only be the case where the band is outside the range for the integral.
		    {
		      sim_num += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		    }
		}
	      sim_denom = qromb (tb_planck, fthresh, fmax, 1.e-4);
	    }


//      printf ("And we are back in the room with PL=%e and Planck=%e\n",sim_num,sim_denom);

	}
    }


/* This completes the calculation of those levels for which we have Topbase x-sections. Now proces ions that use Verner data */

  else if (ion[nvmin].phot_info == 0)
    {			
      n = nvmin;		//just the ground state ioinzation fraction.
      sim_xver = &xphot[ion[n].nxphot];
      fthresh = sim_xver->freq_t;
//      printf("XXXXXXXXX fthresh=%e\n",fthresh);

      fmax = sim_xver->freq_max;	//So at this point these are the maximal allowable
      if (f1 > fthresh)
	fthresh = f1;		//So move fthresh up, if f1 was greater than this
      if (f2 < fmax)		//Move fmax down if we wanted a narrower range
	fmax = f2;
//      printf("We want to do verner, with fmax=%e and fthresh=%e",fmax,fthresh);

// Now integrate only if its in allowable range  && there are ions to recombine
      if (fmax > fthresh)
//      printf ("We are going to verner qromb with ion %i, fthresh=%e and fmax=%e\n",n,fthresh,fmax);
	{
	  sim_denom = qromb (verner_planck, fthresh, fmax, 1.e-4);
	  /* for the power law, we need to carry out several seperate integrations since there is no stricture that the power law description of the spectrum is continuous over the frequency boundaries. This seems to cause seriousw problems in qromb which hangs the code. We must loop over the drequency intervals */


/* There are five cases we need to consider.

1: If the current band contains the threshold frequency and the max frequency, we set 
alpha and w to the values for this band, and integreate from fthresh to fmax.

2: If the threshold frequency is in the band, but the max frequency isn't, we need to 
integrate from fthresh to the maximum frequency in the band

3: If the threshold frequency isnt in the band, but the max frequency is, we need to 
integrate from the minimum frequency in the band to fmax

4: If fthresh is less than the min, and fmax is greater than the max for a band, we need to 
integreate over the whole band. 

5: If the maximum frequency for this band is less than the threshold frequency, or the minuum 
frequency for the band is greater than the maximum frequency, we dont need to do anything

*/

	  for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	    {
	      if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//Case 1
		{
		  xsim_alpha = xxxplasma->pl_alpha[j];
		  xsim_w = xxxplasma->pl_w[j];
		  sim_num += qromb (verner_pow, fthresh, fmax, 1.e-4);
		}
	      else if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j + 1] < fmax)	//case 2 
		{
		  xsim_alpha = xxxplasma->pl_alpha[j];
		  xsim_w = xxxplasma->pl_w[j];
		  sim_num +=
		    qromb (verner_pow, fthresh, geo.xfreq[j + 1], 1.e-4);
		}
	      else if (geo.xfreq[j] > fthresh && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//case 3
		{
		  xsim_alpha = xxxplasma->pl_alpha[j];
		  xsim_w = xxxplasma->pl_w[j];
		  sim_num += qromb (verner_pow, geo.xfreq[j], fmax, 1.e-4);
		}
	      else if (geo.xfreq[j] > fthresh && geo.xfreq[j + 1] < fmax)	// case 4
		{
		  xsim_alpha = xxxplasma->pl_alpha[j];
		  xsim_w = xxxplasma->pl_w[j];
		  sim_num +=
		    qromb (verner_pow, geo.xfreq[j], geo.xfreq[j + 1], 1.e-4);
		}
	      else		//case 5 - should only be the case where the band is outside the range for the integral.
		{
		  sim_num += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		}
	    }
	}
    }
  else
    {
	    xinteg_sim_err++;
	    /* 71 - 111229 - ksl - Suppress this error after 100 instances so program does not bomb */
	    if (xinteg_sim_err<100){
      Error
	("xinteg_sim: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	 nion,ion[nion].z,ion[nion].istate);
	    }
	    else if (xinteg_sim_err==100){
		    Error("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	    }

	   
      sim_frac = 1.;
      return (sim_frac);
    }

  /* End of main loop */

  if (sim_denom < ((sim_num * current_ratio) / max_ratio))
    {
      Log_silent ("xinteg_sim: Sim denom of %e is too low, resetting to %e\n",
		  sim_denom, ((sim_num * current_ratio) / max_ratio));
      sim_denom = ((sim_num * current_ratio) / max_ratio);
    }

  num_store[nion + 1] = sim_num;
  denom_store[nion + 1] = sim_denom;

  sim_frac = sim_num / sim_denom;


  return (sim_frac);
}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

 tb_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
   ionisation cross section is significant. This is the function for ions with a topbase cross section NSH 16/12/10 


  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111217	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/


double
tb_planck (freq)
     double freq;
{
  double answer, bbe;

  bbe = exp ((H * freq) / (BOLTZMANN * sim_te));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot_topbase (sim_xtop, freq);
  answer /= freq;

  return (answer);
}




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
	verner_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
	ionisation cross section is significant. This is the function for ions with a verner cross section NSH 16/12/10 

  Arguments:  


  Returns:

  Notes:


 

  History:
111217	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/


double
verner_planck (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * sim_te));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot (sim_xver, freq);
  answer /= freq;

  return (answer);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
	tb_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation 
	cross section
  	to allow the numerator of the sim correction factor to be calulated for ions with topbase cross sections.

  Arguments:  


  Returns:

  Notes:


 

  History:
111217	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/



double
tb_pow (freq)
     double freq;
{
  double answer;

  answer = xsim_w * (pow (freq, (xsim_alpha - 1.0)));
  answer *= sigma_phot_topbase (sim_xtop, freq);	// and finally multiply by the cross section.

  return (answer);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	verner_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation cross section
   	to allow the numerator of the sim correction factor to be calulated for ions with verner cross sections.

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111227	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/


double
verner_pow (freq)
     double freq;
{
  double answer;

  answer = xsim_w * (pow (freq, (xsim_alpha - 1.0)));
  answer *= sigma_phot (sim_xver, freq);	// and finally multiply by the cross section.

  return (answer);
}




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
 	This little routine solves equation 18 of Sim et at (2008) to get an estimate of the spectral index of a cell
 	Equation 18 seems to be a monotonically increasing function in alpha, so this should work (famous last words)

  Description:	

  Arguments:  


  Returns:

  Notes:
  	Rather than use a solver this is just done by iteration

	The routine is only called by balance_sub.c.  It may not be part of python

	In power.c, alpha is caluculated using zbrent and sim_alpha_func, although
	ksl will ultimately renamve this to something else, e.g mean_freq_to_alpha_func
	which would be more descriptive.


 

  History:
NSH 7/2/11	Coded
111227	ksl	Added standard documentation header to encourage more systematic documentation

 ************************************************************************/


double
sim_alphasolve (ratans, numin, numax)
     double ratans, numin, numax;
{
  double alphamin, alphamax, alphamid;
  double error, ratmid, alpha;


// Set crazy values of range of alpha - if this screws up I might need to put in a check here
  alphamin = -10;
  alphamax = 10;
  alpha = -0.2;
// Set a value of error that gets the loop running
  error = 10.;
  printf ("Solving with ratans=%e, numin=%e, numax=%e\n", ratans, numin,
	  numax);
  printf ("Just out of interest, for alpha=-0.2, we get %e\n",
	  ((alpha + 1) / (alpha + 2)) * ((pow (numax, (alpha + 2)) -
					  pow (numin,
					       (alpha + 2))) / (pow (numax,
								     (alpha +
								      1)) -
								pow (numin,
								     (alpha +
								      1)))));
  while (error > 0.001)
    {
// Work out a mid point
      alphamid = (alphamax + alphamin) / 2.0;
// What En(2)/En(1) is implied by this alpha
      ratmid =
	((alphamid + 1) / (alphamid + 2)) * ((pow (numax, (alphamid + 2)) -
					      pow (numin,
						   (alphamid +
						    2))) / (pow (numax,
								 (alphamid +
								  1)) -
							    pow (numin,
								 (alphamid +
								  1))));
// If actual answer is below ratmid, our solution is in the lower half
      if (ratans < ratmid)
	alphamax = alphamid;
      else
// If actual answer is above ratmin, our solution is in the upper half
	alphamin = alphamid;
      error = fabs (alphamax - alphamin);
      printf
	("Searching for answer, alphamin=%f, alphamax=%f, error=%f, current answer=%e \n",
	 alphamin, alphamax, error, ratmid);
    }
// If we get here, we must have converged.
  return (alphamid);
}





