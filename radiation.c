
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	 radiation(p,ds) updates the radiation field parameters in the wind and reduces 
	 the weight of the photon as a result of the effect free free and photoabsorption.
	 radiation also keeps track of the number of photoionizations of h and he in the
	 cell. 
Arguments:

	PhotPtr p;	the photon
	double ds	the distance the photon has travelled in the cell

Returns:
	Always returns 0.  The pieces of the wind structure which are updated are
	j,ave_freq,ntot,  heat_photo,heat_ff,heat_h,heat_he1,heat_he2, heat_z,
	nioniz, and ioniz[].
 
Description:	 
Notes:
	
	The # of ionizations of a specific ion = (w(0)-w(s))*n_i sigma_i/ (h nu * kappa_tot).  (The # of ionizations
	is just given by the integral of n_i sigma_i w(s) / h nu ds, but if the density is assumed to
	be constant and sigma is also constant [therefy ignoring velocity shifts in a grid cell], then
	n_i sigma and h nu can be brought outside the integral.  Hence the integral is just over w(s),
	but that is just given by (w(0)-w(smax))/kappa_tot.)  The routine calculates the number of ionizations per
	unit volume.
	
	This routine is very time counsuming for our normal file which has a large number of x-sections.  The problem
	is going through all the do loops, not insofar as I can determine anything else.  Trying to reduce the calucatiions
	by using a density criterion is a amall effect.  One really needs to avoid the calculations, either by avoiding
	the do loops altogether or by reducing the number of input levels.  It's possible that if one were clever about
	the thresholds (as we are on the lines that one could figure out a wininning strategy as it is all brute force
        do loops..  

	57h -- This routine was very time consuming.  The time spent is/was well out of proportion to
	the affect free-bound processes have on the overall spectrum, and so signifincant changes have been made
	to the earlier versions of the routine.  The fundamental problem before 57h was that every time one
	entered the routine (which was every for every movement of a photon) in the code.  Basically there were
	3  changes made:
		1. During the detailed spectrum cycle, the code avoids the portions of the routine that are
	only needed during the ionization cycles.
		2. Switches have been installed tha skip the free-bound section altogether if PHOT_DEN_MIN is
	less than zero.  It is very reasonable to do this for the detailed spectrum calculation if one is
	doing a standard 
	were to skip portions of when they were not needed 

	?? Would it be more natural to include electron scattering here in Radiation as well ??
	?? Radiation needs a major overhaul.  A substantial portion of this routine is not needed in the extract 
	?? portion of the calculation.  In addition the do loops go through all of the ions checking one at a time
	?? whether they are above the frequencey threshold.  
	?? The solution I believe is to include some kind of switch that tells the routine when one is doing
	?? the ionization calculation and to skip the unnecessary sections when extract is being carried out.
	?? In addition, if there were a set of ptrs to the photionization structures that were orded by frequency,
	?? similar to the line list, one could then change to loops so that one only had to check up to the
	?? first x-section that had a threshold up to the photon frequency, and not all of the rest.
	?? At present, I have simply chopped the photoionizations being considered to include only thresholds
        ?? shortward of the Lyman limit...e.g. 1 Rydberg, but this makes it impossible to discuss the Balmer jump
History:
 	97jan	ksl	Coded as part of python effort
	98jul	ksl	Almost entirely recoded to eliminate arbitrary split between
			the several preexisting routines.
	98nov	ksl	Added back tracking of number of h and he ionizations
	99jan	ksl	Increased COLMIN from 1.e6 to 1.e10 and added checks so
			that one does not attempt to calculate photoionization
			cross-sections below threshold.  Both of these changes
			are intended to speed this routine.
        01oct   ksl     Modifications begun to incorporate Topbase photoionization
                        x-sections.
        01oct   ksl     Moved fb_cooling to balance_abso.  It's only used by
                        balance and probably should not be there either.
	02jun	ksl	Improved/fixed treatment of calculation of number of ionizations.
	04apr	ksl	SS had modified this routine to allow for macro-atoms in python_47, but his modifications
			left very little for radiation to accomplish.  I have returned to the old version of 
			routine, and moved the little bit that needed to be done in this routine for
			the macro approach to the calling routine.  Once we abandon the old approach this
			routine can be deleted.
	04aug	ksl	Fixed an error which had crept into the program between 51 and 51a that caused the
			disk models to be wrong.  The problem was that there are two places where the
			same frequency limit should have been used, but the limits had been set differently.
	04dec	ksl	54a -- Miniscule change to eliminate -03 warnings.  Also eliminate some variables
			that were not really being used.
	06may	ksl	57+ -- Began modifications to allow for splitting the wind and plasma structures
	06jul	ksl	57h -- Made various changes intended to speed up this routine. These include
			skipping sections of the routine in the spectrum generation
			phase that are not used, allowing control over whether fb opacities
			are calculated at all, and using frequency ordered pointer arrays
			to shorten the loops.
	11aug	nsh	70 changes made to radiation to allow compton cooling to be computed
	11aug	nsh	70 Changed printout of spectra in selected regions so it is always
			the midpoint of the wind
	12may	nsh	72 Added induced compton
	12jun 	nsh	72 Added lines to write out photon stats to a file dirung diagnostics. This is
			to allow us to see how well spectra are being modelled by power law W and alpha

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

#define COLMIN	0.01

int iicount = 0;

int
radiation (p, ds)
     PhotPtr p;
     double ds;
{
  PhotoionizationPtr x_ptr;
  TopPhotPtr x_top_ptr;

  WindPtr one;
  PlasmaPtr xplasma;

  double freq;
  double kappa_tot, frac_tot, frac_ff;
  double frac_z, frac_comp;  /* nsh 1108 added frac_comp - the heating in the cell due to compton heating */
  double frac_ind_comp; /* nsh 1205 added frac_ind_comp - the heating due to induced compton heating */
  double kappa_ion[NIONS];
  double frac_ion[NIONS];
  double density, ft, tau, tau2;
  double energy_abs;
  int n, nion, i; /* nsh 1108 i added as counter for banded j and ave_freq */
  double q, x, z;
  double w_ave, w_in, w_out;
  double den_config ();
  int nconf;
  double weight_of_packet, y;
  int ii,jj;  


  one = &wmain[p->grid];	/* So one is the grid cell of interest */
  xplasma = &plasmamain[one->nplasma];
  check_plasma (xplasma, "radiation");
  freq = p->freq;
  
  kappa_tot = frac_ff = kappa_ff (xplasma, freq);	/* Add ff opacity */
  kappa_tot += frac_comp = kappa_comp (xplasma, freq);    /* 70 NSH 1108 calculate compton opacity, store it in kappa_comp and also add it to kappa_tot, the total opacity for the photon path */
  kappa_tot += frac_ind_comp = kappa_ind_comp (xplasma, freq, ds, p->w);
  frac_tot = frac_z = 0;	/* 59a - ksl - Moved this line out of loop to avoid warning, but notes 
				   indicate this is all disagnostic and might be removed */
//	printf ("In radiation we have ds=%e, W=%e, nu=%e\n",ds,p->w,p->freq);

  if (freq > phot_freq_min)

    {

      if (geo.ioniz_or_extract)	// 57h -- ksl -- 060715
	{			// Initialize during ionization cycles only
	  for (nion = 0; nion < nions; nion++)
	    {
	      kappa_ion[nion] = 0;
	      frac_ion[nion] = 0;
	    }
	}
/* Next section is for photoionization with Topbase.  There may be more
than one x-section associated with an ion, and so one has to keep track
of the energy that goes into heating electrons carefully.  */

/* Next steps are a way to avoid the loop over photoionization x sections when it should not matter */
      if (DENSITY_PHOT_MIN > 0)	// 57h -- ksl -- 060715
	{			// Initialize during ionization cycles only


/* 57h -- 06jul -- ksl -- change loop to use pointers ordered by frequency */
	  for (n = 0; n < ntop_phot; n++)
	    {
	      x_top_ptr = phot_top_ptr[n];
	      ft = x_top_ptr->freq[0];
	      if (ft > freq)
		break;		// The remaining transitions will have higher thresholds
	      if (freq < x_top_ptr->freq[x_top_ptr->np - 1])
		{
/*Need the appropriate density at this point. */

		  nconf = x_top_ptr->nlev;
		  density = den_config (xplasma, nconf);
		  if (density > DENSITY_PHOT_MIN)
		    {
		      kappa_tot += x =
			sigma_phot_topbase (x_top_ptr, freq) * density;
//		printf("Freq=%e, lower level=%i, ion=%i sigma=%e\n",freq,nconf,phot_top_ptr[n]->z,sigma_phot_topbase (x_top_ptr, freq));
/* I believe most of next steps are totally diagnsitic; it is possible if 
statement could be deleted entirely 060802 -- ksl */
		      if (geo.ioniz_or_extract)	// 57h -- ksl -- 060715
			{	// Calculate during ionization cycles only
			  frac_tot += z = x * (freq - ft) / freq;
			  nion = config[nconf].nion;
			  if (nion > 3)
			    {
			      frac_z += z;
			    }
			  frac_ion[nion] += z;
			  kappa_ion[nion] += x;
			}

		    }


		}
					
	    }

// Next section is for photoionization of those ions using VFKY values

	  for (n = 0; n < nxphot; n++)
	    {
	      x_ptr = xphot_ptr[n];
	      nion = x_ptr->nion;

	      if (ion[nion].phot_info == 0)
		{		// Avoid ions with topbase x-sections
		  ft = x_ptr->freq_t;
		  if (ft > freq)
		    break;	// The remaining transitions will have higher thresholds

		  density =
		    xplasma->density[nion] * ion[nion].g /
		    xplasma->partition[nion];

		  if (density > DENSITY_PHOT_MIN)
		    {
		      kappa_tot += x = sigma_phot (x_ptr, freq) * density;

/* Nxt if statment downs to kappa_ion appers to be totally diagnositc -060802 -- ksl */
		      if (geo.ioniz_or_extract)	// 57h -- ksl -- 060715
			{	// Calculate during ionization cycles only
			  frac_tot += z = x * (freq - ft) / freq;

			  if (nion > 3)
			    {
			      frac_z += z;
			    }
			  frac_ion[nion] += z;
			  kappa_ion[nion] += x;
			}
		    }
		}
	    }
	}

    }
  tau = kappa_tot * ds;
  w_in = p->w;

if (sane_check(tau)) {
	Error("Radiation:sane_check CHECKING ff=%e, comp=%e, ind_comp=%e\n",frac_ff,frac_comp,frac_ind_comp);
}
/* Calculate the reduction in the w of the photon bundle along with the average
   weight in the cell */

  if (tau > 0.001)
    {				/* Need differentiate between thick and thin cases */
      x = exp (-tau);
      p->w = w_out = w_in * x;
      energy_abs = w_in - w_out;
      w_ave = (w_in - w_out) / tau;
    }
  else
    {
      tau2 = tau * tau;
      p->w = w_out = w_in * (1. - tau + 0.5 * tau2);	/*Calculate to second order */
      energy_abs = w_in * (tau - 0.5 * tau2);
      w_ave = w_in * (1. - 0.5 * tau + 0.1666667 * tau2);
    }

  /*74a_ksl: 121215 -- Added to check on a problem photon */
  if (sane_check(p->w)){
		  Error("Radiation:sane_check photon weight is %e for tau %e\n",p->w,tau);
		  }

  if (geo.ioniz_or_extract == 0)
    return (0);			// 57h -- ksl -- 060715

/* Everything after this is only needed for ionization calculations */
/* Update the radiation parameters used ultimately in calculating t_r */

  xplasma->ntot++;


  if (HEV*p->freq > 13.6)
	{
  	xplasma->ip+=((w_ave*ds)/(H*p->freq)); 
	/* 70h -- nsh -- 111004 added to try to calculate the IP for the cell. Note that 
	 * this may well end up not being correct, since the same photon could be counted 
	 * several times if it is rattling around.... */
	}

/* NSH 15/4/11 Lines added to try to keep track of where the photons are coming from, 
 * and hence get an idea of how 'agny' or 'disky' the cell is. */
/* ksl - 1112 - Fixed this so it records the number of photon bundles and not the total
 * number of photons.  Also used the PTYPE designations as one should as a matter of 
 * course
 */



//OLD71  if (p->origin == 0)
//OLD71	xplasma->ntot_star+=(w_in/(H*p->freq));
//OLD71  else if (p->origin == 1)
//OLD71	xplasma->ntot_bl+=(w_in/(H*p->freq));
//OLD71  else if (p->origin == 2)
//OLD71	xplasma->ntot_disk+=(w_in/(H*p->freq));
//OLD71  else if (p->origin == 3)
//OLD71	xplasma->ntot_wind+=(w_in/(H*p->freq));
//OLD71  else if (p->origin == 4)
//OLD71	xplasma->ntot_agn+=(w_in/(H*p->freq));


  if (p->origin == PTYPE_STAR)
	xplasma->ntot_star++;
  else if (p->origin == PTYPE_BL)
	xplasma->ntot_bl++;
  else if (p->origin == PTYPE_DISK)
	xplasma->ntot_disk++;
  else if (p->origin == PTYPE_WIND)
	xplasma->ntot_wind++;
  else if (p->origin == PTYPE_AGN)
	xplasma->ntot_agn++;


/*photon weight times distance in the shell is proportional to the mean intensity */
  xplasma->j += w_ave * ds;

/* frequency weighted by the weights and distance       in the shell .  See eqn 2 ML93 */
  xplasma->mean_ds += ds;
  xplasma->n_ds ++;
  xplasma->ave_freq += p->freq * w_ave * ds;

  if (p->freq > xplasma->max_freq)
	xplasma->max_freq = p->freq;


//wind_n_to_if(one->nwind,&ii,&jj);  ??? Not complete. Intended to allow more flexible tracking of photon spectra  ksl 1108
//if ii==ndim/2:  // 110804 - ksl - Adapt to print out a column in the middle no matter what the dimensions
//	if (one->nwind > 59 && one->nwind < 90)
//  1.75e16 <     (2*sqrt(one->xcen[0]*one->xcen[0]+one->xcen[1]*one->xcen[1]) - sqrt(one->x[0]*one->x[0]+one->x[1]*one->x[1]))) 
//	if (geo.wind_type == 9)
//	if (one->nwind==22)
//printf ("PHOTON_DETAILS %3d %3d %3d %8.3e %8.3e %8.3e cell%3d wind cell%3d %e %e %e\n",geo.wcycle,ii,jj,p->freq,w_ave,ds,one->nplasma,one->nwind,p->w,w_in,tau);

  if (diag_on_off == 1 && ncstat > 0)
	{
	for (i=0;i<ncstat;i++)
		{
		if (one->nplasma==ncell_stats[i])
			{
			fprintf (pstatptr,"PHOTON_DETAILS %3d %3d %3d %8.3e %8.3e %8.3e cell%3d wind cell%3d\n",geo.wcycle,ii,jj,p->freq,w_ave,ds,one->nplasma,one->nwind);
			}
		}
	}


/* 1108 NSH  THe next loop updates the banded versions of j and ave_freq, note that the lines above still update 
 * the parameters for all the photons produced. We are simply controlling which photons are used for the sim powewr 
 * law stuff by banding.
 * nxfreq refers to how many frequencies we have defining the bands. So, if we have 5 bands, we have 6 frequencies, 
 * going from xfreq[0] to xfreq[5] Since we are breaking out of the loop when i>=nxfreq, this means the last loop 
 * will be i=nxfreq-1*/

/* 71 - 111229 - ksl - modified to reflect fact that I have moved nxbands and xfreq into the geo structure */

  for (i=0 ; i<geo.nxfreq ; i++)  
	{
	if (geo.xfreq[i] < p->freq && p->freq <= geo.xfreq[i+1])
		{
		xplasma->xave_freq[i] += p->freq * w_ave * ds;  /*1108 NSH/KSL frequency weighted by weight and distance */
		xplasma->xsd_freq[i] += p->freq * p->freq * w_ave * ds; /*1208 NSH imput to allow standard deviation to be calculated */
		xplasma->xj[i] += w_ave * ds;  			/*1108 NSH/KSL photon weight times distance travelled */
		xplasma->nxtot[i] ++; 				/*1108 NSH increment the frequency banded photon counter */

		}
	}



  if (sane_check (xplasma->j) || sane_check (xplasma->ave_freq))
    {
      Error ("radiation:sane_check Problem with j %g or ave_freq %g\n", xplasma->j,
	     xplasma->ave_freq);
    }


  if (kappa_tot > 0)
    {
      //If statement added 01mar18 ksl to correct problem of zero divide
      //  in odd situations where no continuum opacity
      z = (energy_abs) / kappa_tot;
      Log("!!JMHEAT: heat_tot %8.2e  heat_photo %8.2e", xplasma->heat_tot, xplasma->heat_ff );
      xplasma->heat_ff += z * frac_ff;
      xplasma->heat_tot += z * frac_ff;
      xplasma->heat_comp += z * frac_comp; /* NSH 1108 Calculate the heating in the cell due to compton heating */
      xplasma->heat_tot += z * frac_comp;  /* NSH 1108 Add the compton heating to the total heating for the cell */
      xplasma->heat_tot += z * frac_ind_comp; /* NSH 1205 Calculate the heating in the celldue to induced compton heating */
      xplasma->heat_ind_comp += z * frac_ind_comp; /* NSH 1205 Increment the induced compton heating counter for the cell */
      if (freq > phot_freq_min)
	//
//      if (freq > (CR / 100.)) //modified CR to CR/100 - SS June 04
	/* 04aug -- ksl -- 52.  Using CR/100 can speed the program up
	 * somewhat, but the limit here needs to be the same as the
	 * one above.  Setting the two limits differently can cause
	 * unpredictable and serious errors.
	 */
	{
	  xplasma->heat_photo += z * frac_tot;
	  Log("!!JMHEAT: heat_tot %8.2e  heat_photo %8.2e", xplasma->heat_tot, xplasma->heat_ff );
	  xplasma->heat_z += z * frac_z;
	  xplasma->heat_tot += z * frac_tot;	//All of the photoinization opacities
/* Calculate the number of photoionizations per unit volume for H and He */
	  xplasma->nioniz++;
	  q = (z) / (H * freq * one->vol);
/* So xplasma->ioniz for each species is just 
     (energy_abs)*kappa_h/kappa_tot / H*freq / volume
    or the number of photons absorbed in this bundle per unit volume by this ion
*/

	  for (nion = 0; nion < nions; nion++)
	    {
	      xplasma->ioniz[nion] += kappa_ion[nion] * q;
	      xplasma->heat_ion[nion] += frac_ion[nion] * z;
	    }

	}
    }

  /* Now for contribution to inner shell ionization estimators (SS, Dec 08) */
  for (n = 0; n < nauger; n++)
    {
      ft = augerion[n].freq_t;
      //printf("Auger tests: %g %g %g\n", augerion[n].freq_t, freq, p->freq);
      if (p->freq > ft)
 	{
 	  //	  printf("Adding a packet to AUGER via radiation %g \n", freq);
	  
 	  weight_of_packet = w_ave;
 	  x = sigma_phot_verner(&augerion[n], freq); //this is the cross section
 	  y = weight_of_packet * x * ds;
	  
 	  xplasma->gamma_inshl[n] += y / (freq * H * one->vol);
 	}
    }
  
  
  
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: 
	double kappa_ff(w,freq) calculates the free-free opacity allowing for stimulated emission
 	
Arguments:

Returns:
	
Description:	 
	Formula from Allen
	Includes only singly ionized H and doubly ionized he 	

Notes:

History:
	98aug	ksl	Coded as part of python effort
        04Apr   SS      If statement added for cases with hydrogen only.
                        Note: this routine assumes that H I is the first ion
                        and that He II is the fourth ion.
	05may	ksl	57+ -- Modified to eliminate WindPtr in favor
			of PlasmaPtr
   	12oct	nsh	73 -- Modified to use a prefector summed over all ions, calculated prior
			to the photon flight

**************************************************************/


double
kappa_ff (xplasma, freq)
     PlasmaPtr xplasma;
     double freq;
{
  double x;
  double exp ();
  double x1,x2, x3;

if (gaunt_n_gsqrd==0) //Maintain old behaviour
	{
  	if (nelements > 1)
    		{
      		x = x1 =
		3.692e8 * xplasma->ne * (xplasma->density[1] +
				 4. * xplasma->density[4]);
    		}
  	else
    		{
      		x = x1 = 3.692e8 * xplasma->ne * (xplasma->density[1]);
    		}
	}
else
	{
   	x=x1=xplasma->kappa_ff_factor;
	}
  x *= x2 = (1. - exp (-H_OVER_K * freq / xplasma->t_e));
  x /= x3 = (sqrt (xplasma->t_e) * freq * freq * freq);

  return (x);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis: 
	double sigma_phot(x_ptr,freq)	calculates the photionization crossection due to the transition 
	associated with x_ptr at frequency freq
Arguments:

Returns:
	
Description:	 
	sigma_phot uses Verner et al.'s interpolation formulae for the photoionization crossection
	to calculate the bound free (or photoionization) optical depth.  The data must
	have been into the photoionization structures xphot with get_atomic_data and
	the densities of individual ions must have been calculated previously.

Notes:

History:
	98jul	ksl	Coded (actually moved from a subroutine called kappa_ds)
	06jul	ksl	57h -- Added code to avoid recalculating crossection
			if that is possible, i.e. if we are calculating
			a the x-section at the same frequency as last
			time the routine was entired for this ion. 

**************************************************************/

double
sigma_phot (x_ptr, freq)
     struct photoionization *x_ptr;
     double freq;
{
  double ft;
  double x, y;
  double f1, f2, f3;
  double xsection;


  ft = x_ptr->freq_t;		/* threshold frequency */


  if (ft < freq && freq < x_ptr->freq_max)
    {
      if (freq == x_ptr->f_last)
	return (x_ptr->sigma_last);	// Avoid recalculating xsection

      x = freq / x_ptr->freq0 - x_ptr->y0;
      y = sqrt (x * x + x_ptr->y1 * x_ptr->y1);

      /* This was line fixed by CK in 1998 Jul */
      f1 = (x - 1.0) * (x - 1.0) + x_ptr->yw * x_ptr->yw;

      f2 = pow (y, 0.5 * x_ptr->p - 5.5);
      f3 = pow (1.0 + sqrt (y / x_ptr->ya), -x_ptr->p);
      xsection = x_ptr->sigma * f1 * f2 * f3;	// the photoinization xsection

/* Store crossesction for future use */
      x_ptr->sigma_last = xsection;
      x_ptr->f_last = freq;

      return (xsection);
    }
  else
    return (0.0);
}


/***********************************************************
				       Space Telescope Science Institute

 Synopsis:
	double sigma_phot_topbase(x_ptr,freq)	calculates the
	photionization crossection due to the transition associated with
	x_ptr at frequency freq
Arguments:
     struct topbase_phot *x_ptr;
     double freq;

Returns:

Description:
	sigma_phot uses the Topbase x-sections to calculate the bound free
	(or photoionization) xsection.	The data must have been into the
	photoionization structures xphot with get_atomic_data and the
	densities of individual ions must have been calculated previously.

Notes:

History:
	01Oct	ksl	Coded as part of general move to use Topbase data
			(especially partial xsections, which did not exist 
			in the Verner et al prescriptions
	02jul	ksl	Fixed error in the way fraction being applied.
			Sigh! and then modified program to use linterp

**************************************************************/

double
sigma_phot_topbase (x_ptr, freq)
     struct topbase_phot *x_ptr;
     double freq;
{
  int nmax;
  double xsection;
  double frac, fbot, ftop;
  int linterp ();
  int nlast;

  

  if (freq < x_ptr->freq[0])
    return (0.0);		// Since this was below threshold
  if (freq == x_ptr->f)
    return (x_ptr->sigma);	// Avoid recalculating xsection

  if (x_ptr->nlast > -1)
    {
      nlast = x_ptr->nlast;
      if ((fbot = x_ptr->freq[nlast]) < freq
	  && freq < (ftop = x_ptr->freq[nlast + 1]))
	{
	  frac = (freq - fbot) / (ftop - fbot);
	  xsection =
	    (1. - frac) * x_ptr->x[nlast] + frac * x_ptr->x[nlast + 1];
	  //Store the results
	  x_ptr->sigma = xsection;
	  x_ptr->f = freq;
	  return (xsection);
	}
    }

/* If got to here, have to go the whole hog in calculating the x-section */
  nmax = x_ptr->np;
  x_ptr->nlast =
    linterp (freq, &x_ptr->freq[0], &x_ptr->x[0], nmax, &xsection);

  //Store the results
  x_ptr->sigma = xsection;
  x_ptr->f = freq;


  return (xsection);

}

/***********************************************************

  Synopsis: 
 	double sigma_phot_verner(x_ptr,freq)	calculates the photionization crossection due to the transition 
 	associated with x_ptr at frequency freq
 Arguments:
 
 Returns:
 
 Description:	 
        Same as sigma_phot but using the older compitation from Verner that includes inner shells
 
 Notes:
 
 History:
 	08dec	SS	Coded (actually mostly copied from sigma_phot)
 
**************************************************************/

double
sigma_phot_verner (x_ptr, freq)
     struct innershell *x_ptr;
     double freq;
{
  double ft;
  double y;
  double f1, f2, f3;
  double xsection;
  
  ft = x_ptr->freq_t;		/* threshold frequency */
  
  if (ft < freq)
    {
      y = freq/x_ptr->E_0*HEV;
      
      f1 = ((y - 1.0) * (y - 1.0)) + (x_ptr->yw * x_ptr->yw);
      f2 = pow (y, 0.5 * x_ptr->P - 5.5 - x_ptr->l);
      f3 = pow (1.0 + sqrt (y / x_ptr->ya), -x_ptr->P);
      xsection = x_ptr->Sigma * f1 * f2 * f3;	// the photoinization xsection
      
      return (xsection);
    }
  else
    return (0.0);
}


/***********************************************************
				       Space Telescope Science Institute

Synopsis:

double den_config(one,nconf) returns the precalculated density
	of a particular "nlte" level.	If levels are not defined for an ion it
	assumes that all ions of are in the ground state.

Arguments:

Returns:

Description: The first case is when the density of the particular level
is in the wind array The second caseis when the density of the levels
for an ion are not tracked, and hence only the total density for the
ion is known.  In that case, we assume the ion is in it's ground state.


Notes:

History:
	01oct	ksl	Coded as part of effort to add topbase
			xsections and detailed balance to python
	05may	ksl	57+ -- Recoded to use PlasmaPtr

**************************************************************/

double
den_config (xplasma, nconf)
     PlasmaPtr xplasma;
     int nconf;
{
  double density;
  int nnlev, nion;

  nnlev = config[nconf].nden;
  nion = config[nconf].nion;

  if (nnlev >= 0)
    {				// Then a "non-lte" level with a density
      density = xplasma->levden[nnlev] * xplasma->density[nion];
    }
  else if (nconf == ion[nion].firstlevel)
    {
/* Then we are using a Topbase photoionization x-section for this ion, 
but we are not storing any densities, and so we assume it is completely in the
in the ground state */
      density = xplasma->density[nion];
    }
  else
    {
      density = 0;
    }

  return (density);


}


/***********************************************************
                Southampton University

Synopsis: pop_kappa_ff_array populates the multiplicative 	
		factor used in the FF calculation. This depends only
		on the densities of ions in the cell, and the electron
		temperature (which feeds into the gaunt factor) so it
		saves time to calculate all this just the once. This
		is called in python.c, just before the photons are 
		sent thruogh the wind.

Arguments:		

Returns:
 
Description:	

Notes:


History:
   12oct           nsh     coded 
   1212	ksl	Added sane check; note that this routine
   		is poorly documented.  Somewhere this 
		should be discribed better.  
 
**************************************************************/


double
pop_kappa_ff_array ()

{

double gsqrd,gaunt,sum;
int i,j;

sum=0.0;
	for (i=0;i<NPLASMA+1;i++)
		{
		for (j=0;j<nions;j++)
			{
			gsqrd=(ion[j].z*ion[j].z*RYD2ERGS)/(BOLTZMANN*plasmamain[i].t_e);
			gaunt=gaunt_ff(gsqrd);
			sum += plasmamain[i].density[j] * ion[j].z * ion[j].z * gaunt;
			/* 74a_ksl  Added to diagnose problem with kappa_ff_fact producing NaN */
			if (sane_check(sum)){
				Error("pop_kappa_ff_array:sane_check sum is %e this is a problem, possible in gaunt %3\n",sum,gaunt);
			}

			}
		plasmamain[i].kappa_ff_factor =  plasmamain[i].ne*sum*3.692e8;

		}

  return (0);
}

