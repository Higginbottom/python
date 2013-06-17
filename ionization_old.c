

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	ionization routines for the wind one cell at a time
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	The intent is that the routine ion_abundances is the steering routine for
	all calculations of the abundances
	
Notes:

History:
	98jul	ksl	Incorporated Christians mods with my comments 
	98nov	ksl	Created the steering routine ion_abundances
	00jan	ksl	Homogenized calls so that ionization_mod_on_the_spot
			and ionization_mod_on_the_spot_exact do not call the
			entire wind structure.  NOTE THAT THESE LAST ROUTINES
			HAVE NOT BEEN SERIOUSLY RECHECKED, SINCE THEY ARE NOT
			MUCH USED.
	00jan	ksl	Moved ionization_mod_on_the_spot and ionization_mod_on_the_spot_exact
			to "Legacy code" Legacy1.c
	01oct	ksl	Add calls to levels for calculation of level populations.  The
			calls are currently hardwired.
	08aug	ksl	60b	Evidently ksl modified the calls to ion_abundances
				previously, but leaving w as the variable was
				confusing.  Fixed this
	080808	ksl	62  - Removed option 4 for a partial detailed balance
			which seems totally obsolete at this point as
			we have not followed this up.  This was part of
			the cleanup of the ionizaton balance routines
			Also removed option 5, as this was not supported
			later in the program
**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

double sim_numin, sim_numax, sim_meanfreq;	// external variables set up so zbrent can solve for alhpa.



int
ion_abundances (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int ireturn, nion, nelem;
  double zbrent (), sim_alpha_func ();
  double alphamin, alphamax, alphatemp, sim_w_temp;

//      printf("NSH here we are in ion_abundances, I think we are running at mode %i\n",mode);

  if (mode == 0)
    {
/*on-the-spot approximation using existing t_e.   This routine does not attempt 
 * to match heating and cooling in the wind element! */

      if ((ireturn = nebular_concentrations (xplasma, 2)))
	{
	  Error
	    ("ionization_abundances: nebular_concentrations failed to converge\n");
	  Error
	    ("ionization_abundances: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
	     xplasma->j, xplasma->t_e, xplasma->w);
	}
    }
  else if (mode == 1)
    {
      // LTE using t_r  (ksl - checked - 080808

      ireturn = nebular_concentrations (xplasma, 1);
    }
  else if (mode == 2)
    {				//  Hardwired concentrations

      ireturn = fix_concentrations (xplasma, 0);
    }
  else if (mode == 3)
    {
/* On the spot, with one_shot at updating t_e before calculating densities
*/

/* Shift values to old */
      xplasma->dt_e_old = xplasma->dt_e;
      xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
      xplasma->t_e_old = xplasma->t_e;
      xplasma->t_r_old = xplasma->t_r;
      xplasma->lum_rad_old = xplasma->lum_rad;

      ireturn = one_shot (xplasma, mode);

/* Convergence check */
      convergence (xplasma);
    }
  else if (mode == 4)
    {				//LTE with SIM correction this is called from define_wind, sim_alpha and sim_w are set to geo values in define_wind. Not sure this is ever called now, we thought it best to set the values to LTE when in define_wind.
      ireturn = nebular_concentrations (xplasma, 5);
    }
  else if (mode == 5)
    {				// One shot at updating t_e before calculating densities using the SIM correction
/* Shift values to old */

      // This call is after a photon flight, so we *should* have access to j and ave freq, and so we can calculate proper values for W and alpha
      // To avoid problems with solving, we need to find a reasonable range of values within which to search for a solution to eq18. A reasonable guess is that it is around the current value....



      sim_numin = 1.24e15;
      sim_numax = 1.21e19;
      sim_meanfreq = xplasma->ave_freq;

      alphamin = xplasma->sim_alpha - 0.1;
      alphamax = xplasma->sim_alpha + 0.1;
      while (sim_alpha_func (alphamin) * sim_alpha_func (alphamax) > 0.0)
	{
	  alphamin = alphamin - 1.0;
	  alphamax = alphamax + 1.0;
	}


/* We compute temporary values for sim alpha and sim weight. This will allow us to check that they are sensible before reassigning them */

      alphatemp = zbrent (sim_alpha_func, alphamin, alphamax, 0.00001);

/*This next line computes the sim weight using an external function. Note that xplasma->j already contains the volume of the cell and a factor of 4pi, so the volume sent to sim_w is set to 1 and j has a factor of 4PI reapplied to it. This means that the equation still works in balance. It may be better to just implement the factor here, rather than bother with an external call.... */
      sim_w_temp =
	sim_w ((xplasma->j) * 4 * PI, 1, 1, xplasma->sim_alpha, sim_numin,
	       sim_numax);

      printf ("NSH We noew have calculated sim_w_temp=%e\n", sim_w_temp);
      if (sane_check (sim_w_temp))
	{
	  Error
	    ("New sim parameters unreasonable, using existing parameters. Check number of photons in this cell");
	}
      else
	{
	  xplasma->sim_alpha = alphatemp;
	  xplasma->sim_w = sim_w_temp;
	}




      xplasma->sim_ip =
	xplasma->sim_w *
	(((pow (50000 / HEV, xplasma->sim_alpha + 1.0)) -
	  pow (100 / HEV,
	       xplasma->sim_alpha + 1.0)) / (xplasma->sim_alpha + 1.0));
      xplasma->sim_ip *= 16 * PI * PI;
      xplasma->sim_ip /= xplasma->rho * rho2nh;




      xplasma->dt_e_old = xplasma->dt_e;
      xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
      xplasma->t_e_old = xplasma->t_e;
      xplasma->t_r_old = xplasma->t_r;
      xplasma->lum_rad_old = xplasma->lum_rad;

      ireturn = one_shot (xplasma, mode);


      for (nelem = 0; nelem < nelements; nelem++)
	{
	  for (nion = ele[nelem].firstion;
	       nion < ele[nelem].firstion + ele[nelem].nions; nion++)
	    {
//                      printf ("For ion number %i of element %i, Ioniz=%e, Recomb=%e, Density=%e\n",nion-ele[nelem].firstion,ion[nion].z       ,xplasma->ioniz[nion],xplasma->recomb[nion],xplasma->density[nion]);
	    }
	}



/* Convergence check */
      convergence (xplasma);
    }


  else
    {
      Error
	("ion_abundances: Could not calculate abundances for mode %d\n",
	 mode);
      exit (0);
    }

  /* If we want the Auger effect deal with it now. Initially, this is
     put in here, right at the end of the ionization calculation -
     the assumption is that the Auger effect is only for making minor
     ions so that the ionization balance of the other ions is not
     affected in an important way. */

  if (geo.auger_ionization == 1)
    {
      auger_ionization (xplasma);
    }


  return (ireturn);

}


/***********************************************************
              Space Telescope Science Institute

 Synopsis: convergence checks to see whehter a single cell
	is or is not converging
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:

History:
	06may	ksl	57+: Began modifications to reflect new
			structure definition.  Note that it
			is likelely that this entire routine
			will ultimatetely be replaced because
			everything here should only be in the wind
**************************************************************/
int
convergence (xplasma)
     PlasmaPtr xplasma;
{
  int trcheck, techeck, hccheck, whole_check, converging;
  double epsilon;

  trcheck = techeck = hccheck = converging = 0;
  epsilon = 0.05;

  if ((xplasma->converge_t_r =
       fabs (xplasma->t_r_old - xplasma->t_r) / (xplasma->t_r_old +
						 xplasma->t_r)) > epsilon)
    trcheck = 1;
  if ((xplasma->converge_t_e =
       fabs (xplasma->t_e_old - xplasma->t_e) / (xplasma->t_e_old +
						 xplasma->t_e)) > epsilon)
    techeck = 1;
  if ((xplasma->converge_hc =
       fabs (xplasma->heat_tot - xplasma->lum_rad) / (xplasma->heat_tot +
						      xplasma->lum_rad)) >
      epsilon)
    hccheck = 1;

  xplasma->converge_whole = whole_check = trcheck + techeck + hccheck;

  if (xplasma->dt_e_old * xplasma->dt_e < 0
      && fabs (xplasma->dt_e) > fabs (xplasma->dt_e_old))
    converging = 1;
  xplasma->converging = converging;

  if (converging == 1)
    {				// Not converging
      xplasma->gain *= 0.7;
      if (xplasma->gain < 0.1)
	xplasma->gain = 0.1;
    }
  else
    {
      xplasma->gain *= 1.1;
      if (xplasma->gain > 0.8)
	xplasma->gain = 0.8;
    }

  return (whole_check);
}





/***********************************************************
              Space Telescope Science Institute

 Synopsis:
   check_convergence -- Do a global check on how well the wind is converging
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:
Notes--Eventually should absorb some of the calculations in wind_updates

History:
	05apr	ksl	55d -- Trivially eliminated MDIM from this
			in effort to make coordinate system independent.
			Whether something is in the wind or not is
			now totally determeined by the volume.
        06may	ksl	57+ -- Now using plasma structure
**************************************************************/
int
check_convergence ()
{
  int n;
  int nconverge, nconverging, ntot;
  double xconverge, xconverging;

  nconverge = nconverging = ntot = 0;

  for (n = 0; n < NPLASMA; n++)
    {
      ntot++;
      if (plasmamain[n].converge_whole == 0)
	nconverge++;
      if (plasmamain[n].converging == 0)
	nconverging++;

    }

  xconverge = ((double) nconverge) / ntot;
  xconverging = ((double) nconverging) / ntot;
  Log
    ("!!Check_converging: %4d (%.3f) converged and %4d (%.3f) converging of %d cells\n",
     nconverge, xconverge, nconverging, xconverging, ntot);
  Log
    ("Summary  convergence %4d %.3f  %4d  %.3f  %d  #  n_converged fraction_converged  converging fraction_converging total cells\n",
     nconverge, xconverge, nconverging, xconverging, ntot);
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	one_shot calculates new densities of ions in a single element of the wind
	according to equation 11 of Lucy & Mazzali, after having found the
	temperature which matches heating and cooling for the previous
	densities
	
 Arguments:		
	PlasmaPtr xplasma;

Returns:
 
Description:
	
Notes:
	This routine attempts to match heating and cooling in the wind element!
	To do this it calls calc_te.  Based on the returned value of te, the
	routine then calculates the concentrations in the on-the-spot approximation.

	IT SEEMS LIKELY some code could be eliminated by simply having this routine
	call the on-the-spot routine directly.

History:
	98	ksl	Coded as part of python effort
	02jul	ksl	Added mode variable so could try detailed balance
	06may	ksl	57+ -- Switched to use plasma structue

**************************************************************/

int
one_shot (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;

{
  double te_old, te_new, dte;
  double gain;




  gain = xplasma->gain;

  te_old = xplasma->t_e;
  te_new = calc_te (xplasma, 0.7 * te_old, 1.3 * te_old);
  xplasma->t_e = (1 - gain) * te_old + gain * te_new;

  dte = xplasma->dt_e;

  Log ("One_shot: cell %i %13.6f %13.6f mode %i\n", xplasma->nplasma, te_old,
       te_new, mode);


/* Modes in the driving routines are not identical to those in nebular concentrations.
The next lines are an attempt to mediate this problem.  It might be better internally
at least to define a flag for using one shot, and have the modes take on the
meaning in nebular concentrations.
*/

  if (mode == 3)
    mode = 2;
  else if (mode <= 1 || mode >= 6)	/* modification to cope with mode 5 - SIM */
    {

      Error ("one_shot: Sorry, Charlie, don't know how to process mode %d\n",
	     mode);
      exit (0);
    }

  if (xplasma->t_r > 10.)
    {				/* Then modify to an on the spot approx */
      if (nebular_concentrations (xplasma, mode))
	{
	  Error
	    ("ionization_on_the_spot: nebular_concentrations failed to converge\n");
	  Error
	    ("ionization_on_the_spot: j %8.2e t_e %8.2e t_r %8.2e w %8.2e\n",
	     xplasma->j, xplasma->t_e, xplasma->w);
	}
      if (xplasma->ne < 0 || VERY_BIG < xplasma->ne)
	{
	  Error ("ionization_on_the_spot: ne = %8.2e out of range\n",
		 xplasma->ne);
	}
    }
  else
    {
      Error ("ionization_on_the_spot: t_r exceptionally small %g\n",
	     xplasma->t_r);
      mytrap ();
      exit (0);
    }


  return (0);
}


/* calc_te determines and returns the electron temperature in the wind such that the energy emitted
   by the wind is equal to energy emitted.

   Description:
   calc_te does not modify any abundances.  It simply takes the current value of the heating in the
   cell and attempts to find the value of the electron temperature which will result in cooling which
   matches the heating.

   This approach is not entirely self consistent because if te is actually different then the
   abundances will be different and the heating will change as well.

   This routine is a kluge because it does not really deal with what happens if the cooling curve 
   has maxima and minima.

   xxxplasma is just a way to tranmit information to zero_emit

   History:

   98dec        ksl     Updated calls so that tmin and tmax were communicated externally,
			rather than hardwired
   01dec	ksl	Reversed 98dec decision as a result of massive changes for python38
   01dec	ksl	Assured that ww->t_e is updated here
   01dec	ksl	Added capability to modify the desired goal of calc_te from the full
			heating to something intermediate between the current value and the
			ultimate goal
   01dec	ksl	Rewrote to assure that boundaries will be bracketed properly and if
			not calc_te will handle
   04June       SS      Modified so that changes in the heating rate due to changes in the
                        temperature are included for macro atoms.
	06may	ksl	Modified for plasma structue
 */

PlasmaPtr xxxplasma;

double
calc_te (xplasma, tmin, tmax)
     PlasmaPtr xplasma;
     double tmin, tmax;
{
  double heat_tot;
  double z1, z2;
  int macro_pops ();
  xxxplasma = xplasma;
  heat_tot = xplasma->heat_tot;

  xplasma->t_e = tmin;
  z1 = zero_emit (tmin);
  xplasma->t_e = tmax;
  z2 = zero_emit (tmax);

  if ((z1 * z2 < 0.0))
    {				// Then the interval is bracketed 
      xplasma->t_e = zbrent (zero_emit, tmin, tmax, 50.);
    }
  else if (fabs (z1) < fabs (z2))
    {
      xplasma->t_e = tmin;
    }
  else
    xplasma->t_e = tmax;

  /* With the new temperature in place for the cell, get the correct value of heat_tot.
     SS June  04 */

  xplasma->heat_tot -= xplasma->heat_lines_macro;
  xplasma->heat_lines -= xplasma->heat_lines_macro;
  xplasma->heat_lines_macro = macro_bb_heating (xplasma, xplasma->t_e);
  xplasma->heat_tot += xplasma->heat_lines_macro;
  xplasma->heat_lines += xplasma->heat_lines_macro;

  xplasma->heat_tot -= xplasma->heat_photo_macro;
  xplasma->heat_photo -= xplasma->heat_photo_macro;
  xplasma->heat_photo_macro = macro_bf_heating (xplasma, xplasma->t_e);
  xplasma->heat_tot += xplasma->heat_photo_macro;
  xplasma->heat_photo += xplasma->heat_photo_macro;



  return (xplasma->t_e);

}


//double
//calc_te (xplasma, tmin, tmax)
//     PlasmaPtr xplasma;
//     double tmin, tmax;
//{
//  double heat_tot;
//  double z1, z2;
//  int macro_pops ();
//  xxxplasma = xplasma;
//  heat_tot = xplasma->heat_tot;
//  xplasma->t_e = tmin = 1;
//  z1 = zero_emit (tmin);
//  xplasma->t_e = tmax = 1e7;
//  z2 = zero_emit (tmax);
//
//      printf("NSH In calc_te  z1%e, z2=%e, tmin=%f, tmax=%f\n",z1,z2,tmin,tmax); 
//  if ((z1 * z2 < 0.0))
//    {                         // Then the interval is bracketed 
//      xplasma->t_e = zbrent (zero_emit, tmin, tmax, 50.);
//      printf("NSH calc_te Bracketed, t_e set to %f\n",xplasma->t_e);
//    }
//  else if (fabs (z1) < fabs (z2))
//    {
//      xplasma->t_e = tmin;
//      printf("NSH calc_te t_e set to tmin (%f) coz %e<%e\n",xplasma->t_e,fabs (z1) , //fabs (z2));
//    }
//  else
//    {
//   xplasma->t_e = tmax;
//      printf("NSH calc_te t_e set to tmax (%f)\n",xplasma->t_e);
//      }
  /* With the new temperature in place for the cell, get the correct value of heat_tot.
     SS June  04 */

//  xplasma->heat_tot -= xplasma->heat_lines_macro;
//  xplasma->heat_lines -= xplasma->heat_lines_macro;
//  xplasma->heat_lines_macro = macro_bb_heating (xplasma, xplasma->t_e);
//  xplasma->heat_tot += xplasma->heat_lines_macro;
//  xplasma->heat_lines += xplasma->heat_lines_macro;

//  xplasma->heat_tot -= xplasma->heat_photo_macro;
//  xplasma->heat_photo -= xplasma->heat_photo_macro;
//  xplasma->heat_photo_macro = macro_bf_heating (xplasma, xplasma->t_e);
//  xplasma->heat_tot += xplasma->heat_photo_macro;
//  xplasma->heat_photo += xplasma->heat_photo_macro;



//  return (xplasma->t_e);

//





/* This is just a function which has a zero when total energy loss is equal to total energy gain */

double
zero_emit (t)
     double t;
{
  double difference;
  double total_emission ();
  int macro_pops ();
  double macro_bb_heating (), macro_bf_heating ();


  /* NOTE - IMPORTANT

     SS May 04: I'm removing adiabatic cooling for now. I'll need to be put
     back in once the heating/cooling tests are all sorted out. 

     SS June 04: Adiabatic cooling is still not switched on. But it is now 
     switched off in emission.c rather than here.

     SS July 04: Adiabatic cooling is now switched back on. These comments can
     all be deleted once this is tested (soon).
   */

  /* This block is not needed now - SS July 04
     if (xxxplasma->lum_adiabatic > 0)
     {
     Error("zero_emit: adiabatic cooling is switched on.\n");
     }
   */

  /*Original method */
  xxxplasma->t_e = t;


  /* Correct heat_tot for the change in temperature. SS June 04. */
  //macro_pops (xxxplasma, xxxplasma->ne);
  xxxplasma->heat_tot -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines -= xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines_macro = macro_bb_heating (xxxplasma, t);
  xxxplasma->heat_tot += xxxplasma->heat_lines_macro;
  xxxplasma->heat_lines += xxxplasma->heat_lines_macro;

  xxxplasma->heat_tot -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo -= xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo_macro = macro_bf_heating (xxxplasma, t);
  xxxplasma->heat_tot += xxxplasma->heat_photo_macro;
  xxxplasma->heat_photo += xxxplasma->heat_photo_macro;

  //  difference = (xxxplasma->heat_tot - total_emission (xxxplasma, 0., VERY_BIG));
  difference =
    xxxplasma->heat_tot - xxxplasma->lum_adiabatic -
    total_emission (&wmain[xxxplasma->nwind], 0., VERY_BIG);
//     printf("NSH Zero emit returns %e for t_e= %f\n",difference,t);
  return (difference);
}

double
sim_alpha_func (alpha)
     double alpha;
{
  double answer;
  answer =
    ((alpha + 1.) / (alpha + 2.)) * ((pow (sim_numax, (alpha + 2.)) -
				      pow (sim_numin,
					   (alpha + 2.))) / (pow (sim_numax,
								  (alpha +
								   1.)) -
							     pow (sim_numin,
								  (alpha +
								   1.))));

  answer = answer - sim_meanfreq;
//      printf("NSH alpha=%f,f1=%e,f2=%e,meanfreq=%e,ans=%e\n",alpha,sim_numin,sim_numax,sim_meanfreq,answer);
  return (answer);
}
