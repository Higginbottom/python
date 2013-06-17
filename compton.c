


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  These are the routines to handle compton scattering.

  Description:	There are two routines in here:
		kappa_compton (xplasma,freq) this calculates the opacity in the cell xplasma due to compton cooling

  Arguments: 		


  Returns:

  Notes:

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  kappa_compton computes the opacity in the cell due to compton heating.

  Description:	

  Arguments:   xplasma - pointer to the current plasma cell
		freq - the frequency of the current photon packet being tracked

  Returns:   kappa - the compton opacity for the cell.

  Notes:   This implements the equation kappa=(sigmaT*ne*freq)/(me*c^2)
            Initally it is called by radiation

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.
feb 2013 - nsh - approximate KN cross section replaced by correct value

 ************************************************************************/


double
kappa_comp (xplasma, freq)
     PlasmaPtr xplasma;		// Pointer to current plasma cell
     double freq;		// Frequency of the current photon being tracked
{
  double x;			// The opacity of the cell by the time we return it.
  double sigma;			/*The cross section, thompson, or KN if hnu/mec2 > 0.01 */

/*	alpha=1/(1+freq*HRYD*(1.1792e-4+(7.084e-10*freq*HRYD))); NSH 130214 This is the approximate way of doing it.*/

  sigma = klein_nishina (freq);	//NSH 130214 - full KN formula

  x = (sigma * H) / (MELEC * C * C);	//Calculate the constant
  x *= xplasma->ne * freq;	//Multiply by cell electron density and frequency of the packet.
  return (x);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  kappa_ind_compton computes the opacity in the cell due to induced compton heating.

  Description:	

  Arguments:   xplasma - pointer to the current plasma cell
		freq - the frequency of the current photon packet being tracked

  Returns:   kappa - the inducd compton opacity for the cell.

  Notes:   This implements the induced compton term in equation 6.5 in cloudy 10.3

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.
feb 2013 - nsh - approximate KN cross section replaced by correct value

 ************************************************************************/


double
kappa_ind_comp (xplasma, freq, ds, w)
     PlasmaPtr xplasma;		// Pointer to current plasma cell
     double freq;		// Frequency of the current photon being tracked
     double w;			// The weight of the photon packet
     double ds;			//The distance the photon travels
{
  double x;			// The opacity of the cell by the time we return it.
  double sigma;			/*The cross section, thompson, or KN if hnu/mec2 > 0.01 */
  double J, expo;		//The estimated intensity in the cell
  int i;
//      J=(4*PI*w*ds)/(C*xplasma->vol); //Calcuate the intensity NSH This works for a thin shell... Why? Dont know.
  J = 0.0;			/* NSH 130605 to remove o3 compile error */
  if (geo.ioniz_mode == 5 || geo.ioniz_mode == 7)	/*If we are using power law ionization, use PL estimators */
    {
      for (i = 0; i < geo.nxfreq; i++)
	{
	  if (geo.xfreq[i] < freq && freq <= geo.xfreq[i + 1])	//We have found the correct model band
	    {
	      if (xplasma->spec_mod_type[i] < 0)	//Only bother if we have a model in this band
		{
		  J = 0.0;	//THere is no modelin this band, so the best we can do is assume zero J
		}
	      else if (xplasma->spec_mod_type[i] == SPEC_MOD_PL)	//Power law model
		{
		  J = xplasma->pl_w[i] * pow (freq, xplasma->pl_alpha[i]);
		}
	      else if (xplasma->spec_mod_type[i] == SPEC_MOD_EXP)	//Exponential model
		{
		  J =
		    xplasma->exp_w[i] * exp ((-1 * H * freq) /
					     (BOLTZMANN *
					      xplasma->exp_temp[i]));
		}
	      else
		{
		  Error
		    ("kappa_ind_comp - unknown spectral model (%i) in band %i\n",
		     xplasma->spec_mod_type[i], i);
		  J = 0.0;	//Something has gone wrong
		}
	    }
	}
    }
  else				/*Else, use BB estimator of J */
    {
      expo = (H * freq) / (BOLTZMANN * xplasma->t_r);
      J = (2 * H * freq * freq * freq) / (C * C);
      J *= 1 / (exp (expo) - 1);
      J *= xplasma->w;
    }




//      printf("We think ds is %e vol is %e area is %e\n",ds,xplasma->vol,pow((xplasma->vol/(4*PI*ds)),0.5));
//      alpha1=THOMPSON/(1+freq*HRYD*(1.1792e-4+(7.084e-10*freq*HRYD))); //KN cross section NSH 130214 - approximate

  sigma = klein_nishina (freq);	//NSH 130214 - full KN formula


  x = (xplasma->ne) / (MELEC);
  x *= sigma * J;		// NSH 130214 factor of THOMPSON removed, since alpha is now the actual compton cross section
  x *= 1 / (2 * freq * freq);
  if (sane_check (x))
    {
      Error
	("kappa_ind_comp:sane_check - undefined value for Kappa_ind_comp - setting to zero\n");
      return (0.0);
    }
  return (x);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  total_comp computes the cooling in the cell due to compton cooling.

  Description:	

  Arguments:   one - pointer to the current wind cell
		t_e - electron temperature of the current cell

  Returns:   the compton luminosity for the cell.

  Notes:   This implements the equation C=16piVsigmatJ(kTe/me)

  History:
2011	nsh	Coded as part of the effort to include compton scattering in August 2011 at Southampton.

 ************************************************************************/

double
total_comp (one, t_e)
     WindPtr one;		// Pointer to the current wind cell - we need the cell volume, this is not in the plasma structure
     double t_e;		//Current electron temperature of the cell
{
  double x;			//The returned variable
  int nplasma;			//The cell number in the plasma array
  PlasmaPtr xplasma;		//pointer to the relevant cell in the plasma structure


  nplasma = one->nplasma;	//Get the correct plasma cell related to this wind cell
  xplasma = &plasmamain[nplasma];	//copy the plasma structure for that cell to local variable

  x = 16. * PI * THOMPSON * BOLTZMANN / (MELEC * C * C);	//Keep all the constants together
  x *= xplasma->ne * one->vol * xplasma->j * t_e;	//multiply by the volume (from wind) and j (from plasma) and t_e


  return (x);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  klein_nishina computes the KN cross section for a photon of frequency nu.

  Description:	

  Arguments:   nu - photon frequency


  Returns:   the KN cross section.

  Notes:   This implements equation 7.5 in Rybicki and Lightman

  History:
2013	nsh	Coded

 ************************************************************************/

double
klein_nishina (nu)
     double nu;			//The frequency of the photon packet
{
  double x;			//h nu / kt
  double x1, x2, x3, x4;	//variables to store intermediate results.
  double kn;			// the final cross section

  kn = THOMPSON;		/* NSH 130605 to remove o3 compile error */
  x1 = x2 = x3 = x4 = 0.0;	/* NSH 130605 to remove o3 compile error */
  x = (H * nu) / (MELEC * C * C);
  if (x > 0.0001)
    {
      x1 = 1. + x;
      x2 = 1. + (2. * x);
      x3 = log (x2);
      x4 = ((2. * x * x1) / x2) - x3;
      x4 *= x1 / (x * x * x);
      x4 = x4 + x3 / (2. * x);
      x4 = x4 - (1 + 3. * x) / (x2 * x2);
      kn *= 0.75 * THOMPSON * x4;
    }

  return (kn);
}
