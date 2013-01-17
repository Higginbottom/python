
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"





/***********************************************************
                Space Telescope Science Institute

Synopsis:
Python uses a form of stratified sampling in an attempt to assure
that there are photon bundles at (high, generally) frequencies

This is the routine that initializes the bands.  It is hardwired
at present for H and He, and so may not be ideal for O VI, C IV,
or NV.

   
Arguments:		
	mode	0	Use temperature to define a single band
		1	Use f1 and f2 to define a single band
		2	Use t,f1 and f2, and hardwired bands
			to define multiple bands
Returns:
 
 
Description:	

		
Notes:

0812 - In python, bands_init is only called with mode 2. Other diagnostic
routines, e.g balance make use of the other possibilities.  But we 
need to address the problem that with attempts to broaden the
utility of the python, we need more flexibility in to set the bands.

History:
	01dec	ksl	python40
	02jul	ksl	Adapted from an similar routine used to set
			up bands in photon_gen.  This may eventually
			replace that routine
	0810	ksl	67 - Cleaned up slightly during relook
			at balance
	0812	ksl	67c - Changed names to bands_init to make
			it clear other routines created for diagnostic
			purposes had to be updated, as modified 
			to handle yso and cv type models in same
			code
**************************************************************/


/* Actual structures are in python.h.  Here for reference only.

#define NBANDS 10
struct xbands
{
	double f1[NBANDS],f2[NBANDS];
	double min_fraction[NBANDS];
	double nat_fraction[NBANDS];          // The fraction of the accepted luminosity in this band
	double used_fraction[NBANDS];
	double flux[NBANDS];                  //The "luminosity" within a band
	double weight[NBANDS];
	int nphot[NBANDS];
	int nbands;           // Actual number of bands in use
}
xband;


*/



int
bands_init (t, f1, f2, imode, band)
     double t;			// A temperature which can be used to set absolute limits on the bands
     double f1, f2;		// frequency limits that can overide any other limits
     int imode;			// A switch used for determining how the bands are to be populated
     struct xbands *band;

{
  int mode;
  int nband;
  double xx;



  if (imode == -1)
    {
      mode = 2;
      rdint
	("Photon.sampling.approach(0=T,1=(f1,f2),2=cv,3=yso,4=user_defined),",
	 &mode);
    }
  else
    {
      mode = imode;
    }

  if (mode == 0)
    {
      band->nbands = 1;
      band->f1[0] = BOLTZMANN * t / H * 0.05;
      band->f2[0] = BOLTZMANN * t / H * 20.;
      band->min_fraction[0] = 1.0;
    }
  else if (mode == 1)
    {
      band->nbands = 1;
      band->f1[0] = f1;
      band->f2[0] = f2;
      band->min_fraction[0] = 1.0;
    }
  else if (mode == 2)  /* Traditional cv setup */
    {

      band->nbands = 4;
      band->f1[0] = f1;
      band->f2[0] = band->f1[1] = 13.599 / HEV;
      band->f2[1] = band->f1[2] = 24.588 / HEV;
      band->f2[2] = band->f1[3] = 54.418 / HEV;
      band->f2[3] = f2;
      band->min_fraction[0] = 0;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.1;
      band->min_fraction[3] = 0.1;
      if (f1 > band->f2[0])
	{
	  Error ("bands_init: f1 (%e) > 13.599/HEV)\n", f1);
	  exit (0);
	}
      if (f2 < band->f2[2])
	{
	  Error ("bands_init: f2 (%e) < 54.418/HEV)\n", f2);
	  exit (0);
	}

    }
  else if (mode == 3)  /* YSO setup */
    {
      band->nbands = 4;
      band->f1[0] = f1;
      band->f2[0] = band->f1[1] = 1.511 / HEV;
      band->f2[1] = band->f1[2] = 3.3998 / HEV;
      band->f2[2] = band->f1[3] = 6.0000 / HEV;
      band->f2[3] = f2;
      band->min_fraction[0] = 0;
      band->min_fraction[1] = 0.1;
      band->min_fraction[2] = 0.2;
      band->min_fraction[3] = 0.2;
      if (f1 > band->f2[0])
	{
	  Error ("bands_init: f1 (%e) > 13.599/HEV)\n", f1);
	  exit (0);
	}
      if (f2 < band->f2[2])
	{
	  Error ("bands_init: f2 (%e) < 54.418/HEV)\n", f2);
	  exit (0);
	}

    }
  else if (mode == 4)
    {
      rdint ("Num.of.frequency.bands", &band->nbands);
      Log ("Lowest photon energy is ev (freq) is %f (%.2e)\n", f1 * HEV, f1);
      Log ("Highest photon energy is ev (freq) is %f (%.2e)\n", f2 * HEV, f2);
      printf
	("Enter band boundaries in increasing eV, and assure they are between lowest and highest energy\n");

      band->f1[0] = f1;

      for (nband = 0; nband < band->nbands - 1; nband++)
	{
	  rddoub ("Band.boundary(eV)", &xx);
	  band->f2[nband] = band->f1[nband + 1] = xx / HEV;

	}
      band->f2[nband - 1] = f2;

      printf
	("Enter mimimum fraction of photons in each band.  The total must be < or = to 1\n");

      for (nband = 0; nband < band->nbands; nband++)
	{
	  rddoub ("Band.minimum_fraction)", &band->min_fraction[nband]);
	}

    }

  else
    {
      Error ("Init bands: Unknown mode %d\n", mode);
      mytrap ();
    }


  return (0);
}
