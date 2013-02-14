
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************

This is just a test C function to see how mergers work.

Imagine we have 2 workers, JM and OtherPerson

1- JM initially wrote code (copied from photon_ge.c more or less)
2- JM made some random changes

**********************************************************/

double main(r, tstar, freqmin, freqmax, ioniz_or_final, f)
{
  double lumstar;
  double log_g;
  double emit, emittance_bb (), emittance_continuum ();
  int spectype;

  log_g = geo.gstar = log10 (G * geo.mstar / (geo.rstar * geo.rstar *2000.0)); /* JM Edited this */
  lumstar = 8 * PI * STEFAN_BOLTZMANN * r * r * tstar * tstar * tstar; /* JM Edited this */

  if (ioniz_or_final == 1)
    spectype = geo.star_spectype;	/* type for final spectrum */
  else
    spectype = geo.star_ion_spectype;	/*type for ionization calculation */

  if (spectype >= 0)
    {
      /* JM Edited this */
      emit = emittance_continuum (spectype, freqmin, freqmax, tstar, log_g) *100.0; /* JM Edited this */
    }
  else
    {
      emit = emittance_bb (freqmin, freqmax, tstar);
    }

  *f = emit;			/* Calculate the surface flux between freqmin and freqmax */
  *f *= (4. * PI * r * r);

  geo.lum_star = lumstar;
  return (lumstar);
}

/* ALL DONE */
