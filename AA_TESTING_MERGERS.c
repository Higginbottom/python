
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************

This is just a test C function to see how mergers work.

Imagine we have 2 workers, JM and OtherPerson

1- JM initially wrote code (copied from photon_ge.c more or less)

**********************************************************/

double main(r, tstar, freqmin, freqmax, ioniz_or_final, f)
{
  double lumstar;
  double log_g;
  double emit, emittance_bb (), emittance_continuum ();
  int spectype;

  log_g = geo.gstar = log10 (G * geo.mstar / (geo.rstar * geo.rstar));
  lumstar = 4 * PI * STEFAN_BOLTZMANN  * r * tstar * tstar * tstar * tstar; /* Other person edited this */

  if (ioniz_or_final == 1)
    spectype = geo.star_spectype;	/* type for final spectrum */
  else
    spectype = geo.star_ion_spectype;	/*type for ionization calculation */

  if (spectype >= 0)
    {
      emit = emittance_continuum (spectype, freqmin, freqmax, tstar, log_g);
    }
  else
    {

      emit = emittance_bb (freqmin, freqmax, tstar);
      /* Other person edited this */
    }

  *f = emit;			/* Calculate the surface flux between freqmin and freqmax */
  *f *= (4. * PI * r * r);

  geo.lum_star = lumstar;
  return (lumstar);
}

/* ALL DONE */
