#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic.h"
#include "python.h"

#include "models.h"
//




/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

double one_continuum(spectype,t,g,freqmin,freqmax) gets a photon frequency 
from a continuum grid 

Arguments:

	spectype				An index to the continum grid which one wants
						to use for obtaining a photon
	double t,g,				Two parameters defining the model for which 
						a photon is requried, often temperature and 
						gravity, though this is not required
	freqmin,freqmax;			minimum and maximum frequency of interest


Returns:

 
Description:	
	

Notes:
	In versions of python prior to python_52, different routines
	were used to read kurucz models and hubeny models.  This was
	historical, since the routines ksl had to read the kurucz models
	read binary files (from a time when compputers were very slow).

	The new routines are much more general, and are based on routines
	ksl had written for his fitting program kslfit.  The underlying
	routines can interpolate between grids of more than two dimensions
	although here we assume that all of the continuum grids are two
	dimensional.  

History:
	04aug	ksl	created from hub.c usied in all versions
			of python prior to python_52.  
 
**************************************************************/



double old_t, old_g, old_freqmin, old_freqmax;
double jump[] = { 913.8 };

double
one_continuum (spectype, t, g, freqmin, freqmax)
     int spectype;
     double t, g, freqmin, freqmax;
{
  //OLD not used in routine double par[2];              // For python we assume only two parameter models
  double lambdamin, lambdamax;
  double w_local[NCDF],f_local[NCDF];
  double f,y;
  double cdf_get_rand ();
  // int model (), nwav;
  int model ();
  int n,nwave;

  if (old_t != t || old_g != g || old_freqmin != freqmin || old_freqmax != freqmax)
  {                             /* Then we must initialize */
//	  printf ("Initialising coz %e ne %e or %e ne %e or %e ne %e or %e ne %e\n",old_t,t,old_g,g,old_freqmin,freqmin,old_freqmax,freqmax);
	  
      lambdamin = C * 1e8 / freqmax;
      lambdamax = C * 1e8 / freqmin;
	  nwave=0;
	  if (comp[spectype].xmod.w[0] < lambdamin && lambdamin < comp[spectype].xmod.w[comp[spectype].nwaves-1])
	  {
		  w_local[nwave]=lambdamin;
		  linterp(lambdamin, comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
		  f_local[nwave]=y;
		  nwave++;
	  }
	  
	  for (n=0;n<comp[spectype].nwaves;n++)
	  {
//		  printf ("%e %e ",comp[spectype].xmod.w[n],comp[spectype].xmod.f[n]);
		  if (comp[spectype].xmod.w[n] > lambdamin && comp[spectype].xmod.w[n] <= lambdamax)
		  {
//			  printf (" inserted %e %e\n",lambdamin,lambdamax);
			  w_local[nwave]=comp[spectype].xmod.w[n];
			  f_local[nwave]=comp[spectype].xmod.f[n];
			  nwave++;
		  }
//		  else
//		  {
//		  printf ("\n ");
//	  		}
		  
	  }
	  
	  if (comp[spectype].xmod.w[0] < lambdamax && lambdamax < comp[spectype].xmod.w[comp[spectype].nwaves-1])
	  {
		  w_local[nwave]=lambdamax;
		  linterp(lambdamax, comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
		  f_local[nwave]=y;
		  nwave++;
	  }
	  
	  
	  //There are two pathological cases to deal with, when we only have one non zero point, we need to make an extra point just up/dpown from the penultimate/zecond point so we can make a sensible CDF.
	  
	  
	  if (f_local[nwave-2]==0.0) //We have a zero just inside the end
	  {
		  nwave++;
		  w_local[nwave-1]=w_local[nwave-2];
		  f_local[nwave-1]=f_local[nwave-2];
		  w_local[nwave-2]=w_local[nwave-3]/(1.-DELTA_V/(2.*C));
		  linterp(w_local[nwave-2], comp[spectype].xmod.w, comp[spectype].xmod.f, comp[spectype].nwaves, &y, 0);
		  f_local[nwave-2]=y;
	  }
	  
	  
	  
	  
	  
	  
//	  printf ("BLAH min %e max %e\n",lambdamin,lambdamax);
//	  printf ("BLAH wmin %e wmax %e\n",comp[spectype].xmod.w[0],comp[spectype].xmod.w[comp[spectype].nwaves-1]);
	  
//	  for (n=0;n<nwave;n++)
//		  printf ("BLAH %e %e\n",w_local[n],f_local[n]);
	  
    //OLD not used in routine par[0] = t;
    //OLD not used in routine par[1] = g;
    //OLD nwav not used here ksl.  nwav = model (spectype, par);
    /*  Get_model returns wavelengths in Ang and flux in ergs/cm**2/Ang */

    if (cdf_gen_from_array
        (&comp[spectype].xcdf, w_local, f_local, nwave, lambdamin, lambdamax) != 0)
    {
      Error ("In one_continuum after return from cdf_gen_from_array\n");
    }
    old_t = t;
    old_g = g;
    old_freqmin = freqmin;
    old_freqmax = freqmax;
  }

  f = (C * 1.e8 / cdf_get_rand (&comp[spectype].xcdf));
  if (f > freqmax)
  {
    Error ("one_continuum: f too large %e\n");
    f = freqmax;
  }
  if (f < freqmin)
  {
    Error ("one_continuum: f too small %e\n");
    f = freqmin;
  }
  return (f);
}



double
emittance_continuum (spectype, freqmin, freqmax, t, g)
     int spectype;
     double freqmin, freqmax, t, g;
{
  int nwav, n;
  double w, x, lambdamin, lambdamax;
  double dlambda;
  double par[2];
  int model ();

  lambdamin = C / (freqmax * ANGSTROM);
  lambdamax = C / (freqmin * ANGSTROM);
  par[0] = t;
  par[1] = g;
  model (spectype, par);
  nwav=comp[spectype].nwaves;

  if (lambdamax > comp[spectype].xmod.w[nwav- 1] || lambdamin < comp[spectype].xmod.w[0])
  {

    Error ("emittance_continum: Requested wavelengths extend beyond models wavelengths for list %s\n", comp[spectype].name);
    Error ("lambda %f %f  model %f %f\n", lambdamin, lambdamax, comp[spectype].xmod.w[0], comp[spectype].xmod.w[nwav - 1]);

  }
  x = 0;
  for (n = 0; n < nwav; n++)
  {
    w = comp[spectype].xmod.w[n];
    if (n == 0)
    {
      dlambda = comp[spectype].xmod.w[1] - comp[spectype].xmod.w[0];
    }
    else if (n == nwav - 1)
    {
      dlambda = comp[spectype].xmod.w[n] - comp[spectype].xmod.w[n - 1];
    }
    else
    {
      dlambda = 0.5 * (comp[spectype].xmod.w[n + 1] - comp[spectype].xmod.w[n - 1]);
    }
    if (lambdamin < w && w < lambdamax)
    {
      x += comp[spectype].xmod.f[n] * dlambda;
    }
  }
  x *= 4. * PI;
  return (x);
}
