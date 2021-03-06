#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_stellar_wind_params gets input data which is necessary for a Castor and Larmors'
	(1979) description of an O star wind, e.g.
   		v(r)=V_o + (V_infinity-V_o) (1-R/r)**beta

Arguments:		

Returns:
 
Description:	
	The parameters, zdom[ndom].cl...,  obtained here are only used in the routines in stellar_winds.c
	which calculate the velocity and density of the wind during the initialization process.
	Other portions of the structure, zdom[ndom] defined here are more general purpose.

Notes:
	Although it seems natural to assume that the wind starts at the photosphere, this
	is not required (and in Mauche and Raymond's spherical model for a CV wind they set
	the inner radius of the wind at several WD radii. (04jun--ksl-It does not look to me
	as if the code as written really implements the possibility of a wind that starts
	at a radius larger than the star.  This should be checked if one really wants to
	do this.)


History:
 	98dec	ksl	Coded and debugged as part of major change in IO structure required when
 			adding a spherical wind
	99dec	ksl	Corrected coding error, which resulted in geo.cl_rmin not being set
			to geo.wind_rmin.  Note that it is not clear that geo.wind_rmind,
			and goe.cl_rmin are needed separately (and someday this should be corrected.)
	04jun	ksl	Moved from python.c to this file for consistency, and moved reading
			of wind.mdot to this file.  Also added a new variable to geo, namely
			geo.stellar_wind_modt so that could distinguish from stellar and
			disk_wind...which is needed for the yso models.  This implies that
			old .pf files for stellar winds will have to be modified, but there
			are few of these.  Deleted lines which eliminated the disk for 
			stellar models.
	15aug	ksl	Change to accomodate multiple domains
**************************************************************/


int
get_stellar_wind_params (ndom)
     int ndom;
{
  Log ("Creating a wind model for a Star\n");




  zdom[ndom].stellar_wind_mdot = 1.e-6;
  zdom[ndom].rmin = geo.rstar;
  zdom[ndom].cl_beta = 1.0;
  zdom[ndom].cl_rmin = 2.8e9;
  zdom[ndom].cl_v_zero = 200e5;
  zdom[ndom].cl_v_infinity = 3000e5;

  rddoub ("stellar_wind_mdot(msol/yr)", &zdom[ndom].stellar_wind_mdot);
  zdom[ndom].stellar_wind_mdot *= MSOL / YR;

  rddoub ("stellar.wind.radmin(cm)", &zdom[ndom].rmin); /*Radius where wind begins */
  if (zdom[ndom].rmin < geo.rstar)
  {
    Error ("get_stellar_wind_params: It is unreasonable to have the wind start inside the star!\n");
    Log ("Setting geo.rmin to geo.rstar\n");
    zdom[ndom].rmin = geo.rstar;
  }
  zdom[ndom].cl_rmin = zdom[ndom].rmin;
  rddoub ("stellar.wind_vbase(cm)", &zdom[ndom].cl_v_zero);     /* Velocity at base of the wind */
  rddoub ("stellar.wind.v_infinity(cm)", &zdom[ndom].cl_v_infinity);    /* Final speed of wind in units of escape velocity */

  rddoub ("stellar.wind.acceleration_exponent", &zdom[ndom].cl_beta);   /* Accleration scale exponent */

/* Assign the generic parameters for the wind the generic parameters of the wind */
  geo.rmin = zdom[ndom].rmin;   //71 ksl - Not modified this so that we did not waste cells
  zdom[ndom].rmax = geo.rmax;
  zdom[ndom].wind_thetamin = 0.0;
  zdom[ndom].wind_thetamax = 90. / RADIAN;

/* define the the variables that determine the gridding */
  zdom[ndom].wind_rho_min = 0;
  zdom[ndom].wind_rho_max = geo.rmax;
  zdom[ndom].zmax = geo.rmax;


  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = 0.3 * geo.rstar;
    zdom[ndom].zlog_scale = 0.3 * geo.rstar;
  }

  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double stellar_velocity(ndom, x,v) calulates the v the wind at a position 
	x (where both x and v in cartesian coordinates)
Arguments:		
	ndom 			The domain number
	double x[]		the postion where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	
	The model is that of Castor and Lamars 1979  where

	v(r)=V_o + (V_infinity-V_o) (1-R/r)**beta
	
	The values of the individiual constants should all be part of the structure zdom[ndom].

	V_o:  			zdom[ndom].cl_v_zero;		velocity at base of wind 
	V_infinity:		zdom[ndom].cl_v_infinity;	the velocity at infinity
	R			zdom[ndom].cl_rmin	       	the inner radius of the wind
	beta			zdom[ndom].cl_beta		power law exponent for the velocity law

		
Notes:
	v is set to V_o inside of rstar or w_rscale, even though the should not really exist there!

	Since this is a radial wind, it is not affected by the disk at all. In particular there
	are no modifications to this routine if the disk is vertically extended.

History:
 	98dec	ksl	Coded as part of effort to add a stellar wind option to python
	04jun	ksl	Caused routine to return 0 velocity at origin to avoid clutter of irrelevant
			error messages.
 
**************************************************************/

double
stellar_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, speed, zzz;
  double length ();

  if ((r = length (x)) == 0.0)
  {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    return (0.0);
  }


  if (r <= geo.rstar || r <= zdom[ndom].cl_rmin)
    speed = zdom[ndom].cl_v_zero;
  else
  {
    zzz = pow (1. - zdom[ndom].cl_rmin / r, zdom[ndom].cl_beta);
    speed = zdom[ndom].cl_v_zero + (zdom[ndom].cl_v_infinity - zdom[ndom].cl_v_zero) * zzz;
  }
  v[0] = speed * x[0] / r;
  v[1] = speed * x[1] / r;
  v[2] = speed * x[2] / r;

  return (speed);

}

/***********************************************************
		Space Telescope Science Institute

 Synopsis:
	double stellar_rho(x) calculates the density of an stellar_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the density
Returns:
	The density at x is returned in gram/cm**3
	
Description:

	rho=mdot/(4PI*r*r*v);	
		
Notes:

History:
 	98dec	ksl	Modified from sv.c portions of the code.
 
**************************************************************/

double
stellar_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rho, v[3];

  r = length (x);
  if (r < zdom[ndom].cl_rmin)
  {
    rho = zdom[ndom].stellar_wind_mdot / (4. * PI * r * r * zdom[ndom].cl_v_zero);
  }
  else
  {
    rho = zdom[ndom].stellar_wind_mdot / (4. * PI * r * r * stellar_velocity (ndom, x, v));
  }

  return (rho);
}
