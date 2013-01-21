

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	Compute the abundances baset on the power law approximation
	use by Stuart in one of his early AGN papers.
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:

History:
	111126	ksl	This was originally coded by nsh, but it
			was inserted directly into what was intended
			to be the steering routine for all ionization
			calculations.  I have moved it to its own
			routine in python_71
	111227	ksl	Eliminated n as an external variable
			because it really is not one
	111229	ksl	Small modifications made to reflect moving
			nxfreq and xfreq into the geo structure
        120802  nsh	Recoded so all this does is compute the estimators
			this is needed for other PL modes.
        120817  nsh     Changed many things to allow for an exponential model
			to be produced. Changed name of code to spec_estimators
			and changed file it is stored in to spectral_estimators.c

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/* external variables set up so zbrent can solve for alpha using sim_appha as integrand. NSH120817 Changed the names to remove reference to sim */
double spec_numin, spec_numax, spec_numean;	



int
spectral_estimators (xplasma)
     PlasmaPtr xplasma;
{
  int ireturn;
  double pl_alpha_min, pl_alpha_max, pl_alpha_temp, pl_w_temp, j;
  double exp_temp_min,exp_temp_max; /*120817 the 'temperature' range we are going to search for an effective temperature for the exponential model */
  double  exp_temp_temp, exp_w_temp; /*120817 the temporary values for temperature and weight of the exponential model */  
  int n;
  double pl_sd,exp_sd; /*120817 Computed standard decviations for two models for comparison with true value */
  double ALPHAMAX=20.0; /*120817 Something to make it a bit more obvious as to what values of the powerlaw exponent we consider reasonable */
  int plflag,expflag; /*120817 Two flags to say if we have a reasonable PL or EXP model, set to 1 initially, -1 means there has been some failure that means we must not use this model, +1 means it is OK */


/* This call is after a photon flight, so we *should* have access to j and ave freq, and so we can calculate proper values for W and alpha
To avoid problems with solving, we need to find a reasonable range of values within which to search for a solution to eq18. A reasonable guess is that it is around the current value....
*/



/* We loop over all of the bands, the first band is band number 0, and the last is band nxfreq-1 */
/* 71 - 111229 - ksl - Small modification to reflect moving nxfreq, etc in the geo structure */
  for (n = 0; n < geo.nxfreq; n++)

    {
    plflag=expflag=1;  //Both potential models are in the running
      if (xplasma->nxtot[n] == 0)
	{
	  Error
	    ("power_abundances: no photons in band %d which runs from %10.2e(%8.2fev) to %10.2e(%8.2fev) \n",n,geo.xfreq[n],geo.xfreq[n]*HEV,geo.xfreq[n+1],geo.xfreq[n+1]*HEV);
//	  spec_numin = xband.f1[0];	/*NSH 1108 Use the lower bound of the lowest band for sumnumin */
//	  spec_numax = xband.f2[xband.nbands - 1];	/*NSH 1108 and the upper bound of the upperband for max */
//	  spec_numean = xplasma->ave_freq;
//	  j = 0;		// If there are no photons then our guess is that there is no flux in this band
	  xplasma->pl_w[n] = 0;	//We also want to make sure that the weight will be zero, this way we make sure there is no contribution to the ionization balance from this frequency.
	  xplasma->pl_alpha[n] = 999.9; //Give alpha a value that will show up as an error
	  xplasma->exp_w[n] = 0.0;  //Make sure that w is zero, s no chance of mucking up ionization balance
	  xplasma->exp_temp[n] = 1e99;  //Give temp a value that will show up as an error
	  xplasma->spec_mod_type[n] = -1; //This tells the code that we have failed to model the spectrum in this band/cell/
	}
      else
	{
	  spec_numin = geo.xfreq[n];	/*1108 NSH n is defined in python.c, and says which band of radiation estimators we are interested in using the for power law ionisation calculation */
	  if (xplasma->max_freq < geo.xfreq[n + 1])
		{
		Log ("NSH resetting max frequency of band %i from %e to %e due to lack of photons\n",n,geo.xfreq[n+1],xplasma->max_freq);
		spec_numax = xplasma->max_freq;
		}
	  else
		{
	  	spec_numax = geo.xfreq[n + 1];
		}
	  spec_numean = xplasma->xave_freq[n];
	  j = xplasma->xj[n];
	

      Log
	("NSH We are about to calculate w and alpha, j=%10.2e, mean_freq=%10.2e, numin=%10.2e(%8.2fev), numax=%10.2e(%8.2fev), number of photons in band=%i\n",
	 j, spec_numean, spec_numin, spec_numin*HEV, spec_numax ,spec_numax*HEV, xplasma->nxtot[n]);


      /*1108 NSH ?? this could be a problem. At the moment, it relies on sim_alpha being defined at this point. */
      /*We should initialise it somewhere so that it *always* has a value, not just when the power law */
//	printf ("NSH We are starting with alpha=%f\n",xplasma->pl_alpha[n]);
      pl_alpha_min =  - 0.1; /*Lets just start the search around zero*/
      pl_alpha_max =  + 0.1;

      while (pl_alpha_func (pl_alpha_min) * pl_alpha_func (pl_alpha_max) > 0.0)
	{
	  pl_alpha_min = pl_alpha_min - 1.0;
	  pl_alpha_max = pl_alpha_max + 1.0;
	}
//	printf ("NSH PL alpha bracketed between %f (%e) and %f (%e)\n",pl_alpha_min,pl_alpha_func(pl_alpha_min),pl_alpha_max,pl_alpha_func(pl_alpha_max));
      if (sane_check(pl_alpha_func(pl_alpha_min)) || sane_check(pl_alpha_func(pl_alpha_max)))
	{
	Error ("spectral_estimators:sane_check Alpha cannot be bracketed in band %i cell %i- setting w to zero\n",n,xplasma->nplasma);
	xplasma->pl_w[n]=0.0;
	xplasma->pl_alpha[n]=-999.0; //Set this to a value that might let us diagnose the problem
	plflag=-1; 
	} 
      else
	{
/* We compute temporary values for sim alpha and sim weight. This will allow us to 
 * check that they are sensible before reassigning them */

	pl_alpha_temp = zbrent (pl_alpha_func, pl_alpha_min, pl_alpha_max, 0.00001);

      	if (pl_alpha_temp > ALPHAMAX)
		pl_alpha_temp = ALPHAMAX;	//110818 nsh check to stop crazy values for alpha causing problems
      	if (pl_alpha_temp < -1.*ALPHAMAX)
		pl_alpha_temp = -1.*ALPHAMAX;

/*This next line computes the PL weight using an external function. Note that xplasma->j already 
 * contains the volume of the cell and a factor of 4pi, so the volume sent to sim_w is set to 1 
 * and j has a factor of 4PI reapplied to it. This means that the equation still works in balance. 
 * It may be better to just implement the factor here, rather than bother with an external call.... */
//	printf ("NSH calling pl w with alpha=%f, numin=%e numax=%e\n",pl_alpha_temp,spec_numin,spec_numax);
      	pl_w_temp = pl_w (j ,pl_alpha_temp, spec_numin, spec_numax);

      	if (sane_check (pl_w_temp))
		{
	  	Error
	    	("spectral_estimators:sane_check New PL parameters unreasonable, using existing parameters. Check number of photons in this cell\n");
		plflag=-1; // Dont use this model
		xplasma->pl_w[n] = 0.0;
		xplasma->pl_alpha[n] =-999.0;
		}
      	else
		{
		xplasma->pl_alpha[n] = pl_alpha_temp;
	  	xplasma->pl_w[n] = pl_w_temp;
		}
	}

    
	
      exp_temp_min =  xplasma->t_r* 0.9; /*Lets just start the search around the radiation temperature in the cell*/
      exp_temp_max =  xplasma->t_r* 1.1;
      while (exp_temp_func (exp_temp_min) * exp_temp_func (exp_temp_max) > 0.0)
	{
	  exp_temp_min = exp_temp_min * 0.9;
	  exp_temp_max = exp_temp_max * 1.1;
	}
//	printf ("NSH exp_temp bracketed between %f (%e) and %f (%e)\n",exp_temp_min,exp_temp_func(exp_temp_min),exp_temp_max,exp_temp_func(exp_temp_max));
      if (sane_check(exp_temp_func(exp_temp_min)) || sane_check(exp_temp_func(exp_temp_max)))
	{
	Error ("spectral_estimators:sane_check Exponential temperature cannot be bracketed in band %i - setting w to zero\n",n);
	xplasma->exp_w[n]=0.0;
	xplasma->exp_temp[n]=-1e99;
	expflag=-1; //Discount an exponential model
	} 
      else
	{
/* We compute temporary values for sim alpha and sim weight. This will allow us to 
 * check that they are sensible before reassigning them */

	exp_temp_temp = zbrent (exp_temp_func,exp_temp_min, exp_temp_max, 0.00001); /* Solve for the effective temperature */
	exp_w_temp = exp_w  (j, exp_temp_temp, spec_numin, spec_numax); /* Calculate the weight */


	  if (sane_check (exp_w_temp))
		{
	  	Error
	    	("spectral_estimators:sane_check New exponential parameters unreasonable, using existing parameters. Check number of photons in this cell\n");
		expflag=-1; //discount an exponential model
		xplasma->exp_w[n] = 0.0;
		xplasma->exp_temp[n] = -1e99;
		}
      	else
		{
		xplasma->exp_temp[n] = exp_temp_temp;
	  	xplasma->exp_w[n] = exp_w_temp;
		}
	}
	
       	exp_sd=exp_stddev (xplasma->exp_temp[n],spec_numin,spec_numax);
        pl_sd=pl_stddev (xplasma->pl_alpha[n], spec_numin, spec_numax);
	Log ("NSH in this cell band %i PL estimators are w=%10.2e, alpha=%5.3f giving sd=%e compared to %e\n",n,xplasma->pl_w[n],xplasma->pl_alpha[n],pl_sd,xplasma->xsd_freq[n]); 
	Log ("NSH in this cell band %i exp estimators are w=%10.2e, temp=%10.2e giving sd=%e compared to %e\n",n,xplasma->exp_w[n],xplasma->exp_temp[n],exp_sd,xplasma->xsd_freq[n]); 
	exp_sd=fabs((exp_sd-xplasma->xsd_freq[n])/xplasma->xsd_freq[n]);
	pl_sd=fabs((pl_sd-xplasma->xsd_freq[n])/xplasma->xsd_freq[n]);
	/* NSH 120817 These commands decide upon the best model, based upon how well the models predict the standard deviation */
	if (expflag > 0 && plflag > 0)
		{        
		if (exp_sd < pl_sd) xplasma->spec_mod_type[n] = SPEC_MOD_EXP;
		else xplasma->spec_mod_type[n] = SPEC_MOD_PL;
		}
	else if (plflag > 0)
		{
		xplasma->spec_mod_type[n] = SPEC_MOD_PL;
		}
	else if (expflag > 0)
		{
		xplasma->spec_mod_type[n] = SPEC_MOD_EXP;
		}
	else
		{
		xplasma->spec_mod_type[n] = -1; //Oh dear, there is no suitable model
		Error ("No suitable model in band %i cell %i\n",n,xplasma->nplasma);
		}
	Log ("NSH In cell %i, band %i, the best model is %i\n",xplasma->nplasma, n,xplasma->spec_mod_type[n]);
	}  //End of loop that does things if there are more than zero photons in the band.
    }  //End of loop over bands	

 

  return (ireturn);

}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	integrand used with zbrent to calulate alpha from
	a mean frequency and a minimum and maximum frequecy
	for a band
        The function is the integral of nu j(nu) divided by the integral of j(nu).
	The prefector w is lost in the division, so this function should equal the
	mean frequency. Zbrent searches for the value of alpha which makes this
	function equal to the mean frequency computed for the cell from photon 
	packets. So we subtract the mean frequency at the end to make the value of the
	function equal to zero at the correct alpha.
	
 Arguments:		

Returns:
 
Description:	
	
Notes:

History:  NSH 120817 - Tidied up to make it a bit more clear what is going on here. There are
	now two functions, one is pl_mean, which is used elswhere, and ome is pl_alpha_func which
	is very simple, just works out the difference between the computed mean, and the measured
	mean


**************************************************************/

double
pl_alpha_func (alpha)
     double alpha;
{
  double answer;

  answer = pl_mean(alpha,spec_numin,spec_numax) - spec_numean;
 //    printf("NSH alpha=%.3f,f1=%10.2e,f2=%10.2e,meanfreq=%10.2e,ans=%e\n",alpha,spec_numin,spec_numax,spec_numean,answer);
  return (answer);
}

 double
pl_mean (alpha,numin,numax)
	double alpha;
	double numin,numax;
{
 double answer,numerator,denominator;

  numerator =      (pow(numax,(alpha+2.)) - (pow(numin,(alpha+2.))))  /  (alpha+2.); /*NSH 120817 This is the integral of (nu Bnu) */
  denominator =    (pow(numax,(alpha+1.)) - (pow(numin,(alpha+1.))))  /  (alpha+1.); /*NSH 120817 This is the integral of (nu Bnu) */
  answer = numerator/denominator;
  return (answer);
}


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

	This program works out the scaling factor for a power law representation of a spectrum in 
	a cell. Assuming alpha is known, it computes w for E=w x integral from numin to numax of nu^alpha.
	en1 is the band limited mean radiation density. 
        

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
111227	ksl	This routine was written by Nick in 2011, but has very
		little documentations or explanation
120817  nsh     routine moved into its natural home where it is used to 
		compute the prefactor for a power law estimation of a
		spectrum in a cell. 

 ************************************************************************/


double
pl_w (j, alpha, numin, numax)
     double j;		//the band limited spectral density
     double alpha;		//Computed spectral index for the cell
     double numin, numax;	//Range of frequencies we are considering
{
  double w;			//the answer

  double integral;  /*This will hold the unscaled integral of jnu from numin to numax */

  integral=(pow(numax,(alpha+1.0))-pow(numin,(alpha+1.0)))/(alpha+1.0); /* This is the integral */
  w=j/integral;  /*THe scaling factor is the actual band limited J / the integral. */

  return (w);
}




/**************************************************************************
                    Southampton University


  Synopsis:  

	This program works out the standard deviation of a photon distribution modelled
	by a power law. First we work out the integral of nu^2 jnu / int jnu.
	Then we subtract the (mean/j)^2 and then take the square root.
        

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:
120817  nsh     coded

 ************************************************************************/


double
pl_stddev (alpha, numin, numax)
     double alpha;		//Computed spectral index for the cell
     double numin, numax;	//Range of frequencies we are considering
{
  double answer;			//the answer

  double numerator;  /*This will hold the unscaled integral of nu^2 jnu from numin to numax */
  double denominator; /*This is the unscaled integral of jnu from numin to numax */
  double mean;

  numerator=(pow(numax,(alpha+3.0))-pow(numin,(alpha+3.0)))/(alpha+3.0); /* This is the integral for nu^2 bnu*/
  denominator=(pow(numax,(alpha+1.0))-pow(numin,(alpha+1.0)))/(alpha+1.0); /* This is the integral for nbu*/
  mean=pl_mean(alpha,numin,numax);
  answer=sqrt((numerator/denominator)-(mean*mean));


  return (answer);
}


/***********************************************************
                                       Southampton University

 Synopsis:
 	integrand used with zbrent to calulate effective temperature for an
	exponential model from a mean frequency and a minimum and maximum frequecy
	for a band
        The mean is the integral of nu j(nu) divided by the integral of j(nu).
	The prefector w is lost in the division, so this function should equal the
	mean frequency. Zbrent searches for the value of 'temperature' which makes this
	function equal to the mean frequency computed for the cell from photon 
	packets. So we subtract the mean frequency at the end to make the value of the
	function equal to zero at the correct temperature.
	
 Arguments:		

Returns:
 
Description:	
	
Notes:

History:  NSH 120817 - Tidied up to make it a bit more clear what is going on here. There are
	now two functions, one is pl_mean, which is used elswhere, and ome is pl_alpha_func which
	is very simple, just works out the difference between the computed mean, and the measured
	mean


**************************************************************/

double
exp_temp_func (exp_temp)
     double exp_temp;
{
  double answer;

  answer = exp_mean(exp_temp,spec_numin,spec_numax) - spec_numean;
 //    printf("NSH alpha=%.3f,f1=%10.2e,f2=%10.2e,meanfreq=%10.2e,ans=%e\n",alpha,spec_numin,spec_numax,spec_numean,answer);
  return (answer);
}

 double
exp_mean (exp_temp,numin,numax)
	double exp_temp;
	double numin,numax;
{
 double answer,numerator,denominator;
 double exp1; /* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt*/
 double emin,emax; /*The exponential evaluated at numin and numax */
 double t1,t2;

  exp1=(-1.0*H)/(BOLTZMANN*exp_temp);

	t1=exp1*numin;
	t2=exp1*numax; 


  emin = exp(t1);
  emax = exp(t2);


  numerator = ((emax/exp1)*((exp1*numax)-1.0))-((emin/exp1)*((exp1*numin)-1.0)); /*THe integral of nu e^(-hnu/kt) from numin to numax*/


  denominator =    emax-emin ;  /*NSH 120817 This is the integral of (J(nu)=e^(-hnu/jkt)) */
  answer = numerator/denominator; /*NB there is a factor of exp1 top and bottom which I have cancelled out */
  return (answer);
}


/**************************************************************************
                    Southampton university


  Synopsis:  

	This program works out the scaling factor for an exponential representation of a spectrum in 
	a cell. Assuming the effective tremperature is known, it computes w for E=w x integral from numin to numax of e^(-nhu/kt).

        

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:


 

  History:

120817  nsh     coded 

 ************************************************************************/


double
exp_w (j, exp_temp, numin, numax)
     double j;		//the band limited spectral density
     double exp_temp;		//Computed effective temperature for the cell
     double numin, numax;	//Range of frequencies we are considering
{
  double w;			//the answer

 double integral;  /*This will hold the unscaled integral of jnu from numin to numax */
 double exp1; /* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt*/
 double emin,emax; /*The exponential evaluated at numin and numax */

  exp1=(-1.0*H)/( BOLTZMANN*exp_temp);
  
 
  emin = exp(exp1*numin);
  emax = exp(exp1*numax);
  integral=(emax-emin)/exp1 ; /* This is the integral */
  w=j/integral;  /*THe scaling factor is the actual band limited J / the integral. */

  return (w);
}




/**************************************************************************
                    Southampton University


  Synopsis:  

	This program works out the standard deviation of a photon distribution modelled
	by an exponential. First we work out the integral of nu^2 jnu / int jnu.
	Then we subtract the (mean/j)^2 and then take the square root.
        

  Description:	

  Arguments:  	


  Returns:

  Notes:


 

  History:
120817  nsh     coded

 ************************************************************************/


double
exp_stddev (exp_temp, numin, numax)
     double exp_temp;		//Computed spectral index for the cell
     double numin, numax;	//Range of frequencies we are considering
{
  double answer;			//the answer
 double exp1; /* We supply a temperature, but actually we expect the correct function to be of the form e^-hnu/kt, so this will hold -1*h/kt*/
 double emin,emax; /*The exponential evaluated at numin and numax */
 double nmax,nmin; /*The integral is a bit complicated, so we have split it up into value for numin and numax */
 double numerator,denominator;
 double mean;

  exp1=(-1.0*H)/( BOLTZMANN*exp_temp);
  emin = exp(exp1*numin);
  emax = exp(exp1*numax);

  nmax=emax*(((numax*numax)/(exp1))-((2*numax)/(exp1*exp1))+((2.0)/(exp1*exp1*exp1))); /* Is the integral for nu^2 bnu exaluated for numax*/
  nmin=emin*(((numin*numin)/(exp1))-((2*numin)/(exp1*exp1))+((2.0)/(exp1*exp1*exp1))); /* Is the integral for nu^2 bnu exaluated for numin*/
  numerator=nmax-nmin; /*The definite integral */
  denominator= (emax-emin)/exp1; /* This is the integral for nbu*/
  mean=exp_mean(exp_temp,numin,numax);
  answer=sqrt((numerator/denominator)-(mean*mean));


  return (answer);
}

