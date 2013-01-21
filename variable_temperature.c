

/* This file was created to contain all of the routines associated
	with the attempt to calculate ionic abundances using a
	temperature fixed for each pair of ions to make the uncorrected
	abundance ration roughly equal to 1. There are then correction
	factors applied, either based upon a power law model of
	the true radiation field, or as a dilute blackbody.

	12Feb	nsh	Began coding 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"



struct photoionization *xver;	//Verner & Ferland description of a photoionization x-section
struct topbase_phot *xtop;	//Topbase description of a photoionization x-section
PlasmaPtr xxxplasma;
double qromb_temp;			//This is a storage variable for the current electron temperature so it is available for qromb calls
int niterate;			//Make this a variable that all the subroutines cn see, so we can decide if we need to recompute the numerators




/***********************************************************
                                       Southampton University

  Synopsis:   

int
variable_temperature (xplasama, mode)  modifies the densities of ions, levels, and
	partition functions of ions within a cell of the wind based upon the mode, 
	and other data contained within the WindPtr itself, such as t_r, t_e, w, 
	based on the "mode".  

        Unlike nebular_concentrations, this code works on pairs of ions at a time, 
	and uses a temperature calculated to be suitable for that particular pair
	to have a ratio of abundances about equal to 1.
  
  Arguments:		
     PlasmaPtr ww;
     int mode;			// 6=correct using a dilute blackbody
     				   7=correct using a broken power law spectrum


  Returns:
 	variable temperature alters portions of the wind ptr.  Exactly how things 
	are changed depends on the mode.
 
 	variable_temperature returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
  Description:	


 
  Notes:



  History:
	2012Feb	nsh	Coded and debugged as part of QSO effort. 
        1212Dec nsh	Recoded so that the densities are computed in a temporary array, and only 
			committted to the real density structure once we are sure the code converges.
	

**************************************************************/


//#define SAHA 4.82907e15		/* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
//#define MAXITERATIONS	200
//#define FRACTIONAL_ERROR 0.03
//#define THETAMAX	 1e4
//#define MIN_TEMP         100. //0712 moved into python.h

double xxxne,xip;

int
variable_temperature (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   6=correct using dilute blackbody, 7=power law
{
  int nion; // niterate; moved outside the code so all routines can see it
  double xnew, xsaha;
  double theta, x;
  double get_ne();
  double t, t_e,t_r, xtemp, nh, xne, xxne, www;
  double a,b;
  double newden[NIONS]; //NSH 121217 - recoded so density is computed in a temperary array
 // double xip; //ionzation potential of lower ion.
  int nelem,first,last;
  double t_r_part_correct,t_e_part_correct;
  double sum,big;
  double pi_fudge,recomb_fudge,gs_fudge,tot_fudge; /*The three correction factors for photoionization rate, recombination to ground state, and recombination rate */
  xxxplasma = xplasma;    /*Copy xplasma to local plasma varaible, used to communicate w and alpha to the power law correction routine. NSH 120703 - also used to set the denominator calculated for a given ion for this cell last time round*/


  nh = xplasma->rho * rho2nh;	//LTE
  t_e = xplasma->t_e; 
  t_r = xplasma->t_r;
  www = xplasma->w;
//		printf ("HERE WE ARE IN VT t_e=%e, t_r=%e, www=%e\n",t_e,t_r,www);



/* Copy the current densities into the temperary array */

  for (nion = 0; nion < nions; nion++)
	{
 	newden[nion]=xplasma->density[nion];
	}



  /* make an initial estimate of ne based on H alone,  Our guess
     assumes ion[0] is H1.  Note that x below is the fractional
     ionization of H and should vary from 0 to 1

     Note that in a pure H plasma, the left hand side of the Saha 
     equation can be written as (x*x*nh*nh)/((1-x)*nh)  and can
     then be converted into a quadratic.  That is what is done
     below.  

     Since this is an initial estimate g factors have been ignored
     in what follows.
   */
  //  t_e=7.6753e3;
  if (t_e < MIN_TEMP)
    t_e = MIN_TEMP;		/* fudge to prevent divide by zeros */

  xsaha = SAHA * pow (t_e, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t_e)) / nh;

  if (theta < THETAMAX)
    {
      x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
      xne = xxne = xxxne = x * nh;
    }
  else
    xne = xxne = xxxne = nh; /*xxne is just a store so the error can report the starting value of ne. xxxne is the shared variable so the temperature solver routine can access it*/

  if (xne < 1.e-6)
    xne = xxxne = 1.e-6;	   /* fudge to assure we can actually calculate
				   xne the first time through the loop */
  /* At this point we have an initial estimate of ne. */

  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
       /* We loop over all elements */
    for (nelem = 0; nelem < nelements; nelem++)
    	{
      	first = ele[nelem].firstion;	/*first and last identify the position in the array */
      	last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */
      	if (first < 0 || first >= nions || last < 0 || last > nions)
		{
	  	Error
	    	("variable_temperature: Confusion for element %d with first ion %d and last ion %d\n",
	     nelem, first, last);
	  	exit (0); 
		}
	/* and now we loop over all the ions of this element */
      	sum = newden[first] = 1.0; /* set the density of the first ion of this element to 1.0 - this is (always??) the neutral ion */
      	big = pow (10., 250. / (last - first)); /*make sure that we will not overflow the last ion */

      	for (nion = first + 1; nion < last; nion++)  /*nion is the upper ion of the pair, last is one above the actual last ion, so the last ion to be considered is last-1 */
		{
		/* now we need to work out the correct temperature to use */
		xip=ion[nion-1].ip;  //the IP is that from the lower to the upper of the pair
		xtemp=zbrent(temp_func,MIN_TEMP,1e8,10);  //work out correct temperature

/* given this temperature, we need the pair of partition functions for these ions */		
		partition_functions_2 (xplasma, nion, xtemp, 1); //weight of 1 give us the LTE populations.
/* and now we need the saha equation linking the two states at our chosen temp NB the partition function of the electron is included in the SAHA factor*/
		xsaha = SAHA * pow (xtemp, 1.5);
	  	b = xsaha * xplasma->partition[nion]
	    	* exp (-ion[nion - 1].ip / (BOLTZMANN * xtemp)) / (xne *
							   xplasma->partition[nion-1]);
//printf ("PART FUNC FOR 	for element %i, ion %i and %i, xtemp= %e is %f and %f\n",ion[nion-1].z,ion[nion-1].istate,ion[nion].istate,xtemp,partition[nion-1],partition[nion]);
		t_e_part_correct=xplasma->partition[nion-1]/xplasma->partition[nion];
		partition_functions_2 (xplasma,nion, t_e, 0); //Our only real guess here is that the electron temperature might give a good estimate of the partition function
		t_e_part_correct*=(xplasma->partition[nion]/xplasma->partition[nion-1]);
//printf ("PART FUNC FOR 	for element %i, ion %i and %i, dil_tr=%e is %f and %f\n",ion[nion-1].z,ion[nion-1].istate,ion[nion].istate,t_r,partition[nion-1],partition[nion]); 
//		
//printf ("PART FUNC FOR 	for element %i, ion %i and %i, t_e=   %e is %f and %f\n",ion[nion-1].z,ion[nion-1].istate,ion[nion].istate,t_e,partition[nion-1],partition[nion]);
//		t_e_part_correct*=(partition[nion]/partition[nion-1]);
//		printf ("PART_CORRECT t_r %f for element %i, upper ion %i, xtemp=%e, t_r=%e  \n",t_r_part_correct,ion[nion].z,ion[nion].istate,xtemp,t_r);
/* we now correct b to take account of the temperature and photon field 
		t_r and www give the actual radiation field in the cell, xtemp is the temp we used
  		t_e is the actual electron temperature of the cell*/

	      	recomb_fudge = sqrt (t_e / xtemp);
         	if (mode == 6) /* Correct the SAHA equation abundance pair using an actual radiation field modelled as a dilute black body */
			{
			if (t_r/xtemp <0.2) /*If the true radiation temperature is much lower than the temperature at which the ion is expected to exist, we wont see it, so save time and dont bother calculting correction factors*/
				{
				tot_fudge = 0.0;
	 			}
			else
				{
				pi_fudge = bb_correct_2 (xtemp, t_r, www, nion);
				gs_fudge = compute_zeta (t_e, nion -1, 2); /* Calculate the ground state recombination rate correction factor based on the cells true electron temperature.  */
				tot_fudge=pi_fudge*recomb_fudge*(gs_fudge+www*(1-gs_fudge));
				}
			}
		else if (mode == 7) /* correct the SAHA equation abundance pair using an actual radiation field modelled as a broken power law */
			{
			pi_fudge = pl_correct_2 (xtemp, nion);
			gs_fudge = compute_zeta (t_e, nion -1, 2); /* Calculate the ground state recombination rate correction factor based on the cells true electron temperature.  */
			tot_fudge=pi_fudge*recomb_fudge*gs_fudge*t_e_part_correct;
			}

		/* apply correction factors */
		b *= tot_fudge;
		if (b > big)	 
	   	b = big;		//limit step so there is no chance of overflow
	  	a = newden[nion - 1] * b;
	  	sum += newden[nion] = a;
	  	if (newden[nion] < 0.0)
			{
			Error ("variable_temperature: ion %i has a negative density %e\n",nion,newden[nion]);
	    		mytrap ();
			}
	  	if (sane_check (sum))
			Error ("variable_temperature:sane check failed for running density sum=%e, last ion=%i\n",sum,nion);
		}  //end of loop over ions - we now have a full set of ions for this element

      	a = nh * ele[nelem].abun / sum;   //the scaling factor to get the overall abundance right
      	for (nion = first; nion < last; nion++)
		{
	  	newden[nion] *= a;  //apply scaling
	  	if (sane_check (newden[nion]))  //check nothing has gone crazy
			Error ("variable_temperature:sane check failed for density newden=%e, for ion=%i\n",newden[nion],nion);
		}

      	if (geo.macro_ioniz_mode == 1)
		{
  		macro_pops (xplasma, xne);
		}

      /*Set some floor so future divisions are sensible */
      	for (nion = 0; nion < nions; nion++)
		{
	  	if (newden[nion] < DENSITY_MIN)
	    		newden[nion] = DENSITY_MIN;
		}

      	/* Now determine the new value of ne from the ion abundances */
	}   //end of loop over elements
    xnew = get_ne (newden);	/* determine the electron density for this density distribution */
//	printf ("Solver, change in n_e = %e vs FRACTIONAL ERROR of %e in xne of %e\n",fabs ((xne - xnew) / (xnew)) , FRACTIONAL_ERROR,xne);
    if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */
    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	{
	break;
	}
    xne = xxxne = (xnew+xne)/2.;   /*New value of ne */
    niterate++;


    if (niterate == MAXITERATIONS)
    	{
      	Error ("variable_temperature: failed to converge t %.2g nh %.2g xnew %.2g\n",
	     t, nh, xnew);
      	Error ("variable_temperature: xxne %e theta %e\n", xxne, theta);
      	return (-1);
    	}
    }


/* Finally transfer the calculated densities to the real density array. This is only called if the code
 has iterated correctly, if not then the real densities will stay the same as last time*/

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nion] = newden[nion];
	}
    }





  partition_functions (xplasma, 4); /*WARNING fudge NSH 11/5/14 - this is as a test. We really need a better implementation of partition functions and levels for a power law illuminating spectrum. We found that if we didnt make this call, we would end up with undefined levels - which did really crazy things*/
  return (0);
}





/***********************************************************
                                       Southampton University

  Synopsis:   
	bb_correct_2(xtemp, t_r,www,nion) calculates the 
 	correction factor of the number density of two adjacent levels calculated
	at any temperature. LM is a special case of this where the temperatures
	are the same. 
 
  Arguments:		

	double xtemp		the temperature at which the saha abundances were calculated
	double t_r		temperature of the radiation field 
	double www		dilution factor for the radiation field
	int nion;		number of the upper ion in the pair for which the correction factor is 					being calculated.

  Returns:
 	j			the correction factor
 
 
 	
 Description:	


 
  Notes:



  History:
 
	12feb		NSH	Coded as part of the quasar project

**************************************************************/
int bb_correct_err=0;

double
bb_correct_2 (xtemp, t_r, www, nion)
     double t_r, xtemp, www;
     int nion;
{
  double q;
  int ion_lower,n;
  double numerator, denominator;
  int ntmin,nvmin;
  double fthresh,fmax,fmaxtemp;



  /* Initialization of fudges complete */

  ion_lower = nion-1;	/*the call uses nion, which is the upper ion in the pair */


  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;
      nvmin = ion_lower;
    }
  else				
    {
      Error ("bb_correct_2: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }

/*The integration is to correct for photoionization rates from the lower to the upper, so
 all integrals have to be over the lower ion cross section. The integrals use the code in
 power_sub, so the temperature has to be set to an external variable. Annoyingly this is 
 called sim_te....*/

 if (ion[ion_lower].phot_info == 1)  //topbase
    {
      n = ntmin;
      xtop = &phot_top[n];
      fthresh = xtop->freq[0];
      fmaxtemp = xtop->freq[xtop->np - 1];
      fmax = check_fmax (fthresh,fmaxtemp,xtemp);
      if (fthresh > fmax)
	{
	Error("bb_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
	q=1.0;
	return(q);
	}

      qromb_temp=t_r;  //The numerator is for the actual radiation temperature
      
     if (pow(((t_r-xxxplasma->PWntemp[ion_lower])/xxxplasma->PWntemp[ion_lower]),2)>1e-6)
	{
        numerator=www*qromb (tb_planck1, fthresh, fmax, 1.e-4); //and is corrected for W
        xxxplasma->PWnumer[ion_lower]=numerator;
	xxxplasma->PWntemp[ion_lower]=t_r;
	}
      else
	{	
	numerator = xxxplasma->PWnumer[ion_lower];
	}

      if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator


      qromb_temp=xtemp; //The denominator is calculated for the LTE rate at our ideal temp

//	printf ("topbase n=%i,nion=%i,temp=%e,fthresh=%e,fmax=%e,hnu=%e,hnu/kt=%e,fmax=%e\n",n,ion_lower,temp,fthresh,fmax,fthresh*H,(fthresh*H)/(temp*BOLTZMANN),5.879e10*temp);



      if (pow(((xtemp-xxxplasma->PWdtemp[ion_lower])/xxxplasma->PWdtemp[ion_lower]),2)>1e-6)
	{
        denominator=qromb (tb_planck1, fthresh, fmax, 1.e-4);	
        xxxplasma->PWdenom[ion_lower]=denominator;
	xxxplasma->PWdtemp[ion_lower]=xtemp;
	}
      else
	{
	denominator=xxxplasma->PWdenom[ion_lower];
	}


    }
 else if (ion[ion_lower].phot_info == 0)  // verner
    {			
      n = nvmin;		//just the ground state ioinzation fraction.
      xver = &xphot[ion[n].nxphot];
      fthresh = xver->freq_t;
      fmaxtemp=xver->freq_max;
      fmax = check_fmax (fthresh,fmaxtemp,xtemp);
      if (fthresh > fmax)
	{
	Error("bb_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
	q=1.0;
	return(q);
	}
      qromb_temp=t_r;
		


     if (pow(((t_r-xxxplasma->PWntemp[ion_lower])/xxxplasma->PWntemp[ion_lower]),2)>1e-6)
	{
        numerator = www*qromb (verner_planck1, fthresh, fmax, 1.e-4);
        xxxplasma->PWnumer[ion_lower]=numerator;
        xxxplasma->PWntemp[ion_lower]=t_r;
	}
      else
	{
	numerator = xxxplasma->PWnumer[ion_lower];
	}
       if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator
      qromb_temp=xtemp;
      fthresh = xver->freq_t;
      if (pow(((xtemp-xxxplasma->PWdtemp[ion_lower])/xxxplasma->PWdtemp[ion_lower]),2)>1e-6)
	{
        denominator = qromb (verner_planck1, fthresh, fmax, 1.e-4);	
        xxxplasma->PWdenom[ion_lower]=denominator;
	xxxplasma->PWdtemp[ion_lower]=xtemp;
	}
      else
	{
	denominator=xxxplasma->PWdenom[ion_lower];
	}
    }

 else
    {
	    bb_correct_err++;
	    /* 71 - 111229 - ksl - Suppress this error after 100 instances so program does not bomb */
	    if (bb_correct_err<100){
      Error
	("bb_correct: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	 nion,ion[nion].z,ion[nion].istate);
	    }
	    else if (bb_correct_err==100){
		    Error("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	    }

	   
      q = 1.;
      return (q);
    }

      q = numerator / denominator;

    
  return (q);
}

/* This next little bit of code was used to calculate an ideal temperature. This gets a temperature which gives a ratio of n to n+1 of close to 1, but at this temperature, the blackbody function used for the denominator of the correction factor is almost always zero. I've left the routine in, but it doesnt seem to work...*/



double
temp_func (solv_temp)
     double solv_temp;
{
  double answer;
//	printf ("xxxne=%e, xip=%e\n",xxxne,xip);
//	printf ("1=%e,2=%e,3=%e\n",log(4.83e15/xxxne),1.5*log(temp),(xip/(BOLTZMANN*temp)));
  answer = log(4.83e15/xxxne)+1.5*log(solv_temp)-(xip/(BOLTZMANN*solv_temp));

  return (answer);
}


/***********************************************************
                                       Southampton University

  Synopsis:   
	pl_correct_2(xtemp,nion) calculates the 
 	correction factor of the number density of two adjacent levels calculated
	at any temperature based upon a true radiation field modelled by a broken power law
 
  Arguments:		

	double xtemp		the temperature at which the saha abundances were calculated
	int nion;		the upper ion in an ion pair for which we want to correction factyor

  Returns:
 	q		        the correction factor
 	
 Description:	

	 

 
  Notes:
	The weights and alpha values needed for the integrals are communicated to this routine vias xxxplasma, an external varaible which is popullated in the variable_temperature routine.
	The integrals over power law require individual alphas and w, these are communicated via external variables xpl_alpha and xpl_w.	


  History:
 
	12feb	NSH	Coded as part of the quasar project
        12aug   NSH     Added code to allow an ex[ponential model of j to be used.

**************************************************************/
int pl_correct_err=0;
double xpl_alpha,xpl_w;
double xexp_temp,xexp_w;


double
pl_correct_2 (xtemp, nion)
     double xtemp;
     int nion;
{
  double q;
  int ion_lower,n,j;
  double numerator, denominator;
  int ntmin,nvmin;
  double fthresh,fmax,fmaxtemp;
  double f1,f2;
  double exp_qromb,pl_qromb;


exp_qromb=1e-4;
pl_qromb=1e-4;

  ion_lower = nion-1;	/*the call uses nion, which is the upper ion in the pair, we are interested in the PI from the lower ion */	

  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;  /*We only ever use the ground state cross sections. This is for topbase */
      nvmin = ion_lower; /*and this is for verner cross sections*/
    }
  else				
    {
      Error ("bb_correct_2: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }


/*The integration is to correct for photoionization rates from the lower to the upper, so
 all integrals have to be over the lower ion cross section. The integrals use external functions, so the temperature has to be set to an external variable temp. */

 if (ion[ion_lower].phot_info == 1)  //topbase
    {
      n = ntmin;
      xtop = &phot_top[n];
      fthresh = xtop->freq[0];
      fmax = xtop->freq[xtop->np - 1];

      numerator=0;
	  if (niterate==0)  			//first time of iterating this cycle, so calculate the numerator
              {
	      for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	        {
		xpl_alpha = xxxplasma->pl_alpha[j];
		xpl_w = xxxplasma->pl_w[j];
		xexp_temp = xxxplasma->exp_temp[j];
		xexp_w = xxxplasma->exp_w[j];
	//	if (ion_lower==0) printf ("BCA IN band %i, we will be using model %i, pl_a=%e, pl_w=%e, (%e-%e) exp_t=%e, exp_w=%e \n",j,xxxplasma->spec_mod_type[j],xxxplasma->pl_alpha[j],xxxplasma->pl_w[j],xxxplasma->pl_w[j]*pow(geo.xfreq[j],xxxplasma->pl_alpha[j]),xxxplasma->pl_w[j]*pow(geo.xfreq[j+1],xxxplasma->pl_alpha[j]),xxxplasma->exp_temp[j],xxxplasma->exp_w[j]);
                if (xxxplasma->spec_mod_type[j] > 0)   //Only bother doing the integrals if we have a model in this band
			{
			f1=geo.xfreq[j];
//			if (geo.xfreq[j+1] > xxxplasma->max_freq && geo.xfreq[j] < xxxplasma->max_freq) //The maximum frequency seen in this cell is in this band, so we cannot safely use the power law estimators right up to the top of the band. Note, we hope that all the weights in bands above this will be zero!
//				{		
//				f2=xxxplasma->max_freq;
//				}
//			else
//				{
				f2=geo.xfreq[j+1]; //We can safely integrate over the whole band using the estimators for the cell/band 
//				}
	        	if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)	//Case 1- 
		  		{
				if (xxxplasma->spec_mod_type[j]==SPEC_MOD_PL) 
					{
					numerator += qromb (tb_pow1, fthresh, fmax, 1.e-4); 
					}
				else 
					{
					numerator += qromb (tb_exp1, fthresh, fmax, exp_qromb);
//					printf ("BCA BB integrate (1)= %e,num=%e\n",qromb(tb_planck1, fthresh, fmax, 1.e-4),numerator);
					} 
		  		}
	        	else if (f1 < fthresh && fthresh < f2 && f2 < fmax)	//case 2 
		  		{
		  		if (xxxplasma->spec_mod_type[j]==SPEC_MOD_PL) 
					{
					numerator += qromb (tb_pow1, fthresh, f2, 1.e-4);
					}
				else 
					{
					numerator += qromb (tb_exp1,  fthresh, f2, exp_qromb);
//					printf ("BCA BB integrate (2)= %e, num=%e\n",qromb(tb_planck1, fthresh, f2, 1.e-4),numerator);
					}
		  		}
	        	else if (f1 > fthresh && f1 < fmax && fmax < f2)	//case 3
		  		{
		  		if (xxxplasma->spec_mod_type[j]==SPEC_MOD_PL) 
					{
					numerator += qromb (tb_pow1, f1, fmax, 1.e-4);
					}
				else 
					{
					numerator += qromb (tb_exp1,  f1, fmax, exp_qromb);
//					printf ("BCA BB integrate (3)= %e, num=%e\n",qromb(tb_planck1, f1, fmax, 1.e-4),numerator);
					}
		  		}
	        	else if (f1 > fthresh && f2 < fmax)	// case 4
		  		{
		  		if (xxxplasma->spec_mod_type[j]==SPEC_MOD_PL) 
					{
					numerator += qromb (tb_pow1, f1, f2, 1.e-4);
					}
				else 
					{
					numerator += qromb (tb_exp1, f1, f2, exp_qromb);
//					printf ("BCA BB integrate (4)= %e, num=%e\n",qromb(tb_planck1, f1, f2, 1.e-4),numerator);
					}
		  		}
	        	else		//case 5 - should only be the case where the band is outside the range for the integral.
		  		{
		  		numerator += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		  		}
//	if (ion_lower==0) printf ("BCA numerator = %e after band %i, which was model %i \n",numerator,j,xxxplasma->spec_mod_type[j]);
	        	} //End of loop to only integrate in this band if there is power
		} //End of loop over all bands, at this point we have the numerator
              xxxplasma->PWnumer[ion_lower]=numerator;   // Store the calculated numerator for this cell - it wont change during one ionization cycle
	      } //End of if statement to decide if this is the first iteration
	   else
	      {
	      numerator=xxxplasma->PWnumer[ion_lower];   // We are on the second iteration, so use the stored value
	      }
      if (numerator==0.0)
	    {
	    return (0.0); //There is no need to waste time evaluating the denominator
            }
     fmaxtemp = xtop->freq[xtop->np - 1];
      fmax = check_fmax (fthresh,fmaxtemp,xtemp);
      if (fthresh > fmax)
	{
	Error("pl_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
	q=1.0;
	return(q);
	}
      qromb_temp=xtemp; //The denominator is calculated for the LTE rate at our ideal temp. If we get this wrong, then divide by zeros abound!

      if (pow(((xtemp-xxxplasma->PWdtemp[ion_lower])/xxxplasma->PWdtemp[ion_lower]),2)>1e-6) //If our guess temperature hasnt changed much, use denominator from last time
	{
        denominator=qromb (tb_planck1, fthresh, fmax, 1.e-4);
        xxxplasma->PWdenom[ion_lower]=denominator;
	xxxplasma->PWdtemp[ion_lower]=xtemp;
	}
      else
	{
	denominator=xxxplasma->PWdenom[ion_lower];
	}


    }
 else if (ion[ion_lower].phot_info == 0)  // verner
    {			
      n = nvmin;		//just the ground state ionization fraction.

      xver = &xphot[ion[n].nxphot];
      fthresh = xver->freq_t;
      fmax=xver->freq_max;

//	if (ion[nion].z==26) printf ("VERNER ion%i, integrating from %e to %e \n",ion[nion].istate,fthresh,fmax);

	numerator=0;
	  if (niterate==0)  			//first time of iterating this cycle, so calculate the numerator
              {
              for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
		 {
		 xpl_alpha = xxxplasma->pl_alpha[j]; //Set the variables needed for integration
		 xpl_w = xxxplasma->pl_w[j];
		 xexp_temp = xxxplasma->exp_temp[j];
		 xexp_w = xxxplasma->exp_w[j];
                 if (xxxplasma->spec_mod_type[j] > 0)   //Only bother doing the integrals if we have a model in this band
			{
			f1=geo.xfreq[j];
//			if (geo.xfreq[j+1] > xxxplasma->max_freq && geo.xfreq[j] < xxxplasma->max_freq) //The maximum frequency seen in this cell is in this band, so we cannot safely use the power law estimators right up to the top of the band. Note, we hope that all the weights in bands above this will be zero!
//				{
//				f2=xxxplasma->max_freq;
//				}
//			else
//				{
				f2=geo.xfreq[j+1]; //We can safely integrate over the whole band using the estimators for the cell/band 
//				}
	      		if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)	//Case 1
		  		{
		    		if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL) 
					{
					numerator += qromb (verner_pow1, fthresh, fmax, 1.e-4);
					}
				else 
{				
					numerator += qromb (verner_exp1, fthresh, fmax, 1.e-4);
					
					}
		  		}
	        	else if (f1 < fthresh && fthresh < f2 && f2 < fmax)	//case 2 
		  		{
		    		if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL) numerator += qromb (verner_pow1, fthresh, f2, 1.e-4);
				else numerator += qromb (verner_exp1, fthresh, f2, 1.e-4);
		  		}
	        	else if (f1 > fthresh && f1 < fmax && fmax < f2)	//case 3 
		  		{
		    		if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL) numerator += qromb (verner_pow1, f1, fmax, 1.e-4);
				else numerator += qromb (verner_exp1, f1, fmax, 1.e-4);
		  		}
	        	else if (f1 > fthresh && f2 < fmax)	// case 4
		  		{
		    		if (xxxplasma->spec_mod_type[j] == SPEC_MOD_PL) numerator += qromb (verner_pow1, f1, f2, 1.e-4);
				else numerator += qromb (verner_exp1, f1, f2, 1.e-4);
		  		}
	        	else		//case 5 - should only be the case where the band is outside the range for the integral.
		  		{
		    		numerator += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		  		}

			} //End of if loop to only do any integrations if the band weight is non zero     	
		}  //End of loop over all the bands, we now have the integration
            xxxplasma->PWnumer[ion_lower]=numerator;  //Store the numerator fromn this time - it wont change at all during one ionization cycle	 
	    } //End of if statement to decide if this is the first iteration
	  else   //we have already done one interation, so we can just use the stored numberator
            {
            numerator = xxxplasma->PWnumer[ion_lower];
	    }

       if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator
/* Denominator is integral at LTE of our chosen temperature. */
      fmaxtemp=xver->freq_max;
      fmax = check_fmax (fthresh,fmaxtemp,xtemp);
	
      if (fthresh > fmax)
	{
	Error("pl_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
	q=1.0;
	return(q);
	}




      qromb_temp=xtemp;

     if (pow(((xtemp-xxxplasma->PWdtemp[ion_lower])/xxxplasma->PWdtemp[ion_lower]),2)>1e-6)
	{
        denominator = qromb (verner_planck1, fthresh, fmax, 1.e-4);
        xxxplasma->PWdenom[ion_lower]=denominator;
	xxxplasma->PWdtemp[ion_lower]=xtemp;
	}
      else
	{
	denominator=xxxplasma->PWdenom[ion_lower];
	}
    }	 //End of if statement for verner cross section

 else
    {
	    pl_correct_err++; /* If we get here, there are no cross sections available */
	    if (pl_correct_err<100){
      Error
	("pl_correct: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	 nion,ion[nion].z,ion[nion].istate);
	    }
	    else if (pl_correct_err==100){
		    Error("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	    }

	   
      q = 1.; /*This is really bad actually, this will leave the abundances all wrong.... */
      return (q);
    }

   q = numerator / denominator; /*Just work out the correction factor - hope it isnt infinite! */





  return (q);
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

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.

 ************************************************************************/


double
tb_planck1 (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot_topbase (xtop, freq);
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

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.

 ************************************************************************/


double
verner_planck1 (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot (xver, freq);
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

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
 ************************************************************************/



double
tb_pow1 (freq)
     double freq;
{
  double answer;

  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));
  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.

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

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
 ************************************************************************/


double
verner_pow1 (freq)
     double freq;
{
  double answer;

  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));
  answer *= sigma_phot (xver, freq);	// and finally multiply by the cross section.

  return (answer);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	verner_exp is the function to allow integration of an exponential photon distribution multiplied by the inoisation cross section
   	to allow the numerator of the sim correction factor to be calulated for ions with verner cross sections.

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Aug NSH - Written as part of the effort to improve the modelling
 ************************************************************************/


double
verner_exp1 (freq)
     double freq;
{
  double answer;

  answer = xexp_w * exp ((-1.0*H*freq)/(BOLTZMANN*xexp_temp));
  answer *= sigma_phot (xver, freq);	// and finally multiply by the cross section.
  answer /= freq;
  return (answer);
}




/**************************************************************************
                    Southampton University


  Synopsis:  

  Description:	
	tb_exp is the function to allow integration of an exponential photon distribution multiplied by the inoisation 
	cross section
  	to allow the numerator of the correction factor to be calulated for ions with topbase cross sections.

  Arguments:  


  Returns:

 Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Aug Written by NSH as part of the effort to improve spectral modelling
 ************************************************************************/



double
tb_exp1 (freq)
     double freq;
{
  double answer;
	
  answer = xexp_w * exp ((-1.0*H*freq)/(BOLTZMANN*xexp_temp));
  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.
  answer /= freq;
  return (answer);
}

