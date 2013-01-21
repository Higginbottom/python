#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

int
sim_pl (xplasma)
     PlasmaPtr xplasma; 			
{
/*	PhotoionizationPtr x_ptr; 	*/
/*	int first,last; */
	double freq1,freq2,weight,t_e;
/*	double integ_vkpl(); */
/*	for (nion = 0; nion < nions; nion++)
		{
		kappa_ion[nion] = 0;
		frac_ion[nion] = 0;
	    	} 
*/
freq1=1e14;
freq2=1e16;	
      t_e = xplasma->t_e;
      weight = xplasma->w;

/*printf("We have %i elements",nelements);

/*  for (nion = 0; nion < nions; nion++)
    {
	pow_integ=integ_vkpl(x_ptr,freq1,freq2);
	printf("the integral is %e",pow_integ);
    }
*/
	printf ("weight = %e",weight);
	printf ("t_e=%e",t_e);	
	printf ("We got to the new code OK - nothing here tho\n");



return(0);
}




