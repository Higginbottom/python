#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       University of Southampton

Synopsis:
	matom_emiss_report is a routine which reports the matom level and kpkt emissivities summed
	over all cells

Arguments: 

Returns: 
	simply logs the level emissivites and kpkt emissivites in macro atom mode.
 
Description:	 
	prints out the level emissivities in macro atoms, summed over each cell. 
    Should only be used in macro atom mode in the spectral cycles 
    (when geo.matom_radiation = 1). 
	
History:
	130609 JM started routine 

**************************************************************/

int matom_emiss_report()
{
  //FILE *matomfile;
  int n, m;
  double emiss_sum, abs_sum;

  /* Cycle over macro atom levels and log emissivities */
  emiss_sum=0.0;
  abs_sum=0.0;

  for (m = 0; m < nlevels_macro + 1; m++)
  {
    for (n = 0; n < NPLASMA; n++)
    {
      emiss_sum += macromain[n].matom_emiss[m];
      abs_sum += macromain[n].matom_abs[m];
    }
    
    Log("MACRO LEVEL EMISSIVITIES (summed over cells): n %d matom_abs %8.4e matom_emiss %8.4e\n",
         n, abs_sum, emiss_sum); 
  }
 


  /* Log kpkt emissivities as well */
  emiss_sum=0.0;
  abs_sum=0.0;

  for (n = 0; n < NPLASMA; n++)
  {
    emiss_sum += plasmamain[n].kpkt_emiss;
    abs_sum += plasmamain[n].kpkt_abs;
    //Log("KPKT EMISSIVITIES: n %d macromain[n].kpkt_abs %8.4e macromain[n].kpkt_emiss %8.4e\n",
    //     n, plasmamain[n].kpkt_abs, plasmamain[n].kpkt_emiss);
  }

  Log("KPKT EMISSIVY (summed over cells): kpkt_abs %8.4e kpkt_emiss %8.4e\n",
       abs_sum, emiss_sum); 

  /* Log totals */ 
  Log("TOTALS: f_matom %le f_kpkt %le\n", geo.f_matom, geo.f_kpkt);
  
  return (0);
}

/***********************************************************
                                       University of Southampton

Synopsis:
	matom_trap is a simple trap and exit for debugging

**************************************************************/

int matom_trap()
{
  Log("matom_trap!! EXITING...\n");
  exit(0);
  return (0);
}
     
