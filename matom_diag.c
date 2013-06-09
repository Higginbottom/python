#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       University of Southampton

Synopsis:
	matom_emiss_report is a routine which reports the matom level and kpkt emissivities 

Arguments: 

Returns: 
	simply logs the level emissivites and kpkt emissivites in macro atom mode.
 
Description:	 
	prints out the level emissivities in macro atoms. Should only be used in macro atom mode
	in the spectral cycles (when geo.matom_radiation =1). 
	
History:
	130609 JM started routine 

**************************************************************/

int matom_emiss_report()
{
  //FILE *matomfile;
  int n, m;

  /* Cycle over macro atom levels and log emissivities */
  for (n = 0; n < NPLASMA; n++)
  {
    for (m = 0; m < nlevels_macro + 1; m++)
    {
      Log("MACRO LEVEL EMISSIVITIES: n %d m %d macromain[n].matom_abs[m] %8.4e macromain[n].matom_emiss[m] %8.4e\n",
           n, m, macromain[n].matom_abs[m], macromain[n].matom_emiss[m]);
    }
  }


  /* Log kpkt emissivities as well */
  for (n = 0; n < NPLASMA; n++)
  {
    Log("KPKT EMISSIVITIES: n %d macromain[n].kpkt_abs %8.4e macromain[n].kpkt_emiss %8.4e\n",
         n, plasmamain[n].kpkt_abs, plasmamain[n].kpkt_emiss);
  }

  /* Log totals */ 
  Log("TOTALS: f_matom %le f_kpkt %le\n", geo.f_matom, geo.f_kpkt);
  
  return (0);
}
     



/***********************************************************
                                       University of Southampton

Synopsis:
	matom_emiss_report is a routine which reports the level emissivities 

Arguments: 

Returns: 
	prints out the level emissivities in macro atoms. 
 
Description:	 
	prints out the level emissivities in macro atoms. Should only be used in macro atom mode
	in the spectral cycles (when geo.matom_radiation =1). Mostly written 
	
History:
	130609 JM started routine to report level emissivities

**************************************************************/
/*
int matom_emiss_report()
{
  FILE *matomfile;


*/
