


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	partition_functions calculates the partition functions
	for a single ceell in the grid

Arguments:


Returns:
	
Description:

Mode here is identical to that used by nebular_concentrations, e.g

0 = LTE using t_r
1 = LTE using t_e
2 = Lucy and Mazzali



Notes:

	0800802 - This is a rewritten version of two routines, do_partitions
	and partition.   It is still true that levels is almost identical
	to some of the code here and therefore it is unclear why it is needed.
	levels is always called from partition_functions


History:
	080802	ksl	60b -- Brought levels routine into
			do_partitions since they are always
			called one after the other.  The routines
			are almost identical and it remains
			unclear why this was coded twice
	080802	ksl	Combined two routines do_partitions
			and partition into one in order to
			make this routine very similar to
			levels
	080804	ksl	Removed hubeny calculation of g from
			this code.  The advantage of the 
			Hubeny calculation was that it had
			a correction for density that we
			no longer had, but it was not being
			accessed by the code at all.

**************************************************************/



int
partition_functions (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int nion;
  double partition ();
  double t, weight;

  int n, m;
  int m_ground;			/* added by SS Jan 05 */
  double z, kt;


  if (mode == 0)
    {
      //LTE using t_r
      t = xplasma->t_r;
      weight = 1;
    }
  else if (mode == 1)
    {
      //LTE using t_e
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == 2)
    {
      //Non LTE calculation with radiative weights
      t = xplasma->t_r;
      weight = xplasma->w;
    }
  else
    {
      Error ("partition_functions: Unknown mode %d\n", mode);
      exit (0);
    }

  /* Calculate the partition function for each ion in turn */
  kt = BOLTZMANN * t;

  for (nion = 0; nion < nions; nion++)
    {

      if (ion[nion].nlevels > 0)
	//Calculate data on levels using a weighed BB assumption
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //Makes explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlevels; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else if (ion[nion].nlte > 0)
	//Calculate using "non-lte" levels
	{
	  m = ion[nion].first_nlte_level;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //This statement makes an explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlte; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else
	{
	  z = ion[nion].g;
	}


      xplasma->partition[nion] = z;
    }

  levels (xplasma, mode);	//levels from reduced BB

  return (0);
}
