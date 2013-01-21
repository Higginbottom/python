
/***********************************************************
                                       Southampton

 Synopsis:
	 open_diagfile sets up diagnostic files for use in python when the extra.diagnostics flag is 
Arguments:

	PhotPtr p;	the photon
	double ds	the distance the photon has travelled in the cell

Returns:
	Always returns 0.  .
 
Description:	 
Notes:
	

History:

	12jun 	nsh	72 Added lines to set up a file for outputting photons in given cells. The cells are read in from a file called diag_cells.dat. The existance of this file defined wether the diagnostic takes place.

**************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

int eplinit = 0;
int pstatinit = 0; /*To say if we have checked to see if we need to log photons */

int
open_diagfile ()
{
FILE *cellfile; /*File that may or may not exist, pointing to cells we want to write out photon stats for*/
int cell; /*Temporary storage of cell to use */

  if (eplinit == 0)
    {
      epltptr = fopen ("python.ext", "w");
      eplinit = 1;
    }

  ncstat=0; /*Zero the counter for the number of cells to be tracked */
  if (pstatinit == 0) /* Check we havent already done this */
	{
	cellfile = fopen ("diag_cells.dat","r"); /*This is the file containing cells to track*/
	if (cellfile!=NULL) /*If there actually *is* a file read it */
		{
		while (fscanf (cellfile, "%d", &cell) == 1) /*If the line contains only one integer number read it in, otherwise quit reading */
			{
			printf ("We have a cell - %i, ncstat=%i, NCSTAT=%i\n",cell,ncstat,NCSTAT);
			if (-1 < cell && cell < geo.nplasma && ncstat<NCSTAT) /*if the cells are real */
				{
				printf ("Accepted\n");
				ncell_stats[ncstat]=cell;
				ncstat=ncstat+1;
				}
			else
				{
				Error("open_diagfile: %i is an unacceptable cell number for photon tracking\n",cell);
				}
			}
		fclose(cellfile);
		pstatptr=fopen("cell_phot_stats.dat","w");
		}
	else
		{
		Log("open_diagfile: We have no file of cells to track, so we wont be doing any cell tracking\n");
		}
	pstatinit = 1; /* We have initialised this routine */
	}
	
//  diag_on_off = 0;            // 0=off everything else is on
  return (0);
}
