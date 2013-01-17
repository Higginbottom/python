
/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	py_wind is a program which can be used to display various parameters of a wind
		as calculated by python.  This is the main routine.

Arguments:		

	py_wind [-h] [-s] [root]

	where
		-h 	prints out a short help file and exits (see help routine below)
		-s	causes certain parameters to be printed out in individual files
			after which the program exits
		root	optional root name of wind_save file.  If this is not given,
			the user is queried for this


Returns:
 
Description:	
	
	Py_wind simply reads and then displays the wind file created by python.  It can
	select various parameters from the wind and write them to a file so that one
	can create displays of them.  
	
	The file can contain either the original gridding which was used by python, in which 
	case the file prefix will be "x.", or "z", in which case it will be regridded to a 
	linear array.  This option is intended so one can create a contour plot more easily.


		
Notes:

History:
 	97jun	ksl	Coding on py_wind began.
 	97sep13	ksl	Modified slightly to reflect the fact that EOF exits program through rdpar
 	97oct2	ksl	Modified so that t_e and t_r can be printed separately
 	98feb	ksl	Modified to add the amount of energy and momentum flux absorbed per cell
 				as output possibilities. Also in py_ion added the ability to print out the
 				luminosity in CIV if the wind was assumed to be a thin plasma.
	98jul	ck	Modified so the number of ionzations and recombinations can be checked
			against one another
	98oct	ksl	Modified abs_sum and lum_suim so that the various processes which go
				into heating and cooling the plasma can be displayed individually.
	99nov	ksl	Added an additional option to create a fixed set of output files. 
	00mar	ksl	Updated main routine and write_array subroutine to allow one to choose
			whether to relinearize the grid or to leave it in as it was in the windfile
			Added output files for velocity.
	01sep	ksl	Updated ion_summary so that either concentrations or ion fractions could be
			examined.
	02apr	ksl	Fixed so could read wind_save files with various sizes 
	02jun	ksl	Fixed and so could read diagnostic outputs when one had printed out 
			results at each iteration of the ionization cycle
	04nov	ksl     Made major modification to eliminate duplicative code, by adding new
			task display.  Also introduced variable py_wind_project to switch
			between raw mode (regardless of coord system), and a mode which projects
			on to a yz plane the system of coordinates.
	04nov	ksl	Added a choice to explore what is the situation at a specific position
			in the grid.  This was added to see how different coordinate systems
			compared.
	05jan	ksl	54f -- Modified			
	05apr	ksl	56 -- Eliminated MDIM from the routines in this file.  
			MDIM is isolated to routines like display.  Added swithches to
			change choices regarding what, if anything is printed to a 
			file.
	05jul	ksl	56d -- Moved all of the subroutines to a separate routine py_wind_sub.c
			This is intended to facilitate the create of templates by Makefile.
	06jul	ksl	57g -- Made modifications to reflect the change in the way structures
			had been defined.
	090117	ksl	68a -- Added command line functionality so could run more easily
			in a non interactive mode
	090125	ksl	Changed functionality of I and i switches to allow for printing
			out information about the number of scatters, etc, by an ion
			in the wind.  This reflects a desire to be able to understand
			better where in the wind ions are "active" in creating the
			spectrum.

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"



int
main (argc, argv)
     int argc;
     char *argv[];
{

  WindPtr w;

  int n, i, istate;
  int ochoice;
  char c;

  char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH],
    specfile[LINELENGTH];
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char photfile[LINELENGTH];
  double lambda, freq;
  int interactive;
  int iswitch;





  // py_wind uses rdpar, but only in an interactive mode. As a result 
  // there is no associated .pf file

  interactive = 1;		/* Default to the standard operating mofe for py_wind */

  if (argc == 1)
    {
      printf ("Root for wind file :");
      fgets (input, LINELENGTH, stdin);
      get_root (root, input);
    }
  else
    {
      for (i = 1; i < argc; i++)
	{
	  if (strcmp (argv[i], "-h") == 0)
	    {
	      py_wind_help ();
	    }
	  if (strcmp (argv[i], "-s") == 0)
	    {
	      interactive = 0;
	    }
	  else if (strncmp (argv[i], "-", 1) == 0)
	    {
	      Error ("py_wind: unknown switch %s\n", argv[i]);
	      py_wind_help ();
	    }
	}

      strcpy (input, argv[argc - 1]);
      get_root (root, input);
    }


  printf ("Reading data from file %s\n", root);

  /* Now create the names of all the files which will be written */

  strcpy (wspecfile, root);
  strcpy (specfile, root);
  strcpy (windradfile, root);
  strcpy (windsavefile, root);
  strcpy (photfile, root);

  strcat (wspecfile, ".spec_tot");
  strcat (specfile, ".spec");
  strcat (windradfile, ".wind_rad");
  strcat (windsavefile, ".wind_save");
  strcat (photfile, ".phot");
//  strcat (parfile, ".pf");


  /* Initialize other variables here */

  py_wind_project = 1;		// The default is to try to project onto a yz plane 

/* Read in the wind file */

/* Note that wind_read allocates the space for the WindPtr array.  The
reason that this is done in wind_read is because wind_read also reads
the geo structure.  The geo struc contains the dimensions of the wind
array and so until it is read the space for the WindPtr structure cannot
be allocated.  Howver, one cannot assign w in a subroutine.  Therefore
the old call, namely wind_read(w,windsavefile, will not work.  The statement
w=wmain is a way around this since wmain is external and therefore can
be assigned.  Then w can be set to this ptr and all is restored. The
use of w is endemic in the program. and it is always called through main.
I did not change this now.  Though it could be done.  02apr ksl */

  wind_read (windsavefile);
  w = wmain;

/* aaa is used to store variable for writing to files for the purpose of plotting*/
  aaa = calloc (sizeof (freq), NDIM2);

  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


/* Produce a standard set of output files and exit*/
  if (interactive == 0)
    {
	zoom (1);  /* This affects the logfile */
      ochoice=1;
      complete_file_summary (w, root, ochoice);
      exit (0);
    }

/* Choices */
  ochoice = 0;
  rdint ("Make_files(0=no,1=original,2=regrid_to_linear)", &ochoice);
  c = 'i';
  zoom (1);
  iswitch=0;


a:printf
    ("\nn=ne,  v=vel,  i or I=ion info, j=ave_tau, f=ave_freq, p=nphot\n");
  printf
    ("r=t_r, t=t_e,  w=rad_weight, s=vol,     l=lum,     C=cooling/heating,  b=adiabatic cooling\n");
  printf
    ("a=abs, c=c4,   g=photo,      h=recomb,  k=tau H,   l=lum,     m=F_rad,   x=total,y=mod_te,\n");
  printf ("o=overview,    e=everything, P=Partial emission meas\n");
  printf
    ("W=wind_regions, D=dvds_ave, X=position summary, M=macro atom info, G=inner shell\n");
  printf
    ("z=Zoom,u=unZoom,Z=switch to/from raw and yz projected modes, F=Create files, A=Change file write defaults\n");

  printf ("EOF=quit\n");
  printf ("Model %s   :\n", root);
  rdchar ("Choice", &c);
  switch (c)
    {
    case 'a':			/* Energy absorbed */
      abs_summary (w, root, ochoice);
      break;
    case 'A':			// Change the file defaults
      rdint ("Make_files(0=no,1=original,2=regrid_to_linear)", &ochoice);
      break;

      break;
    case 'b':			/*Adiabatic cooling */
      adiabatic_cooling_summary (w, root, ochoice);
      break;
    case 'c':			/*C4 emission */
      line_summary (w, n, istate, root, ochoice);
      break;
    case 'C':			/*the ratio cooling to heating */
      coolheat_summary (w, root, ochoice);
      break;
    case 'D':			/* dvds summary */
      dvds_summary (w, root, ochoice);
      break;
    case 'e':			/* print out everything about an element */
      wind_element (w);
      break;
    case 'f':			/* Electron summary */
      freq_summary (w, root, ochoice);
      break;
    case 'F':			/* Complete file summary */
      complete_file_summary (w, root, ochoice);
      break;
    case 'g':			/*n photo */
      photo_summary (w, root, ochoice);
      break;
    case 'G':			/* inner shell summary */
      inner_shell_summary (w, root, ochoice);
      break;
    case 'h':			/*n photo */
      recomb_summary (w, root, ochoice);
      break;
    case 'i':			/* Allow user to display information about ions in the wind */
    case 'I':			/* Allow user to display information about ions in the wind */

      rdint("Ion_info_type(0=fraction,1=density,2=scatters,3=abs",&iswitch);

      n = 6;
      istate = 4;

      while (rdint ("element(0=return)", &n) != EOF)
	{
	  if (n <= 0)
	    break;
	  rdint ("ion", &istate);
	  ion_summary (w, n, istate, iswitch, root, ochoice);	// 0 implies ion fractions
	}
      break;

    case 'j':			/* Calculate the average tau at the center of a cell */
      n = 6;
      istate = 4;
      lambda = 1550;

      rddoub ("wavelength", &lambda);
      freq = C / (lambda * 1.e-8);

      while (rdint ("element(0=return)", &n) != EOF)
	{
	  if (n <= 0)
	    break;
	  rdint ("ion", &istate);
	  tau_ave_summary (w, n, istate, freq, root, ochoice);
	}
      break;
    case 'k':			/* tau at H edge */
      tau_h_summary (w, root, ochoice);
      break;
    case 'l':			/* Lum of shell */
      lum_summary (w, root, ochoice);
      break;
    case 'm':			/* Radiation force */
      mo_summary (w, root, ochoice);
      break;
    case 'M':
      macro_summary (w, root, ochoice);
      break;
    case 'n':			/* Electron summary */
      electron_summary (w, root, ochoice);
      break;
    case 'o':			/* overview */
      overview (w, root);
      break;
    case 'p':			/* nphot summary */
      nphot_summary (w, root, ochoice);
      break;
    case 'P':			/* Allow user to display information about the wind */

      n = 6;
      istate = 4;

      while (rdint ("element(0=return)", &n) != EOF)
	{
	  if (n <= 0)
	    break;
	  rdint ("ion", &istate);
	  partial_measure_summary (w, n, istate, root, ochoice);
	}
      break;
    case 'r':			/* Temp summary */
      temp_rad (w, root, ochoice);
      break;
    case 's':			/* Volume summary */
      vol_summary (w, root, ochoice);
      break;
    case 't':			/* Temp summary */
      temp_summary (w, root, ochoice);
      break;
    case 'v':			/* Velocity summary */
      velocity_summary (w, root, ochoice);
      break;
    case 'w':			/* inten weight summary */
      weight_summary (w, root, ochoice);
      break;
    case 'W':			/*Show regions in the wind */
      wind_reg_summary (w, root, ochoice);
      break;
    case 'x':			/*Total emission */
      total_emission_summary (w, root, ochoice);
      break;
    case 'X':			/* Position summary */
      position_summary (w);
      break;
    case 'y':			/* Recalculate temperatures */
      modify_te (w, root, ochoice);
      break;
    case 'z':			/* inspect a specific region */
      zoom (0);
      break;
    case 'Z':			/* Switch between raw and projected mode */
      if (py_wind_project == 0)
	{
	  py_wind_project = 1;
	  Log ("Switching to raw display, whatever the coordinate system");
	}
      else if (py_wind_project == 1)
	{
	  py_wind_project = 0;
	  Log ("Switching to projected y z display");
	}
      break;
    case 'u':			/* Go back to full image */
      zoom (1);
      break;
    case 'q':			/* quit */
      exit (0);
      break;

    }

  goto a;

}



int
py_wind_help ()
{

  char *some_help;

/* Beginning of a string to provide some help for py_wind */
  some_help = "\
\n\
This program reads a wind save file created by python and examine the wind structure\\
\n\
	Usage: py_wind [-h] [-s] [root] \n\
\n\
	where\n\
		-h 	prints out a short help file and exits (see help routine below) \n\
		-s	causes certain parameters to be printed out in individual files \n\
			after which the program exits \n\
		root	optional root name of wind_save file.  If this is not given, \n\
			the user is queried for this \n\
\n\
";

  printf ("%s\n", some_help);

  exit (0);
}
