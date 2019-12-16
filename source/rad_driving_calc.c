
/***********************************************************/
/** @file  windsave2table.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  A standalone routine for writing a standard set of data
 * from a windsavefile to ascii tables all readable with as astropy.io.ascii
 *
 * This routine is run from the command line, as follows
 *
 * windsave2table  windsave_root
 *
 * where windsave_root is the root name for a python run, or more precisely
 * the rootname of a windsave file, as the .pf file is not read.
 *
 * The routine reads the windsavefile and then writes out a selected 
 * set of variables into a variety of number of files, each of which
 * are readable as astropy tables, where each row corresponds to one
 * element of the wind.   
 *
 * All of the files begin
 * with the rootname and the domain number.
 * If Python run that created the windsave file
 * has multiple domains a sepearate output file is created for each
 * domain
 *
 * All of the files, begin with
 *
 *        x        z    i    j inwind 
 *
 * where x and z correspond to the position of the cell, i and j correspond to 
 * the cell number and inwind says whether the cell was in the wind, so that it is
 * fairly straightforward to create plotting routines for various parameters
 * that are contained in the remaining columns.  
 *
 * The files include a so-called masterfile which records basic parameters
 * like n_e, velocity,  rho, t_e, t_r, and ionization fractions
 * for a few key ions in each cell when the windsave file was written
 *
 * ### Notes ###
 *
 * Whereas py_wind is intended to be run interactively, windsave2table is
 * entirely hardwired so that it produces a standard set of output
 * files.  To change the outputs one has to modify the routine
 *
 * This file just contains the driving routine.  All of the 
 * real work is carried out in windsave2table_sub.c  Indeed so 
 * little of the real work is done here that it might be sensible
 * to move some portion of that code here.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"




/**********************************************************/
/** 
 * @brief      windsave2table writes key variables in a windsave file 
 * to an astropy table calculated Python.  This is the  main routine.
 *
 * @param [in] int  argc   The number of argments in the command line
 * @param [in] char *  argv[]   The command line
 * @return     Always returns 0, unless the windsave file is not
 * found in which case the routine will issue and error before
 * exiting.  
 *
 * @details
 * argc and argv[] are the standard variables provided to main 
 * in a c-program.  Only the first command line argument is 
 * parsed and this should be rootname of the windsave file
 *
 * ### Notes ###
 *
 * This routine is a supervisory routine. The real work 
 * is in do_windsave2table
 *
 * The routine does not read the .pf file.  It does read  
 * the windsave file and the associated atomic data file
 * before calling do_windsave2table.
 *
 **********************************************************/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  char root[LINELENGTH], input[LINELENGTH],wspecfile[LINELENGTH], specfile[LINELENGTH];; 
  char windradfile[LINELENGTH], windsavefile[LINELENGTH];
  char parameter_file[LINELENGTH];
  
  char photfile[LINELENGTH];
  int interactive;
  
  
  int i,j,ii,nwind,nplasma;
  double vol;
  char *fgets_rc;
  
  FILE *fptr, *fptr2, *fptr3, *fptr4, *fptr5, *fopen ();        /*This is the file to communicate with zeus */
  


  interactive = 1;              /* Default to the standard operating mofe for py_wind */
  strcpy (parameter_file, "NONE");

  /* Next command stops Debug statements printing out in py_wind */
  Log_set_verbosity (3);

  /* Parse the command line.  It is Important that py_wind_help
   * below is updated if this changes
   */

  if (argc == 1)
  {
    printf ("Root for wind file :");
    fgets_rc = fgets (input, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      Error ("Input rootname is NULL or EOF\n");
      exit (1);                 // Exit if NULL returned
    }
    get_root (root, input);
  }
  else
  {
    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-h") == 0)
      {
        printf ("No help here buddy\n");
      }
      else if (strcmp (argv[i], "-d") == 0)
      {
        interactive = -1;
      }

      else if (strcmp (argv[i], "-s") == 0)
      {
        interactive = 0;
      }
      else if (strcmp (argv[i], "-p") == 0)
      {
        interactive = 0;
        i = i + 1;
        strcpy (parameter_file, argv[i]);
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        Error ("py_wind: unknown switch %s\n", argv[i]);
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


/* Read in the wind file */

  if (wind_read (windsavefile) < 0)
  {
    Error ("py_wind: Could not open %s", windsavefile);
    exit (0);
  }


  printf ("Read wind_file %s\n", windsavefile);

  get_atomic_data (geo.atomic_filename);

  printf ("Read Atomic data from %s\n", geo.atomic_filename);


  Log ("Outputting heatcool file for connecting to zeus\n");
  fptr = fopen ("py_heatcool.dat", "w");
  fptr2 = fopen ("py_flux.dat", "w");
  fptr3 = fopen ("py_ion_data.dat", "w");
  fptr4 = fopen ("py_spec_data.dat", "w");
  fptr5 = fopen ("py_pcon_data.dat", "w");

  fprintf (fptr,
           "i j rcen thetacen x y z vol xcen ycen zcen temp xi ne v_x v_y v_z rho n_h rad_f_w rad_f_phi rad_f_z bf_f_w bf_f_phi bf_f_z\n");
  fprintf (fptr2, "i j F_vis_x F_vis_y F_vis_z F_UV_x F_UV_y F_UV_z F_Xray_x F_Xray_y F_Xray_z\n");   //directional flux by band

  fprintf (fptr3, "nions %i\n", nions);
  for (i = 0; i < nions; i++)
  {
    fprintf (fptr3, "ion %i %s %i %i\n", i, ele[ion[i].nelem].name, ion[i].z, ion[i].istate);
  }
  fprintf (fptr3, "nplasma %i\n", NPLASMA);

  fprintf (fptr4, "nbands %i\n", geo.nxfreq);
  fprintf (fptr4, "nplasma %i\n", NPLASMA);
  for (i = 0; i < geo.nxfreq + 1; i++)
    fprintf (fptr4, "%e ", geo.xfreq[i]);     //hard wired band edges
  fprintf (fptr4, "\n ");

  fprintf (fptr5, "nplasma %i\n", NPLASMA);
  
  
  printf ("We are going to run from %i to %i\n",zdom[0].nstart,zdom[0].nstop);
  
  for (nwind = zdom[0].nstart; nwind < zdom[0].nstop; nwind++)
  {
    if (wmain[nwind].vol > 0.0)
    {
      nplasma = wmain[nwind].nplasma;
      wind_n_to_ij (0, plasmamain[nplasma].nwind, &i, &j);
      vol = wmain[plasmamain[nplasma].nwind].vol;
      fprintf (fptr, "%d %d %e %e ", i, j, wmain[plasmamain[nplasma].nwind].rcen, wmain[plasmamain[nplasma].nwind].thetacen / RADIAN);        //output geometric things
      fprintf (fptr, "%e %e %e %e ", wmain[plasmamain[nplasma].nwind].x[0],wmain[plasmamain[nplasma].nwind].x[1],wmain[plasmamain[nplasma].nwind].x[2], vol);        //output geometric things
      fprintf (fptr, "%e %e %e ", wmain[plasmamain[nplasma].nwind].xcen[0],wmain[plasmamain[nplasma].nwind].xcen[1],wmain[plasmamain[nplasma].nwind].xcen[2]);        //output geometric things
      
      fprintf (fptr, "%e %e %e ", plasmamain[nplasma].t_e, plasmamain[nplasma].xi, plasmamain[nplasma].ne);   //output temp, xi and ne to ease plotting of heating rates
      fprintf (fptr, "%e %e %e ", wmain[plasmamain[nplasma].nwind].v[0],wmain[plasmamain[nplasma].nwind].v[1],wmain[plasmamain[nplasma].nwind].v[2]); //density

      fprintf (fptr, "%e ", plasmamain[nplasma].rho); //density
      fprintf (fptr, "%e ", plasmamain[nplasma].rho * rho2nh);        //hydrogen number density
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[0]);     //electron scattering radiation force in the w(x) direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[1]);     //electron scattering radiation force in the phi(rotational) directionz direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_es[2]);     //electron scattering radiation force in the z direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_bf[0]);     //bound free scattering radiation force in the w(x) direction
      fprintf (fptr, "%e ", plasmamain[nplasma].rad_force_bf[1]);     //bound free scattering radiation force in the phi(rotational) direction
      fprintf (fptr, "%e \n", plasmamain[nplasma].rad_force_bf[2]);   //bound free scattering radiation force in the z direction
      fprintf (fptr2, "%d %d ", i, j);        //output geometric things               
      fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_vis[0], plasmamain[nplasma].F_vis[1], plasmamain[nplasma].F_vis[2]); //directional flux by band
      fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_UV[0], plasmamain[nplasma].F_UV[1], plasmamain[nplasma].F_UV[2]);    //directional flux by band
      fprintf (fptr2, "%e %e %e ", plasmamain[nplasma].F_Xray[0], plasmamain[nplasma].F_Xray[1], plasmamain[nplasma].F_Xray[2]);      //directional flux by band
      fprintf (fptr2, "\n");
      fprintf (fptr3, "%d %d ", i, j);        //output geometric things               
      for (ii = 0; ii < nions; ii++)
        fprintf (fptr3, "%e ", plasmamain[nplasma].density[ii]);
      fprintf (fptr3, "\n");

      fprintf (fptr4, "%d %d ", i, j);        //output geometric things       
      for (ii = 0; ii < geo.nxfreq; ii++)
        fprintf (fptr4, "%e %e %i %e %e %e %e ",
                 plasmamain[nplasma].fmin_mod[ii], plasmamain[nplasma].fmax_mod[ii], plasmamain[nplasma].spec_mod_type[ii],
                 plasmamain[nplasma].pl_log_w[ii], plasmamain[nplasma].pl_alpha[ii], plasmamain[nplasma].exp_w[ii],
                 plasmamain[nplasma].exp_temp[ii]);
      fprintf (fptr4, "\n ");

      fprintf (fptr5, "%i %i %e %e %e %e %e %e %e\n", i, j, plasmamain[nplasma].t_e, plasmamain[nplasma].rho,
               plasmamain[nplasma].rho * rho2nh, plasmamain[nplasma].ne, plasmamain[nplasma].t_opt, plasmamain[nplasma].t_UV,
               plasmamain[nplasma].t_Xray);
    }
  }



  fclose (fptr);
  fclose (fptr2);
  fclose (fptr3);
  fclose (fptr4);
  fclose (fptr5);
}
  
  
  
  
  
  
 
