


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

  This is roughly what the Makefile should look like. Uncomment the mv 
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

simple: simple.o
	gcc simple.o  $(LDFLAGS) -o simple
#	mv $@ $(BIN)

 

  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"

#define LINELENGTH 132

int
main (argc, argv)
     int argc;
     char *argv[];
{
  int i;
  int restart_stat;
  double tmax;
  char root[132];
  char dummy[132];
  char diagfile[LINELENGTH];



  strcpy (root, "foo");
  restart_stat = 0;
  tmax = -1;


  /* The last command line variable is always the .pf file */

  strcpy (dummy, argv[argc - 1]);
  get_root (root, dummy);

  strcpy (diagfile, root);
  strcat (diagfile, ".diag");

  for (i = 1; i < argc - 1; i++)
    {
      if (strcmp (argv[i], "-r") == 0)
	{
	  printf ("Restarting %s\n", root);
	  restart_stat = 1;
	}
      else if (strcmp (argv[i], "-t") == 0)
	{
	  if (sscanf (argv[i + 1], "%lf", &tmax) != 1)
	    {
	      Error ("python: Expected time after -t switch\n");
	      exit (0);
	    }
	  i++;

	}
    }



  if (restart_stat)
    {				// Then we are restarting
      xsignal (root, "%-10s %d %d\n", "RESTART", 1, -15);
      Log_append (diagfile);
    }
  else
    {				// Then we are simply running from a new model
      xsignal (root, "%-10s %d %d\n", "START", 1, -15);
      Log_init (diagfile);
    }

  if (tmax > 0)
    {
      set_max_time (root, tmax);
    }


  xsignal (root, "%-10s %s\n", "NOK", "Working on something");
  sleep (5);			// Sleep 10 seconds
  xsignal (root, "%-10s %s\n", "OK", "Saved everything");
  xsignal (root, "%-10s  %d %d\n", "CHECKTIME", 2, -2);
  check_time (root);

  xsignal (root, "%-10s %s\n", "NOK", "Working on something else");
  sleep (5);			// Sleep 10 seconds
  check_time (root);

  xsignal (root, "%-10s  %d %d\n", "COMPLETE", 1, -15);


}
