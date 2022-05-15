#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int
main (void)
{

  FILE *gnuplot = popen ("gnuplot -persistent", "w");

  fprintf (gnuplot, "set title \"Simple Plots\" font \",20\"\n");
  fprintf (gnuplot, "set key left bottom box opaque font \",20\"\n");
  fprintf (gnuplot, "set xrange [ 0 : 2*pi ] noreverse writeback\n");

  fprintf (gnuplot,
           "plot "
           "'-' title \"sin(x)\" with lines, "
           "'-' title \"cos(x)\" with lines\n");

  // plot sine
  for (double x = 0; x <= 2 * M_PI; x += 2 * M_PI / 100)
    fprintf (gnuplot, "%g %g\n", x, sin (x));
  fprintf (gnuplot, "e\n");

  // plot cosine
  for (double x = 0; x <= 2 * M_PI; x += 2 * M_PI / 100)
    fprintf (gnuplot, "%g %g\n", x, cos (x));
  fprintf (gnuplot, "e\n");

  fflush (gnuplot);

  return EXIT_SUCCESS;
}
