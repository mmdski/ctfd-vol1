/**
 * @file prob3.7.c
 * @author Marian Domanski (mmdski@gmail.com)
 * @brief Problem 3.7
 * @version 0.1
 * @date 2022-05-14
 *
 * @details Ported from fortran in Fletcher, 1998 (Fig. 1.13)
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.1415927

/**
 * @brief Solves 1D transient heat conduction equation using FTCS scheme
 *
 * @param jmax the number of points along the rod
 * @param maxex the number of terms in the exact solution
 * @param nmax the maximum number of time steps
 * @param alph the thermal diffusivity
 * @param s alph*delt/delx/delx
 * @param tmax the minimum time
 * @return int
 */
int
diff (FILE  *file,
      int    jmax,
      int    maxex,
      int    nmax,
      double alph,
      double s,
      double tmax)
{
  // td - dimensional temperature
  // tn -w nondimensional temperature

  double tn[41];
  double dum[41];
  double td[41];
  double x[41];
  double te[41];

  size_t jmap = jmax - 1;

  int    ajm  = jmap;
  double delx = 1. / ajm;
  double delt = delx * delx * s / alph;

  fprintf (file,
           " JMAX=%5d  MAXEX=%5d  NMAX=%5d  TMAX=%8.2f\n",
           jmax,
           maxex,
           nmax,
           tmax);
  fprintf (file,
           " S=%5.3f  ALPH =%10.3e  DELT =%10.3e  DELX =%10.3e\n\n",
           s,
           alph,
           delt,
           delx);
  fprintf (file, " FTCS(EXPLICIT) SCHEME S = %5.3f\n\n", s);

  // set initial conditions
  for (size_t j = 0; j < jmap; j++)
    {
      tn[j] = 0;
    }
  double n  = 0;
  double t  = 0;
  double sj = 1.0 - 2.0 * s;

  size_t j;

  // set boundary conditions
  while (t < tmax)
    {
      if (t < 0.01)
        {
          tn[0]        = 0.5;
          tn[jmax - 1] = 0.5;
        }
      else
        {
          tn[0]        = 1;
          tn[jmax - 1] = 1;
        }
      td[0]        = 100 * tn[0];
      td[jmax - 1] = 100 * tn[jmax - 1];

      // compute f.d. solution
      j      = 1;
      dum[j] = 11. / 12. * s * tn[j - 1] + (1. - 5. / 3. * s) * tn[j] +
               0.5 * s * tn[j + 1] + 1. / 3. * s * tn[j + 2] -
               1. / 12. * s * tn[j + 3];

      for (j = 2; j < jmap - 1; j++)
        dum[j] = -1. / 12. * s * tn[j - 2] + 4. / 3. * s * tn[j - 1] +
                 (1. - 2.5 * s) * tn[j] + 4. / 3. * s * tn[j + 1] -
                 1. / 12. * s * tn[j + 2];

      j      = jmap - 1;
      dum[j] = -1. / 12. * s * tn[j - 3] + 1. / 3. * s * tn[j - 2] +
               0.5 * s * tn[j - 1] + (1. - 5. / 3. * s) * tn[j] +
               11. / 12. * s * tn[j + 1];

      for (size_t j = 1; j < jmap; j++)
        tn[j] = dum[j];

      for (size_t j = 1; j < jmap; j++)
        td[j] = 100 * tn[j];

      t = t + delt;
      fprintf (file, " T= %5.0f  TD=", t);
      for (size_t j = 0; j < jmax; j++)
        fprintf (file, "%6.2f", td[j]);
      fprintf (file, "\n");
    }

  // obtain exact solution and compare
  double sum = 0;
  double dam;
  double dxm;
  double dtm;
  double aj;
  double am;
  for (size_t j = 0; j < jmax; j++)
    {
      aj    = (double) j;
      x[j]  = delx * aj;
      te[j] = 100;

      for (size_t m = 1; m <= maxex; m++)
        {
          am  = (double) m;
          dam = (2. * am - 1.);
          dxm = dam * PI * x[j];
          dtm = -alph * dam * dam * PI * PI * t;

          // limit the argument size of exp(dtm)
          if (dtm < -87.)
            dtm = -87.;

          te[j] = te[j] - 400. / dam / PI * sin (dxm) * exp (dtm);
        }
      sum = sum + pow ((te[j] - td[j]), 2);
    }

  fprintf (file, "\n T= %5.0f  TE=", t);
  for (size_t j = 0; j < jmax; j++)
    fprintf (file, "%6.2f", te[j]);
  fprintf (file, "\n\n");

  // rms is the rms error
  double avs = sum / (1. + ajm);
  double rms = sqrt (avs);
  fprintf (file, "  RMS DIF = %11.4e\n\n", rms);

  return 0;
}

int
main (void)
{
  int    jmax  = 0;
  int    maxex = 0;
  int    nmax  = 0;
  double alph  = 0;
  double s     = 0;
  double tmax  = 0;

  int status = EXIT_SUCCESS;

  const char *path = "prob3.7.dat";
  FILE       *dat;
  if ((dat = fopen (path, "r")) == NULL)
    {
      fprintf (stderr, "Unable to open file: %s\n", path);
      return EXIT_FAILURE;
    }

  do
    {
      int n = fscanf (dat,
                      "%5d%5d%5d%10le%5lf%5lf\n",
                      &jmax,
                      &maxex,
                      &nmax,
                      &alph,
                      &s,
                      &tmax);
      if (n == 6)
        diff (stdout, jmax, maxex, nmax, alph, s, tmax);
      else if (n != EOF)
        {
          fputs ("Failed to match input\n", stderr);
          status = EXIT_FAILURE;
          break;
        }
      else
        {
          break;
        }
    }
  while (1);

  if (fclose (dat) == EOF)
    {
      fputs ("Failed to close file\n", stderr);
      status = EXIT_FAILURE;
    }

  return status;
}
