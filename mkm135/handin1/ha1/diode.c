#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define q_ele 1.60e-19
#define k_blz 1.38e-23

/**
 * contains data points for voltage (x) and
 * current for diode a,b,c,ideal (y)
 */
typedef struct diode_point_t {
  double x;
  double y[4];
} dpoint;

static void diode_fill_iv(dpoint    *data,
		          size_t    data_points);
static void diode_plot2d(const char     *logscale);
static void diode_save_data(const char  *outfile,
			    dpoint      *data,
			    size_t      len);
static void diode_linspace(double   *a,
			   double   start,
		           double   stop,
		           size_t   length);

int main(int argc, char **argv)
{
  char axistype[10];
  size_t data_points;
  if (argc > 2) {
    if (strncmp(argv[1], "linlog", 10) == 0) {
      strncpy(axistype, "y", 10);
    }
    else if (strncmp(argv[1], "loglin", 10) == 0) {
      strncpy(axistype, "x", 10);
    }
    else if (strncmp(argv[1], "loglog", 10) == 0) {
      strncpy(axistype, "xy", 10);
    }
    else if (strncmp(argv[1], "linlin", 10) == 0) {
      strncpy(axistype, "", 10);
    }
    else {
      fprintf(stdout, "Using default axis scaling (%s)\n", "linlin");
      strncpy(axistype, "", 10);
    }

    data_points = (size_t) strtol(argv[2], &argv[2], 10);
    if (data_points < 2) {
      fprintf(stdout, "Using default number of samples (%d)\n", 101);
      data_points = 101;
    }
  }
  else {
    fprintf(stdout, "Usage: %s <axis-type> <sample-count>\n", argv[0]);
    return EXIT_FAILURE;
  }
  
  dpoint data[data_points];
  diode_fill_iv(data, data_points);

  diode_save_data("iv-data.dat", data, data_points);
  diode_plot2d(axistype);

  return EXIT_SUCCESS;
}

static void diode_fill_iv(dpoint *data, size_t data_points)
{
  double N_epi = 1e21;                      /* epi doping [m^-3] */
  double L_epi = 400e-6;                    /* epi diffusion length [m] */
  double N_sub[] = {1e21, 1e22, 1e23};      /* sub doping [m^-3] */
  double t_p = 200e-9;                      /* p+ region thickness [m] */
  double t_epi = 4e-6;                      /* epi region thickness [m] */
  double t_sub = 250e-6;                    /* sub region thickness [m] */
  double T = 300;                           /* diode temperature [K] */
  double A = 1e-7;                          /* diode cross-section [m^2] */
  double V_T = k_blz*T/q_ele;

  double n_i = 1e16;                        /* intrinsic doping [m^-3] */
  double mu_n[] = {1.5e-1, 1.3e-1, 8e-2};   /* electron mobility */
  double mu_p[] = {5e-2, 4.5e-2, 3e-2};     /* hole mobility */
  
  double R_sub[3];
  for (size_t i = 0; i < 3; ++i) {
    R_sub[i] = t_sub / (A * mu_n[i] * N_sub[i] * q_ele);
  }
  double R_epi = t_epi / (A * mu_n[0] * N_epi * q_ele);
  
  double D_P = V_T * mu_p[0];                   /* N_epi in diode */
  
  /* diode I-V data */
  double V[data_points]; diode_linspace(V, 0, 2, data_points);
  double I[4][data_points];
  for (size_t i = 0; i < 4; ++i) {              /* Diode a, b, c and ideal */
    double vdiv = i != 3 ? R_epi/(R_epi + R_sub[i]) : 1;
    for (size_t j = 0; j < data_points; ++j) {
      I[i][j] = A * q_ele*(n_i*n_i)*(D_P / (L_epi * N_epi))
	* (exp((V[j] * vdiv) / V_T) - 1);
    }
  }

  /* fill I-V data */
  for (size_t i = 0; i < data_points; ++i) {
    data[i].x = V[i];
    for (size_t j = 0; j < 4; ++j) {
      data[i].y[j] = I[j][i];
    }
  }

  /* breakdown voltage, 300 V for ideal diode @ doping 10^21 */
  for (int i = 0; i < 3; ++i) {
    fprintf(stdout, "Breakdown voltage: %5.0lf V\n",
	    300.0 * ((R_epi + R_sub[i])/R_epi));
  }
}

static void diode_linspace(double *a, double start, double stop, size_t length)
{
  double da = (stop - start) / (length - 1);
  double val = start;
  for (unsigned i = 0; i < length; ++i) {
    a[i] = val;
    val += da;
  }
}

static void diode_save_data(const char *outfile, dpoint *data, size_t len)
{
  FILE *plotdata = fopen(outfile, "w");
  for (unsigned i = 0; i < len; i++) {
    fprintf(plotdata, "%g %g %g %g %g\n", data[i].x,
	    data[i].y[0], data[i].y[1], data[i].y[2], data[i].y[3]);
  }
  fclose(plotdata);
}

static void diode_plot2d(const char *logscale)
{
  FILE *gnuplot, *fp;
  if ((gnuplot = popen("gnuplot", "w")) == NULL) {
    fprintf(stderr, "gnuplot could not be invoked with "
	    "\'gnuplot\'\n");
    return;
  }
    
  if (strlen(logscale) > 0) {
    fprintf(gnuplot, "set logscale %s 10\n", logscale);
    fprintf(gnuplot, "set format y \"10^{%%T}\"\n");
  }

  if ((fp = fopen("../plotting.txt", "r")) != NULL) {
    for (int i = 0; (i = fgetc(fp)) != EOF;) {
      fputc(i, gnuplot);
    }
    fclose(fp);
    fflush(gnuplot);
    fgetc(stdin);                 /* keep gnuplot open */
    pclose(gnuplot);
  }
  else {
    fprintf(stderr, "gnuplot settings file not found at "
	    "\'../plotting.txt\'\n");
  }
}
