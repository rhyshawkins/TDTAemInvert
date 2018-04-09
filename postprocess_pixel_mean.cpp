//
//    AEM Invert : Software for inversion of AEM data using the
//    trans-dimensional tree method and forward modelling code
//    written by Ross Brodie from Geoscience Australia. See
//
//      R Hawkins, R Brodie and M Sambridge, "Bayesian trans-dimensional inversion of
//    Airborne Electromagnetic 2D Conductivity profiles", Exploration Geophysics, 2017
//    https://doi.org/10.1071/EG16139
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

#include "global.hpp"
#include "chainhistory_pixel.hpp"

static const double CREDIBLE_INTERVAL = 0.95;

static int histogram_index(double v, double vmin, double vmax, int bins);

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins);
static double median_from_histogram(int *hist, double vmin, double vmax, int bins);
static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);
static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);

static char short_options[] = "d:l:i:o:v:D:t:s:m:M:c:C:g:b:z:Z:Lh";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"variance", required_argument, 0, 'v'},
  {"stddev", required_argument, 0, 'D'},
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"mode", required_argument, 0, 'm'},
  {"median", required_argument, 0, 'M'},
  {"credible-min", required_argument, 0, 'c'},
  {"credible-max", required_argument, 0, 'C'},
  {"histogram", required_argument, 0, 'g'},

  {"bins", required_argument, 0, 'b'},
  {"vmin", required_argument, 0, 'z'},
  {"vmax", required_argument, 0, 'Z'},

  {"log", no_argument, 0, 'L'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  char *input_file;
  char *output_file;
  char *variance_file;
  char *stddev_file;

  int degree_depth;
  int degree_lateral;

  int thin;
  int skip;

  char *mode_file;
  char *median_file;
  char *credible_min;
  char *credible_max;
  char *histogram;

  int bins;
  double vmin;
  double vmax;

  FILE *fp_out;

  int credible_drop;
  bool logimage;

  /*
   * Default values
   */
  fp_out = NULL;
  degree_depth = 5;
  degree_lateral = 8;
  
  input_file = NULL;
  output_file = NULL;
  variance_file = NULL;
  stddev_file = NULL;
  
  mode_file = NULL;
  median_file = NULL;
  credible_min = NULL;
  credible_max = NULL;
  histogram = NULL;
  
  bins = 1000;
  vmin = 2.0;
  vmax = 4.0;

  logimage = false;
  vmin = 0.001;
  vmax = 1.0;
  
  thin = 0;
  skip = 0;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {
    case 'd':
      degree_depth = atoi(optarg);
      if (degree_depth < 1) {
	fprintf(stderr, "error: invalid degree\n");
	return -1;
      }
      break;

    case 'l':
      degree_lateral = atoi(optarg);
      if (degree_lateral < 1) {
	fprintf(stderr, "error: invalid lateral degree\n");
	return -1;
      }
      break;

    case 'i':
      input_file = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'v':
      variance_file = optarg;
      break;

    case 'D':
      stddev_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 's':
      skip = atoi(optarg);
      break;
      
    case 'm':
      mode_file = optarg;
      break;
      
    case 'M':
      median_file = optarg;
      break;

    case 'c':
      credible_min = optarg;
      break;

    case 'C':
      credible_max = optarg;
      break;

    case 'g':
      histogram = optarg;
      break;

    case 'b':
      bins = atoi(optarg);
      if (bins < 1) {
	fprintf(stderr, "error: bins must be 1 or greater\n");
	return -1;
      }
      break;

    case 'z':
      vmin = atof(optarg);
      break;

    case 'Z':
      vmax = atof(optarg);
      break;

    case 'L':
      logimage = true;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
      
    }
  }

  if (input_file == NULL) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  if (output_file == NULL) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  ChainHistoryPixel *ch = ChainHistoryPixel::load(input_file);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }

  int width = ch->columns;
  int height = ch->rows;
  int size = width * height;

  double *image = new double[size];

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {
      if (logimage) {
	image[j * width + i] = exp(ch->initial_image[j * width + i]);
      } else {
	image[j * width + i] = ch->initial_image[j * width + i];
      }
    }
  }

  double *mean = new double[width * height];
  memset(mean, 0, sizeof(double) * width * height);

  double *variance = new double[width * height];
  memset(variance, 0, sizeof(double) * width * height);

  int *hist = new int[width * height * bins];
  memset(hist, 0, sizeof(int) * width * height * bins);

  int i = 0;
  int counter = 0;

  for (auto &pp : ch->history) {

    if (pp.accepted) {
      if (logimage) {
	image[pp.idx] = exp(pp.newvalue);
      } else {
	image[pp.idx] = pp.newvalue;
      }	
    }

    if (i > skip) {

      if (thin <= 1 || (i - skip) % thin == 0) {
	/*
	 * Update mean/variance calculation
	 */
	counter ++;
	for (int j = 0; j < size; j ++) {
	  double delta = image[j] - mean[j];
	  
	  mean[j] += delta/(double)(counter);
	  variance[j] += delta * (image[j] - mean[j]);
	}
	
	/*
	 * Update the histogram
	 */
	for (int j = 0; j < size; j ++) {
	  int hi = histogram_index(image[j], vmin, vmax, bins);
	  
	  hist[j * bins + hi] ++;
	}
      }
    }

    i ++;
  }
    
  /*
   * Mean output
   */
  fp_out = fopen(output_file, "w");
  if (fp_out == NULL) {
    fprintf(stderr, "error: failed to open mean file\n");
    return -1;
  }

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {
      fprintf(fp_out, "%10.6f ", mean[j*width + i]);
    }
    fprintf(fp_out, "\n");
  }

  fclose(fp_out);

  for (i = 0; i < size; i ++) {
    variance[i] /= (double)(counter - 1);
  }

  /*
   * Variance output
   */
  if (variance_file != NULL) {
    fp_out = fopen(variance_file, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open variance file\n");
      return -1;
    }
    
    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	fprintf(fp_out, "%10.6f ", variance[j*width + i]);
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }

  /*
   * Std. Deviation output
   */
  if (stddev_file != NULL) {
    fp_out = fopen(stddev_file, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open stddev file\n");
      return -1;
    }
    
    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	fprintf(fp_out, "%10.6f ", sqrt(variance[j*width + i]));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }

  /*
   * Mode output
   */
  if (mode_file != NULL) {
    fp_out = fopen(mode_file, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open mode file\n");
      return -1;
    }

    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	fprintf(fp_out, "%10.6f ", mode_from_histogram(hist + (j*width + i) * bins,
						       vmin, vmax, bins));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }
  
  /*
   * Median output
   */
  if (median_file != NULL) {
    fp_out = fopen(median_file, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open median file\n");
      return -1;
    }

    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	fprintf(fp_out, "%10.6f ", median_from_histogram(hist + (j*width + i) * bins,
							 vmin, vmax, bins));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }

  /*
   * Credible Min
   */
  credible_drop = (int)(((double)counter * (1.0 - CREDIBLE_INTERVAL))/2.0);
  
  if (credible_min != NULL) {
    fp_out = fopen(credible_min, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open credible min file\n");
      return -1;
    }

    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	fprintf(fp_out, "%10.6f ", head_from_histogram(hist + (j*width + i) * bins,
						       vmin, vmax, bins,
						       credible_drop));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }
  
  /*
   * Credible Max
   */
  if (credible_max != NULL) {
    fp_out = fopen(credible_max, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open credible max file\n");
      return -1;
    }

    for (int j = 0; j < height; j ++) {
      for (int i = 0; i < width; i ++) {
	fprintf(fp_out, "%10.6f ", tail_from_histogram(hist + (j*width + i) * bins,
						       vmin, vmax, bins,
						       credible_drop));
      }
      fprintf(fp_out, "\n");
    }

    fclose(fp_out);
  }

  if (histogram != NULL) {
    fp_out = fopen(histogram, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to open histogram file\n");
      return -1;
    }

    fprintf(fp_out, "%d %d\n", size, bins);
    fprintf(fp_out, "%.6f %.6f\n", vmin, vmax);

    for (int j = 0; j < size; j ++) {
      for (int i = 0; i < bins; i ++) {

	fprintf(fp_out, "%d ", hist[j * bins + i]);

      }
      fprintf(fp_out, "\n");
    }
    
    fclose(fp_out);
  }

  
  return 0;
}


static int histogram_index(double v, double vmin, double vmax, int bins)
{
  int i;
  
  i = (int)((double)bins * (v - vmin)/(vmax - vmin));

  if (i < 0) {
    return 0;
  }

  if (i > (bins - 1)) {
    return bins - 1;
  }

  return i;
}

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins)
{
  int i;
  int m;
  int mi;

  m = 0;
  mi = -1;

  for (i = 0; i < bins; i ++) {
    if (hist[i] > m) {
      m = hist[i];
      mi = i;
    }
  }
  
  if (mi < 0) {
    return 0.0;
  }

  return ((double)mi + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double median_from_histogram(int *hist, double vmin, double vmax, int bins)
{
  int i;
  int j;
  int ci;
  int cj;

  i = 0;
  j = bins - 1;
  ci = 0;
  cj = 0;

  while (i != j) {
    if (ci < cj) {
      ci += hist[i];
      i ++;
    } else {
      cj += hist[j];
      j --;
    }
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
  int i;
  int ci;

  i = 0; 
  ci = 0;
  while(i < bins && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i ++;
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
  int i;
  int ci;

  i = bins - 1; 
  ci = 0;
  while(i > 0 && ci < drop) {
    if (hist[i] + ci >= drop) {
      break;
    }

    ci += hist[i];
    i --;
  }

  return ((double)i + 0.5)/(double)bins * (vmax - vmin) + vmin;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  " -o|--output <file>               Output mean model file\n"
	  " -v|--variance <file>             Output variance model file\n"
	  " -D|--stddev <file>               Output std dev. file\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -m|--mode <file>                 Output mode\n"
	  " -M|--median <file>               Output median\n"
	  " -c|--credible-min <file>         Output credible min file\n"
	  " -C|--credible-max <file>         Output credible max file\n"
	  " -g|--histogram <file>            Output histogram file\n"
	  "\n"
	  " -b|--bins <int>                  No. histogram bins\n"
	  " -z|--vmin <float>                Lower range for histogram\n"
	  " -Z|--vmax <float>                Upper range for histogram\n"
	  "\n"
	  " -L|--log                         Chain is in log space (exp to invert)\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
