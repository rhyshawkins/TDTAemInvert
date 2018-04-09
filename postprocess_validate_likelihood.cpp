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

extern "C" {
#include "chain_history.h"
  
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
#include "generic_lift.h"
  
};

#include "global.hpp"

struct user_data {
  int stepcounter;
  int thincounter;
  int thin;
  int skip;
  int max;
  
  int counter;
  double maxerror;
  
  Global *global;
};

static const double CREDIBLE_INTERVAL = 0.95;

static int process(int i,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v);

static char short_options[] = "O:S:H:d:l:D:i:t:s:m:w:W:h";
static struct option long_options[] = {
  {"observations", required_argument, 0, 'O'},
  {"stm", required_argument, 0, 'S'},
  {"hierarchical", required_argument, 0, 'H'},
  
  {"degree-depth", required_argument, 0, 'd'},
  {"degree-lateral", required_argument, 0, 'l'},
  {"depth", required_argument, 0, 'D'},
  
  {"input", required_argument, 0, 'i'},

  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},
  {"max", required_argument, 0, 'm'},
  
  {"wavelet-vertical", required_argument, 0, 'w'},
  {"wavelet-horizontal", required_argument, 0, 'W'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  chain_history_t *ch;

  char *observations;
  std::vector<std::string> stm_files;
  std::vector<std::string> hierarchical_files;
  
  char *input_file;

  int degree_depth;
  int degree_lateral;
  double depth;
  
  int thin;
  int skip;
  int maxcheck;
  int maxsteps;

  FILE *fp_in;

  struct user_data data;
  multiset_int_double_t *S_v;

  int waveleth;
  int waveletv;

  /*
   * Default values
   */
  observations = nullptr;
  
  fp_in = NULL;

  degree_depth = 5;
  degree_lateral = 8;
  depth = 200.0;
  
  input_file = NULL;
  
  thin = 0;
  skip = 0;
  maxcheck = 1000;
  
  maxsteps = 1000000;

  waveleth = 0;
  waveletv = 0;

  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {

    case 'O':
      observations = optarg;
      break;

    case 'S':
      stm_files.push_back(optarg);
      break;

    case 'H':
      hierarchical_files.push_back(optarg);
      break;
      
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

    case 'D':
      depth = atof(optarg);
      if (depth <= 0.0) {
	fprintf(stderr, "error: depth must be greater than 0\n");
	return -1;
      }
      break;

    case 'i':
      input_file = optarg;
      break;

    case 't':
      thin = atoi(optarg);
      break;

    case 's':
      skip = atoi(optarg);
      break;

    case 'm':
      maxcheck = atoi(optarg);
      if (maxcheck <= 0) {
	fprintf(stderr, "error: maxcheck must be greater than 0\n");
	return -1;
      }
      break;
      
    case 'w':
      waveletv = atoi(optarg);
      if (waveletv < 0 || waveletv > Global::WAVELET_MAX) {
	fprintf(stderr, "error: vertical wavelet must be between 0 and %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'W':
      waveleth = atoi(optarg);
      if (waveleth < 0 || waveleth > Global::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be between 0 and %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
      
    }
  }

  if (observations == NULL) {
    fprintf(stderr, "error: required parameter observations file missing\n");
    return -1;
  }

  if (input_file == NULL) {
    fprintf(stderr, "error: required parameter input file missing\n");
    return -1;
  }

  ch = chain_history_create(maxsteps);
  if (ch == NULL) {
    fprintf(stderr, "error: failed to create chain history\n");
    return -1;
  }

  data.stepcounter = 0;
  data.thincounter = 0;
  data.thin = thin;
  data.skip = skip;
  data.max = maxcheck;

  data.counter = 0;
  data.maxerror = 0.0;

  data.global = new Global(observations,
			   stm_files,
			   nullptr,
			   nullptr,
			   degree_lateral,
			   degree_depth,
			   depth,
			   hierarchical_files,
			   0,
			   1000,
			   false,
			   waveleth,
			   waveletv);

  S_v = multiset_int_double_create();
  if (S_v == NULL) {
    fprintf(stderr, "error: failed to create multiset\n");
    return -1;
  }
  
  fp_in = fopen(input_file, "r");
  if (fp_in == NULL) {
    fprintf(stderr, "error: failed to open input file\n");
    return -1;
  }

  /*
   * Process the chain history
   */
  while (!feof(fp_in)) {

    if (chain_history_read(ch,
			   (ch_read_t)fread,
			   fp_in) < 0) {
      if (feof(fp_in)) {
	break;
      }
      
      fprintf(stderr, "error: failed to read chain history\n");
      return -1;
    }

    if (chain_history_replay(ch,
			     S_v,
			     (chain_history_replay_function_t)process,
			     &data) < 0) {
      fprintf(stderr, "error: failed to replay\n");
      return -1;
    }
  }
  printf("Checked %d/%d(%d) records\n", data.counter, data.thincounter, data.stepcounter);
  printf("Max. Error: %.6g\n", data.maxerror);
  fclose(fp_in);

  chain_history_destroy(ch);
  multiset_int_double_destroy(S_v);

  delete data.global;

  return 0;
}

static int process(int stepi,
		   void *user,
		   const chain_history_change_t *step,
		   const multiset_int_double_t *S_v)
{
  struct user_data *d = (struct user_data *)user;
  

  if (step->header.accepted) {
    
    if ((d->counter < d->max) &&
	(d->stepcounter >= d->skip) &&
	(d->thin <= 1 || (d->thincounter % d->thin) == 0)) {

      if (wavetree2d_sub_set_from_S_v(d->global->wt, S_v) < 0) {
	fprintf(stderr, "process: failed to set wavetree (sub)\n");
	return -1;
      }
      
      d->global->lambda_scale = step->header.hierarchical;

      double log_normalization;
      double like = d->global->likelihood(log_normalization);
      double error = fabs(step->header.likelihood - like); 
      printf("Step %d, %d: %d %10.6f Stored %10.6f Computed (delta %.6g)\n",
	     d->stepcounter,
	     d->thincounter,
	     (int)step->header.type,
	     step->header.likelihood,
	     like,
	     error);
      
      if (error > d->maxerror) {
	d->maxerror = error;
      }
      
      d->counter ++;
      
    }
    d->thincounter ++;
  }
  
  d->stepcounter ++;
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -O|--observations <filename>     Input observations\n"
	  " -S|--stm <filename>              Input stm file\n"
	  " -H|--hierarchical <filename>     Input hierarchical file\n"
	  "\n"
	  " -d|--degree-depth <int>    Number of layers as power of 2\n"
	  " -l|--degree-lateral <int>  Number of horizontal samples as power of 2\n"
	  "\n"
	  " -i|--input <file>                Input ch file\n"
	  "\n"
	  " -t|--thin <int>                  Only processing every ith sample\n"
	  " -s|--skip <int>                  Skip n samples from beginning\n"
	  "\n"
	  " -w|--wavelet-vertical <int>      Wavelet for vertical direction\n"
	  " -W|--wavelet-horizontal <int>    Wavelet for horizontal direction\n"
	  "\n"
	  " -h|--help            Show usage\n"
	  "\n",
	  pname);
}

 
