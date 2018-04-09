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

#include <gmp.h>

extern "C" {
#include "wavetree2d_sub.h"
#include "wavetreepp.h"
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"

#include "slog.h"
};

#include "global.hpp"
#include "birth.hpp"
#include "death.hpp"
#include "value.hpp"
#include "hierarchical.hpp"

#include "aemutil.hpp"

static char short_options[] = "i:I:s:M:o:d:l:D:t:S:H:L:k:B:Pw:W:v:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial", required_argument, 0, 'I'},
  {"stm", required_argument, 0, 's'},
  {"prior-file", required_argument, 0, 'M'},
  {"output", required_argument, 0, 'o'},
  
  {"degree-depth", required_argument, 0, 'd'},
  {"degree-lateral", required_argument, 0, 'l'},

  {"depth", required_argument, 0, 'D'},

  {"total", required_argument, 0, 't'},
  {"seed", required_argument, 0, 'S'},

  {"hierarchical", required_argument, 0, 'H'},
  {"lambda-std", required_argument, 0, 'L'},

  {"kmax", required_argument, 0, 'k'},

  {"birth-probability", required_argument, 0, 'B'},

  {"posteriork", 0, 0, 'P'},

  {"wavelet-vertical", required_argument, 0, 'w'},
  {"wavelet-horizontal", required_argument, 0, 'W'},

  {"verbosity", required_argument, 0, 'v'},

  {"help", no_argument, 0, 'h'},
  
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Parameters
  //
  char *input_obs;
  char *initial_model;
  std::vector<std::string> stm_files;
  std::vector<std::string> hierarchical_files;
  char *prior_file;
  char *output_prefix;

  int degreex;
  int degreey;

  double depth;

  int total;
  int seed;

  double lambda_std;
  int kmax;

  double Pb;

  bool posteriork;

  int wavelet_v;
  int wavelet_h;

  int verbosity;
  
  //
  // Defaults
  //

  input_obs = nullptr;
  initial_model = nullptr;
  prior_file = nullptr;
  output_prefix = nullptr;

  degreex = 10;
  degreey = 5;

  depth = 500.0;

  total = 10000;
  seed = 983;

  lambda_std = 0.0;
  kmax = 100;
  
  Pb = 0.05;

  posteriork = false;

  wavelet_v = 0;
  wavelet_h = 0;

  verbosity = 1000;

  //
  // Command line parameters
  //
  option_index = 0;
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input_obs = optarg;
      break;

    case 'I':
      initial_model = optarg;
      break;

    case 's':
      stm_files.push_back(optarg);
      break;

    case 'M':
      prior_file = optarg;
      break;

    case 'o':
      output_prefix = optarg;
      break;

    case 'd':
      degreey = atoi(optarg);
      if (degreey < 1 || degreey > 16) {
	fprintf(stderr, "error: degree y must be between 1 and 16 inclusive\n");
	return -1;
      }
      break;

    case 'l':
      degreex = atoi(optarg);
      if (degreex < 1 || degreex > 16) {
	fprintf(stderr, "error: degree x must be between 1 and 16 inclusive\n");
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
      
    case 't':
      total = atoi(optarg);
      if (total < 1) {
	fprintf(stderr, "error: total must be greater than 0\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;

    case 'H':
      hierarchical_files.push_back(optarg);
      break;

    case 'L':
      lambda_std = atof(optarg);
      if (lambda_std <= 0.0) {
	fprintf(stderr, "error: lambda std dev must be greater than 0\n");
	return -1;
      }
      break;
      
    case 'k':
      kmax = atoi(optarg);
      if (kmax < 1) {
	fprintf(stderr, "error: kmax must be greater than 0\n");
	return -1;
      }
      break;

    case 'B':
      Pb = atof(optarg);
      if (Pb < 0.0 || Pb > 0.5) {
	fprintf(stderr, "error: birth probability must be between 0 and 0.5\n");
	return -1;
      }
      break;

    case 'P':
      posteriork = true;
      break;

    case 'w':
      wavelet_v = atoi(optarg);
      if (wavelet_v < 0 || wavelet_v > Global::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'W':
      wavelet_h = atoi(optarg);
      if (wavelet_h < 0 || wavelet_h > Global::WAVELET_MAX) {
	fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)Global::WAVELET_MAX);
	return -1;
      }
      break;

    case 'v':
      verbosity = atoi(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (input_obs == nullptr) {
    fprintf(stderr, "error: required input parameter input observations missing\n");
    return -1;
  }

  if (stm_files.size() == 0) {
    fprintf(stderr, "error: need at least on stm file\n");
    return -1;
  }

  if (prior_file == nullptr) {
    fprintf(stderr, "error: required prior file parameter missing\n");
    return -1;
  }

  if (stm_files.size() != hierarchical_files.size()) {
    fprintf(stderr, "error: no. of stm and hierarchical files must match %d stm %d hierarchical \n",
	    (int)stm_files.size(),
	    (int)hierarchical_files.size());
    return -1;
  }

  Global global(input_obs,
		stm_files,
		initial_model,
		prior_file,
		degreex,
		degreey,
		depth,
		hierarchical_files,
		seed,
		kmax,
		posteriork,
		wavelet_h,
		wavelet_v);

  Birth birth(global);
  Death death(global);
  Value value(global);

  Hierarchical *hierarchical = nullptr;
  if (lambda_std > 0.0) {
    hierarchical = new Hierarchical(global, lambda_std);
  }

  global.current_likelihood = global.likelihood(global.current_log_normalization);

  printf("Initial Likelihood: %f\n", global.current_likelihood);

  int *khistogram = new int[kmax];
  for (int i = 0; i < kmax; i ++) {
    khistogram[i] = 0;
  }

  FILE *fp_ch = NULL;
  if (chain_history_initialise(global.ch,
			       wavetree2d_sub_get_S_v(global.wt),
			       global.current_likelihood,
			       global.temperature,
			       global.lambda_scale) < 0) {
    fprintf(stderr, "error: failed to initialise chain history\n");
    return -1;
  }
  
  std::string filename = mkfilename(output_prefix, "ch.dat");
  fp_ch = fopen(filename.c_str(), "w");
  if (fp_ch == NULL) {
    fprintf(stderr, "error: failed to create chain history file\n");
    return -1;
  }

  for (int i = 0; i < total; i ++) {

    double u = global.random.uniform();

    if (u < Pb) {

      //
      // Birth
      //
      if (birth.step() < 0) {
	fprintf(stderr, "error: failed to do birth step\n");
	return -1;
      }

    } else if (u < (2.0 * Pb)) {

      //
      // Death
      //
      if (death.step() < 0) {
	fprintf(stderr, "error: failed to do death step\n");
	return -1;
      }

    } else {

      //
      // Value
      //
      if (value.step() < 0) {
	fprintf(stderr, "error: failed to do value step\n");
	return -1;
      }

    }

    if (chain_history_full(global.ch)) {
      
      /*
       * Flush chain history to file
       */
      if (chain_history_write(global.ch,
			      (ch_write_t)fwrite,
			      fp_ch) < 0) {
	fprintf(stderr, "error: failed to write chain history segment to file\n");
	return -1;
      }
      
      if (chain_history_reset(global.ch) < 0) {
	fprintf(stderr, "error: failed to reset chain history\n");
	return -1;
      }
      
    }

      
    chain_history_change_t step;
    
    if (wavetree2d_sub_get_last_perturbation(global.wt, &step) < 0) {
      fprintf(stderr, "error: failed to get last step\n");
      return -1;
    }
    
    step.header.likelihood = global.current_likelihood;
    step.header.temperature = global.temperature;
    step.header.hierarchical = global.lambda_scale;
    if (chain_history_add_step(global.ch, &step) < 0) {
      fprintf(stderr, "error: failed to add step to chain history\n");
      return -1;
    }
  
    //
    // Hierarchical
    //

    if (hierarchical != nullptr) {

      if (hierarchical->step() < 0) {
	fprintf(stderr, "error: failed to do hierarchical step\n");
	return -1;
      }

      hierarchical->get_last_step(&step);

      step.header.likelihood = global.current_likelihood;
      step.header.temperature = global.temperature;
      step.header.hierarchical = global.lambda_scale;
      
      if (chain_history_full(global.ch)) {
	
	/*
	 * Flush chain history to file
	 */
	if (chain_history_write(global.ch,
				(ch_write_t)fwrite,
				fp_ch) < 0) {
	  ERROR("error: failed to write chain history segment to file\n");
	  return -1;
	}
	
	if (chain_history_reset(global.ch) < 0) {
	  ERROR("error: failed to reset chain history\n");
	  return -1;
	}
	  
      }
      
      if (chain_history_add_step(global.ch, &step) < 0) {
	ERROR("error: failed to add step to chain history\n");
	return -1;
      }
    }

    int current_k = wavetree2d_sub_coeff_count(global.wt);
    
    if (verbosity > 0 && (i + 1) % verbosity == 0) {

      INFO("%6d: %f %d dc %f lambda %f:\n",
	   i + 1,
	   global.current_likelihood,
	   current_k,
	   wavetree2d_sub_dc(global.wt),
	   global.lambda_scale);

      INFO(birth.write_long_stats().c_str());
      INFO(death.write_long_stats().c_str());
      INFO(value.write_long_stats().c_str());
      if (hierarchical != nullptr) {
	INFO(hierarchical->write_long_stats().c_str());
      }
    }

    khistogram[current_k - 1] ++;

  }

  filename = mkfilename(output_prefix, "khistogram.txt");
  FILE *fp = fopen(filename.c_str(), "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create khistogram file\n");
    return -1;
  }
  for (int i = 0; i < kmax; i ++) {
    fprintf(fp, "%d %d\n", i + 1, khistogram[i]);
  }
  fclose(fp);

  /*
   * If there are remaining steps to save
   */
  if (chain_history_nsteps(global.ch) > 1) {
    /*
     * Flush chain history to file
     */
    if (chain_history_write(global.ch,
			    (ch_write_t)fwrite,
			    fp_ch) < 0) {
      fprintf(stderr, "error: failed to write chain history segment to file\n");
      return -1;
    }
  }
  fclose(fp_ch);

  filename = mkfilename(output_prefix, "acceptance.txt");
  fp = fopen(filename.c_str(), "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create acceptance file\n");
    return -1;
  }
  fprintf(fp, birth.write_long_stats().c_str());
  fprintf(fp, "\n");
  fprintf(fp, death.write_long_stats().c_str());
  fprintf(fp, "\n");
  fprintf(fp, value.write_long_stats().c_str());
  fprintf(fp, "\n");
  if (hierarchical != nullptr) {
    fprintf(fp, hierarchical->write_long_stats().c_str());
    fprintf(fp, "\n");
  }
  
  fclose(fp);
  
  filename = mkfilename(output_prefix, "final_model.txt");
  if (wavetree2d_sub_save(global.wt, filename.c_str()) < 0) {
    fprintf(stderr, "error: failed to save final model\n");
    return -1;
  }

  filename = mkfilename(output_prefix, "residuals.txt");
  int nres = global.get_residual_size();
  const double *res = global.get_mean_residuals();
  fp = fopen(filename.c_str(), "w");
  if (fp == NULL) {
    ERROR("Failed to create residuals file");
    return -1;
  }
  for (int i = 0; i < nres; i ++) {
    fprintf(fp, "%.9g\n", res[i]);
  }
  fclose(fp);
  
  filename = mkfilename(output_prefix, "residuals_normed.txt");
  res = global.get_mean_normed_residuals();
  fp = fopen(filename.c_str(), "w");
  if (fp == NULL) {
    ERROR("Failed to create residuals squared file");
    return -1;
  }
  for (int i = 0; i < nres; i ++) {
    fprintf(fp, "%.9g\n", res[i]);
  }
  fclose(fp);
  
  filename = mkfilename(output_prefix, "residuals_hist.txt");
  if (!global.save_residual_histogram(filename.c_str())) {
    ERROR("Failed to save residual histogram");
    return -1;
  }
  
  filename = mkfilename(output_prefix, "residuals_cov.txt");
  if (!global.save_residual_covariance(filename.c_str())) {
    ERROR("Failed to save residual covariance");
    return -1;
  }

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <file>               Input observations file\n"
	  " -I|--initial <file>             Starting model file\n"
	  " -s|--stm <file>                 Forward model information file (may be more than 1)\n"
	  " -M|--prior-file <file>          Prior/Proposal file\n"
	  " -o|--output <path>              Output prefix for output files\n"
	  "\n"
	  " -d|--degree-depth <int>         Number of vertical layers expressed as power of 2\n"
	  " -l|--degree-lateral <int>       Number of horizontal points expressed as power of 2\n"
	  "\n"
	  " -D|--depth <float>              Depth to half-space (m)\n"
	  "\n"
	  " -t|--total <int>                Total number of iterations\n"
	  " -S|--seed <int>                 Random number seed\n"
	  "\n"
	  " -H|--hierarchical <filename>    Hierarchical model filename (one for each stm file)\n"
	  " -L|--lambda-std <float>         Std deviation for lambda scaling sampling\n"
	  "\n"
	  " -k|--kmax <int>                 Max. no. of coefficients\n"
	  "\n"
	  " -B|--birth-probability <float>  Birth probability\n"
	  " -P|--posteriork                 Posterior k simulation\n"
	  "\n"
	  " -w|--wavelet-vertical <int>     Wavelet basis to use for vertical direction\n"
	  " -W|--wavelet-horizontal <int>   Wavelet basis to use for horizontal direction\n"
	  "\n"
	  " -v|--verbosity <int>            Number steps between status printouts (0 = disable\n"
	  " -h|--help                       Show usage information\n"
	  "\n",
	  pname);
}
