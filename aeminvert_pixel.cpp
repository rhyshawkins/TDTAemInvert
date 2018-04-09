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

#include <mpi.h>

#include "global_pixel.hpp"
#include "value_pixel.hpp"

#include "aemutil.hpp"

extern "C" {
  #include "slog.h"
};

static char short_options[] = "i:I:s:o:d:l:D:t:S:H:L:v:p:P:r:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial-model", required_argument, 0, 'I'},
  {"stm", required_argument, 0, 's'},
  {"output", required_argument, 0, 'o'},
  
  {"degree-depth", required_argument, 0, 'd'},
  {"degree-lateral", required_argument, 0, 'l'},

  {"depth", required_argument, 0, 'D'},

  {"total", required_argument, 0, 't'},
  {"seed", required_argument, 0, 'S'},

  {"hierarchical", required_argument, 0, 'H'},
  {"lambda", required_argument, 0, 'L'},

  {"verbosity", required_argument, 0, 'v'},

  {"prior-min", required_argument, 0, 'p'},
  {"prior-max", required_argument, 0, 'P'},
  {"proposal-stddev", required_argument, 0, 'r'},

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
  char *output_prefix;

  int degreex;
  int degreey;

  double depth;

  int total;
  int seed;

  int hierarchical;
  std::vector<double> initial_lambda;

  double prior_min;
  double prior_max;
  double proposal_stddev;

  int verbosity;

  std::string filename;
  FILE *fp;

  int mpi_rank;
  int mpi_size;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  //
  // Defaults
  //

  input_obs = nullptr;
  initial_model = nullptr;
  output_prefix = nullptr;

  degreex = 10;
  degreey = 5;

  depth = 500.0;

  total = 10000;
  seed = 983;

  hierarchical = 0;

  prior_min = -3.0;
  prior_max = 0.5;
  proposal_stddev = 0.1;

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
      hierarchical = atoi(optarg);
      if (hierarchical < 0 || hierarchical > 1) {
	fprintf(stderr, "error: hierarchical model must be 0 or 1\n");
	return -1;
      }
      break;
      
    case 'L':
      if (!loadhierarchicallambda(optarg, initial_lambda)) {
	fprintf(stderr, "error: failed to load hierarchical lambda file\n");
	return -1;
      }
      break;

    case 'v':
      verbosity = atoi(optarg);
      break;

    case 'p':
      prior_min = atof(optarg);
      break;

    case 'P':
      prior_max = atof(optarg);
      break;

    case 'r':
      proposal_stddev = atof(optarg);
      if (proposal_stddev <= 0.0) {
	fprintf(stderr, "error: proposal std dev must be greater than 0\n");
	return -1;
      }
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

  GlobalPixel global(input_obs,
		     stm_files,
		     initial_model,
		     prior_min,
		     prior_max,
		     proposal_stddev,
		     degreex,
		     degreey,
		     depth,
		     hierarchical,
		     initial_lambda,
		     seed);

  ValuePixel value(global);
  PixelPerturbation pb;

  if (mpi_size > 1) {
    global.initialize_mpi(MPI_COMM_WORLD);
    value.initialize_mpi(MPI_COMM_WORLD);

    global.current_likelihood = global.likelihood_mpi();
  } else {
    global.current_likelihood = global.likelihood();
  }
    

  if (mpi_rank == 0) {
    INFO("Initial Likelihood: %f\n", global.current_likelihood);
  }

  for (int i = 0; i < total; i ++) {

    //
    // Value
    //
    if (value.step(pb) < 0) {
      fprintf(stderr, "error: failed to do value step\n");
      return -1;
    }

    if (mpi_rank == 0) {
      if (verbosity > 0 && (i + 1) % verbosity == 0) {
	
	INFO("%6d: %f: %s",
	       i + 1,
	     global.current_likelihood,
	     value.write_long_stats().c_str());
	
      }
      global.chainhistory->history.push_back(pb);
    }
  }

  if (mpi_rank == 0) {
    filename = mkfilename(output_prefix, "acceptance.txt");
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL) {
      fprintf(stderr, "error: failed to create acceptance file\n");
      return -1;
    }
    
    fprintf(fp, "%s\n", value.write_long_stats().c_str());
    
    fclose(fp);

    filename = mkfilename(output_prefix, "ch.dat");
    if (!global.chainhistory->save(filename.c_str())) {
      fprintf(stderr, "error: failed to write chain history\n");
      return -1;
    }
  }

  MPI_Finalize();
	    
  
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
	  " -L|--lambda <float>             Fixed noise level\n"
	  "\n"
	  " -p|--prior-min <float>          Uniform prior min value\n"
	  " -P|--prior-max <float>          Uniform prior max value\n"
	  " -r|--proposal-stddev <float>    Std dev for value proposals\n"
	  "\n"
	  " -v|--verbosity <int>            Number steps between status printouts (0 = disable\n"
	  " -h|--help                       Show usage information\n"
	  "\n",
	  pname);
}
