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

#include <getopt.h>

#include "global.hpp"

extern "C" {
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"

#include "wavetree2d_sub.h"
}

static char short_options[] = "i:o:s:D:d:l:w:W:H:L:nh";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"observations", required_argument, 0, 'o'},
  
  {"stm", required_argument, 0, 's'},
  {"depth", required_argument, 0, 'D'},
  
  {"degree-depth", required_argument, 0, 'd'},
  {"degree-lateral", required_argument, 0, 'l'},

  {"wavelet-vertical", required_argument, 0, 'w'},
  {"wavelet-horizontal", required_argument, 0, 'W'},

  {"hierarchical", required_argument, 0, 'H'},

  {"lambda", required_argument, 0, 'L'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input_model;
  char *input_obs;
  std::vector<std::string> stm_files;
  std::vector<std::string> hierarchical_files;
  int degreex;
  int degreey;

  double depth;

  int wavelet_v;
  int wavelet_h;

  double lambda_scale;

  //
  // Defaults
  //

  input_model = nullptr;
  input_obs = nullptr;

  degreex = 10;
  degreey = 5;

  depth = 500.0;

  wavelet_v = 0;
  wavelet_h = 0;

  lambda_scale = 1.0;

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
      input_model = optarg;
      break;
      
    case 'o':
      input_obs = optarg;
      break;

    case 's':
      stm_files.push_back(optarg);
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
      
    case 'H':
      hierarchical_files.push_back(optarg);
      break;

    case 'L':
      lambda_scale = atof(optarg);
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

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  // if (input_model == nullptr) {
  //   fprintf(stderr, "error: required input parameter input model missing\n");
  //   return -1;
  // }
  
  if (input_obs == nullptr) {
    fprintf(stderr, "error: required input parameter input observations missing\n");
    return -1;
  }

  if (stm_files.size() == 0) {
    fprintf(stderr, "error: need at least on stm file\n");
    return -1;
  }

  if (stm_files.size() != hierarchical_files.size()) {
    fprintf(stderr, "error: mismatch in size of hierarchical and stm lists\n");
    return -1;
  }


  Global global(input_obs,
		stm_files,
		input_model,
		nullptr,
		degreex,
		degreey,
		depth,
		hierarchical_files,
		0,
		100,
		false,
		wavelet_h,
		wavelet_v);

  global.lambda_scale = lambda_scale;

  double log_normalization;
  double like = global.likelihood(log_normalization);
  printf("Likelihood: %g (%g)\n", like, log_normalization);

  return 0;
}
  
static void
usage(const char *pname)
{
  fprintf(stderr, "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  "-i|--input <filename>          Input wavetree model file\n"
	  "-o|--observations <filename>   Input observations file\n"
	  "-s|--stm <filename>            Input STM file\n"
	  "\n"
	  "-D|--depth <float>             Depth in metres\n"
	  "-d|--degree-depth <int>        No. layers as a power of 2\n"
	  "-l|--degree-lateral <int>      No. lateral points as power of 2\n"
	  "\n"
	  "-w|--wavelet-vertical <int>    Wavelet in depth direction\n"
	  "-W|--wavelet-horizontal <int>  Wavelet in lateral direction\n"
	  "\n"
	  "-H|--hierarhical <filename>    Hierachical model filename (one for each stm)\n"
	  "\n"
	  "-h|--help                      Usage\n"
	  "\n",
	  pname);
}
