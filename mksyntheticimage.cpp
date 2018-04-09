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
#include <math.h>

#include <getopt.h>

#include <string>

#include "rng.hpp"

#include "aemimage.hpp"

#include "constants.hpp"

static char short_options[] = "W:H:D:m:o:O:b:c:lh";
static struct option long_options[] = {
  {"horizontal-samples", required_argument, 0, 'W'},
  {"depth-samples", required_argument, 0, 'H'},
  {"depth", required_argument, 0, 'D'},
  {"model", required_argument, 0, 'm'},
  
  {"output", required_argument, 0, 'o'},
  {"output-image", required_argument, 0, 'O'},

  {"background-conductivity", required_argument, 0, 'b'},
  {"conductivity", required_argument, 0, 'c'},

  {"list", no_argument, 0, 'l'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

typedef aemimage *(*mkimage_t)(int, int, double, double, double);

static aemimage *mkconstantimage(int hsamples,
				 int dsamples,
				 double depth,
				 double background_conductivity,
				 double conductivity);

static aemimage *mkdettmerimage(int hsamples,
				int dsamples,
				double depth,
				double background_conductivity,
				double conductivity);

static aemimage *mkdettmerpatternimage(int hsamples,
				       int dsamples,
				       double depth,
				       double background_conductivity,
				       double conductivity);

static struct _imagetable {
  std::string name;
  mkimage_t mkimage;
} imagetable[] = {

  {"constant", mkconstantimage},

  {"dettmer", mkdettmerimage},

  {"dettmerpattern", mkdettmerpatternimage},

  {"", 0}
};

static bool ispositivepower2(int i);

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Options
  //
  int hsamples;
  int dsamples;
  double depth;

  double background_conductivity;
  double conductivity;

  char *model_name;
  char *output_file;
  char *output_image;

  //
  // Defaults
  //
  hsamples = 1024;
  dsamples = 32;
  depth = 150.0;

  model_name = nullptr;
  output_file = nullptr;
  output_image = nullptr;

  background_conductivity = 0.050;
  conductivity = 0.200;

  //
  // Cmd line arguments
  //
  option_index = -1;
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
      
    case 'W':
      hsamples = atoi(optarg);
      if (!ispositivepower2(hsamples)) {
	fprintf(stderr, "error: horizontal samples must be a power of 2 greater than 0\n");
	return -1;
      }
      break;

    case 'H':
      dsamples = atoi(optarg);
      if (!ispositivepower2(dsamples)) {
	fprintf(stderr, "error: vertical samples must be a power of 2 greater than 0\n");
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

    case 'm':
      model_name = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'O':
      output_image = optarg;
      break;

    case 'b':
      background_conductivity = atof(optarg);
      if (background_conductivity < CONDUCTIVITY_MIN ||
	  background_conductivity > CONDUCTIVITY_MAX) {
	fprintf(stderr, "error: background conductivity must be between %.3f and %.3f\n",
		CONDUCTIVITY_MIN, CONDUCTIVITY_MAX);
	return -1;
      }
      break;
      
    case 'c':
      conductivity = atof(optarg);
      if (conductivity < CONDUCTIVITY_MIN ||
	  conductivity > CONDUCTIVITY_MAX) {
	fprintf(stderr, "error: background conductivity must be between %.3f and %.3f\n",
		CONDUCTIVITY_MIN, CONDUCTIVITY_MAX);
	return -1;
      }
      break;

    case 'l':
      fprintf(stderr, "Available models:\n");
      for (int i = 0; imagetable[i].name.length() != 0; i ++) {
	printf("  %s\n", imagetable[i].name.c_str());
      }
      return -1;
	

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (model_name == nullptr) {
    fprintf(stderr, "error: required parameter model name missing\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }

  
  mkimage_t mkimage = nullptr;
  for (int i = 0; imagetable[i].name.length() != 0; i ++) {
    if (imagetable[i].name == model_name) {
      mkimage = imagetable[i].mkimage;
      break;
    }
  }

  if (mkimage == nullptr) {
    fprintf(stderr, "error: no model name %s\n", model_name);
    return -1;
  }

  aemimage *image = mkimage(hsamples, dsamples, depth, background_conductivity, conductivity);
  if (image == nullptr) {
    fprintf(stderr, "error: failed to create image\n");
    return -1;
  }

  if (!image->save(output_file)) {
    fprintf(stderr, "error: failed to save image\n");
    return -1;
  }

  if (output_image != nullptr) {
    if (!image->save_image(output_image)) {
      fprintf(stderr, "error: failed to save raw image\n");
      return -1;
    }
  }

  delete image;

  return 0;
}

static bool ispositivepower2(int i)
{
  if (i <= 0) {
    return false;
  } else {

    return (i & (i - 1)) == 0;

  }
}

static aemimage *mkconstantimage(int hsamples,
				 int dsamples,
				 double depth,
				 double background_conductivity,
				 double conductivity)
{
  return new aemimage(dsamples, hsamples, depth, background_conductivity);
}

static aemimage *mkdettmerimage(int hsamples,
				int dsamples,
				double depth,
				double background_conductivity,
				double conductivity)
{
  aemimage *image = new aemimage(dsamples, hsamples, depth, 0.0);

  int vscale = dsamples/8;
  int hscale = hsamples/8;

  for (int j = 0; j < dsamples; j ++) {

    int jj = j/vscale;

    for (int i = 0; i < hsamples; i ++) {

      int ii = i/hscale;

      double c = background_conductivity;

      if (ii >= 4 && ii <= 6 &&
	  jj >= 4 && jj <= 6) {
	c = conductivity;
      }

      image->conductivity[j * hsamples + i] = c;
    }
  }

  return image;
}

static aemimage *mkdettmerpatternimage(int hsamples,
				       int dsamples,
				       double depth,
				       double background_conductivity,
				       double conductivity)
{
  aemimage *image = new aemimage(dsamples, hsamples, depth, 0.0);

  int vscale = dsamples/16;
  int hscale = hsamples/16;

  for (int j = 0; j < dsamples; j ++) {

    int jj = j/vscale;

    for (int i = 0; i < hsamples; i ++) {

      int ii = i/hscale;

      double c = background_conductivity;

      if (ii >= 8 && ii <= 13 &&
	  jj >= 8 && jj <= 13) {

	if (ii == 8 || ii == 13 ||
	    jj == 8 || jj == 13) {
	  c = conductivity;
	} else if (ii == jj || (ii == 12 && jj == 11) || (ii == 11 && jj == 12)) {
	  c = conductivity;
	}
      }

      image->conductivity[j * hsamples + i] = c;
    }
  }

  return image;
}

static void usage(const char *pname)
{
  fprintf(stderr,
          "usage: %s [options]\n"
          "where options is one or more of:\n"
          "\n"
	  " -o|--output <filename>                Output file to write (required)\n"
	  " -m|--model <model name>               Model name to generate (required)\n"
	  "\n"
	  " -W|--horizontal-samples <int>         Horizontal samples (default = 1024) must be power of 2\n"
	  " -H|--depth-samples <int>              Vertical samples (default = 32) must be power of 2\n"
	  " -D|--depth <float>                    Depth in metres (default 500.0)\n"
	  "\n"
	  " -b|--background-conductivity <float>  Background or average conductivity\n"
	  " -c|--conductivity <float>             Conductivity of synthetic anomaly\n"
	  "\n"
	  " -l|--list                             List available models and exit\n"
	  "\n"
	  " -h|--help                             Usage information\n"
	  "\n",
	  pname);
}

