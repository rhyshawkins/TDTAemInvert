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

#include "rng.hpp"

#include "aemobservations.hpp"

static char short_options[] = "N:e:E:p:P:r:R:x:X:z:Z:S:o:h";
static struct option long_options[] = {
  {"nsamples", required_argument, 0, 'N'},
  {"height-mean", required_argument, 0, 'e'},
  {"height-std", required_argument, 0, 'E'},
  {"pitch-mean", required_argument, 0, 'p'},
  {"pitch-std", required_argument, 0, 'P'},
  {"roll-mean", required_argument, 0, 'r'},
  {"roll-std", required_argument, 0, 'R'},
  {"dx-mean", required_argument, 0, 'x'},
  {"dx-std", required_argument, 0, 'X'},
  {"dz-mean", required_argument, 0, 'z'},
  {"dz-std", required_argument, 0, 'Z'},
  {"seed", required_argument, 0, 'S'},
  
  {"output", required_argument, 0, 'o'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static double random_walk_init(Rng &random,
			       double mu,
			       double sigma);

static double random_walk_step(Rng &random,
			       double x0,
			       double mu,
			       double sigma,
			       double scale = 10.0);

static bool ispositivepower2(int i);

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  int N;
  
  double height_mean;
  double height_std;

  double pitch_mean;
  double pitch_std;

  double roll_mean;
  double roll_std;

  double dx_mean;
  double dx_std;

  double dz_mean;
  double dz_std;

  int seed;

  char *output_file;
  
  //
  // State
  //
  double height;
  double roll;
  double pitch;
  double dx;
  double dz;

  //
  // Defaults
  //

  N = 1024;

  height_mean = 100.0;
  height_std = 5.0;

  roll_mean = 0.0;
  roll_std = 2.0;

  pitch_mean = 0.0;
  pitch_std = 1.0;

  dx_mean = -100.0;
  dx_std = 2.0;

  dz_mean = -40.0;
  dz_std = 2.5;

  seed = 983;

  output_file = nullptr;

  option_index = 0;
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'N':
      N = atoi(optarg);
      if (!ispositivepower2(N)) {
	fprintf(stderr, "error: no. samples must be a power of 2 greater than 0\n");
	return -1;
      }
      break;
      
    case 'e':
      height_mean = atof(optarg);
      break;

    case 'E':
      height_std = atof(optarg);
      if (height_std < 0.0) {
	fprintf(stderr, "error: height std must be 0 or greater\n");
	return -1;
      }
      break;

    case 'p':
      pitch_mean = atof(optarg);
      break;

    case 'P':
      pitch_std = atof(optarg);
      if (pitch_std < 0.0) {
	fprintf(stderr, "error: pitch std must be 0 or greater\n");
	return -1;
      }
      break;

    case 'r':
      roll_mean = atof(optarg);
      break;

    case 'R':
      roll_std = atof(optarg);
      if (roll_std < 0.0) {
	fprintf(stderr, "error: roll std must be 0 or greater\n");
	return -1;
      }
      break;

    case 'x':
      dx_mean = atof(optarg);
      break;

    case 'X':
      dx_std = atof(optarg);
      if (dx_std < 0.0) {
	fprintf(stderr, "error: dx std must be 0 or greater\n");
	return -1;
      }
      break;
      
    case 'z':
      dz_mean = atof(optarg);
      break;

    case 'Z':
      dz_std = atof(optarg);
      if (dz_std < 0.0) {
	fprintf(stderr, "error: dz std must be 0 or greater\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;
      
    case 'o':
      output_file = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required parameter output file missing\n");
    return -1;
  }
  
  Rng random(seed);
  aemobservations obs;

  if (height_std == 0.0) {
    height = height_mean;
  } else {
    height = height_mean + random.normal(height_std);
  }

  height = random_walk_init(random, height_mean, height_std);
  
  roll = random_walk_init(random, roll_mean, roll_std);
  pitch = random_walk_init(random, pitch_mean,pitch_std);
  dx = random_walk_init(random, dx_mean, dx_std);
  dz = random_walk_init(random, dz_mean, dz_std);

  for (int i = 0; i < N; i ++) {

    obs.points.push_back(aempoint(height,
				  roll, pitch, 0.0,
				  dx, 0.0, dz,
				  roll, pitch, 0.0));

    height = random_walk_step(random, height, height_mean, height_std);
    roll = random_walk_step(random, roll, roll_mean, roll_std);
    pitch = random_walk_step(random, pitch, pitch_mean, pitch_std);
    dx = random_walk_step(random, dx, dx_mean, dx_std);
    dz = random_walk_step(random, dz, dx_mean, dz_std);

  }

  if (!obs.save(output_file)) {
    fprintf(stderr, "error: failed to save output file\n");
    return -1;
  }

  return 0;
}

static double random_walk_init(Rng &random,
			       double mu,
			       double sigma)
{
  if (sigma > 0.0) {
    return mu + random.normal(sigma);
  } else {
    return mu;
  }
}

static double random_walk_step(Rng &random,
			       double x0,
			       double mu,
			       double sigma,
			       double scale)
{
  if (sigma > 0.0) {
    double x = x0 + random.normal(sigma/scale);
    double u = log(random.uniform());
    
    while (u > ((x0 - mu)*(x0 - mu)/(2.0 * sigma * sigma) - (x - mu)*(x - mu)/(2.0 * sigma * sigma))) {
      x = x0 + random.normal(sigma/scale);
      u = log(random.uniform());
    }

    return x;
  } else {
    return x0;
  }
}

static bool ispositivepower2(int i)
{
  if (i <= 0) {
    return false;
  } else {

    return (i & (i - 1)) == 0;

  }
}


static void usage(const char *pname)
{
  fprintf(stderr,
          "usage: %s [options]\n"
          "where options is one or more of:\n"
          "\n"
	  " -o|--output <filename>             Output file to write (required)\n"
	  "\n"
	  " -N|--nsamples <int>                Number of horizontal samples (must be power of 2)\n"
	  "\n"
	  " -e|--height-mean <float>           Height (m)\n"
	  " -E|--height-std <float>            Height std dev (0 = none)\n"
	  " -p|--pitch-mean <float>            Pitch\n"
	  " -P|--pitch-std <float>             Pitch std dev (0 = none)\n"
	  " -r|--roll-mean <float>             Roll\n"
	  " -R|--roll-std <float>              Roll std dev (0 = none)\n"
	  " -x|--dx-mean <float>               dx\n"
	  " -X|--dx-std <float>                dx std dev (0 = none)\n"
	  " -z|--dz-mean <float>               dz\n"
	  " -Z|--dz-std <float>                dz std dev (0 = none)\n"
	  "\n"
	  " -S|--seed <int>                    Random seed\n"
	  "\n"
	  " -h|--help                          Usage information\n"
	  "\n",
	  pname);
}
