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

#include <getopt.h>

#include "aemobservations.hpp"
#include "aemimage.hpp"
#include "hierarchicalmodel.hpp"
#include "aemutil.hpp"

#include "rng.hpp"

#include "tdemsystem.h"
#include "general_types.h"

static char short_options[] = "i:I:S:o:O:N:s:h";

static struct option long_options[] = {
  {"input-image", required_argument, 0, 'i'},
  {"input-path", required_argument, 0, 'I'},
  {"input-stm", required_argument, 0, 'S'},

  {"output", required_argument, 0, 'o'},
  {"output-true", required_argument, 0, 'O'},

  {"noise", required_argument, 0, 'N'},
  {"seed", required_argument, 0, 's'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};



static void usage(const char *pname);

  
int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input_image;
  char *input_path;

  std::vector<std::string> input_stm;
  std::vector<std::string> input_noise;

  char *output_file;
  char *output_true;

  int seed;

  //
  // Defaults
  //
  input_image = nullptr;
  input_path = nullptr;

  output_file = nullptr;
  output_true = nullptr;

  seed = 983;
  
  
  while (true) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input_image = optarg;
      break;

    case 'I':
      input_path = optarg;
      break;

    case 'S':
      input_stm.push_back(optarg);
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'O':
      output_true = optarg;
      break;

    case 'N':
      input_noise.push_back(optarg);
      break;

    case 's':
      seed = atoi(optarg);
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;
    }
  }

  if (input_image == nullptr) {
    fprintf(stderr, "error: required input image parameter missing\n");
    return -1;
  }

  if (input_path == nullptr) {
    fprintf(stderr, "error: required input flight path parameter missing\n");
    return -1;
  }

  if (input_stm.size() == 0) {
    fprintf(stderr, "error: required input stm file parameter missing\n");
    return -1;
  }

  if (input_noise.size() != input_stm.size()) {
    fprintf(stderr, "error: not enough noise parameters for stm's\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }
      
  aemobservations obs(input_path);

  aemimage image;
  if (!image.load(input_image)) {
    fprintf(stderr, "error: failed to load image file\n");
    return -1;
  }

  if (image.columns != (int)obs.points.size()) {
    fprintf(stderr, "error: mismatch between image columns and number of observations: %d %d\n",
	    image.columns,
	    (int)obs.points.size());
    return -1;
  }
	    

  std::vector<cTDEmSystem*> forwardmodel;
  std::vector<double*> forwardmodel_time;
  
  for (auto &s : input_stm) {

    cTDEmSystem *p = new cTDEmSystem(s);

    forwardmodel.push_back(p);

    double *centre_time = new double[p->WinSpec.size()];

    int t = 0;
    for (auto &w : p->WinSpec) {
      centre_time[t] = (w.TimeLow + w.TimeHigh)/2.0;
      t ++;
    }

    forwardmodel_time.push_back(centre_time);
  }

  std::vector<hierarchicalmodel*> noise;
  std::vector<double> lambda;

  for (auto &s : input_noise) {

    hierarchicalmodel *h = hierarchicalmodel::load(s.c_str());
    if (h == nullptr) {
      fprintf(stderr, "error: failed to load hierarchical noise model\n");
      return -1;
    }

    noise.push_back(h);
  }
  
  cEarth1D earth1d;

  earth1d.conductivity.resize(image.rows);
  earth1d.thickness.resize(image.rows - 1);

  for (int i = 0; i < (image.rows - 1); i ++) {
    earth1d.thickness[i] = image.layer_thickness[i];
  }

  std::vector<aempoint>::iterator i;
  int j;
  for (i = obs.points.begin(), j = 0; i != obs.points.end(); i ++, j ++) {
    printf("  Computing %6d\n", j);
    
    //
    // Set conductivity from image, thicknesses already set outside of loop
    //
    for (int k = 0; k < image.rows; k ++) {
      earth1d.conductivity[k] = image.conductivity[k * image.columns + j];
    }

    //
    // Reset observations X, Z
    //
    aempoint &p = *i;

    //
    // Construct geometry
    //
    cTDEmGeometry geometry(p.tx_height,
			   p.tx_roll,
			   p.tx_pitch,
			   p.tx_yaw,
			   p.txrx_dx,
			   p.txrx_dy,
			   p.txrx_dz,
			   p.rx_roll,
			   p.rx_pitch,
			   p.rx_yaw);

    for (auto &f: forwardmodel) {

      cTDEmResponse response;
    
      f->forwardmodel(geometry, earth1d, response);

      //
      // Add in the response
      //
      aemresponse r(aemresponse::DIRECTION_Z);

      for (auto &i: response.SZ) {
	r.response.push_back(i);
      }

      p.responses.push_back(r);
    }
  }

  printf("  Done\n");

  //
  // Output true observations is required
  //
  if (output_true) {

    if (!obs.save(output_true)) {
      fprintf(stderr, "error: failed to save true observations\n");
      return -1;
    }
  }

  //
  // Add the noise
  //
  Rng random(seed);
    
  for (auto &a: obs.points) {

    int ni = 0;
    
    for (auto &r: a.responses) {

      hierarchicalmodel *noisemodel = noise[ni];
      double *time = forwardmodel_time[ni];
      int di = 0;

      for (auto &d: r.response) {
	double t = time[di];

	double sigma = noisemodel->noise(d, t, 1.0);
	d = d + random.normal(sigma);
	di ++;
      }

      ni ++;
    }
  }

  if (!obs.save(output_file)) {
    fprintf(stderr, "error: failed to save noisy observations\n");
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
	  " -i|--input-image <filename>        Input true image\n"
	  " -I|--input-path <filename>         Input flight path\n"
	  " -S|--input-stm <filename>          Input Source Parameters file\n"
	  "\n"
	  " -o|--output <filename>             Output file to write (required)\n"
	  " -O|--output-true <filename>        Output true observations to file\n"
	  "\n"
	  " -n|--noise <float>                 Std dev of Gaussian noise\n"
	  " -s|--seed <int>                    Random seed\n"
	  "\n"
	  " -h|--help                          Usage information\n"
	  "\n",
	  pname);
}

  
