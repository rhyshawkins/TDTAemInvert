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


static char short_options[] = "i:I:o:r:s:D:d:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"image", required_argument, 0, 'I'},
  {"output", required_argument, 0, 'o'},
  {"response", required_argument, 0, 'r'},
  {"stm", required_argument, 0, 's'},

  {"depth", required_argument, 0, 'D'},
  {"degree-depth", required_argument, 0, 'd'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

static int load_image(const char *filename, int width, int height, double *img);


int main(int argc, char *argv[])
{
  int c;
  int option_index;

  char *input_obs;
  char *input_image;
  char *output;
  char *output_response;
  
  std::vector<std::string> stm_files;

  int degree_depth;
  double depth;

  input_obs = nullptr;
  input_image = nullptr;
  output = nullptr;
  output_response = nullptr;

  degree_depth = 5;
  depth = 200.0;

  option_index = 0;
  while (1) {
    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch(c) {

    case 'i':
      input_obs = optarg;
      break;

    case 'I':
      input_image = optarg;
      break;

    case 's':
      stm_files.push_back(optarg);
      break;

    case 'o':
      output = optarg;
      break;

    case 'r':
      output_response = optarg;
      break;

    case 'D':
      depth = atof(optarg);
      if (depth <= 0.0) {
	fprintf(stderr, "error: depth must be greater than 0\n");
	return -1;
      }
      break;

    case 'd':
      degree_depth = atoi(optarg);
      if (degree_depth < 0) {
	fprintf(stderr, "error: degree depth must be greater than 0\n");
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      
      return -1;
    }
  }

  if (input_obs == NULL) {
    fprintf(stderr, "error: require the input of a observations file\n");
    return -1;
  }

  if (input_image == NULL) {
    fprintf(stderr, "error: require the input of an image file\n");
    return -1;
  }

  if (stm_files.size() == 0) {
    fprintf(stderr, "error: require at least on stm file\n");
    return -1;
  }

  //
  // First load the observations
  //
  aemobservations obs(input_obs);

  int width = (int)obs.points.size();
  int height = 1 << degree_depth;
  
  printf("%d observations\n", width);
  printf("%d layers\n", height);

  aemimage image(height, width, depth, 1.0);

  if (load_image(input_image, width, height, image.conductivity) < 0) {
    fprintf(stderr, "error: failed to load image\n");
    return -1;
  }

  std::vector<cTDEmSystem*> forwardmodel;
  for (auto &s : stm_files) {
    cTDEmSystem *p = new cTDEmSystem(s);
    forwardmodel.push_back(p);
  }

  cEarth1D earth1d;

  earth1d.conductivity.resize(image.rows);
  earth1d.thickness.resize(image.rows - 1);

  for (int i = 0; i < (image.rows - 1); i ++) {
    earth1d.thickness[i] = image.layer_thickness[i];
  }

  FILE *fp_out = stdout;
  if (output != nullptr) {
    fp_out = fopen(output, "w");
    if (fp_out == NULL) {
      fprintf(stderr, "error: failed to create output file\n");
      return -1;
    }
  }

  FILE *fp_response = nullptr;
  if (output_response != nullptr) {
    fp_response = fopen(output_response, "w");
    if (fp_response == NULL) {
      fprintf(stderr, "error: failed to create reponse output file\n");
      return -1;
    }
  }

  for (int i = 0; i < image.columns; i ++) {

    aempoint &p = obs.points[i];

    cTDEmGeometry geometry(p.tx_height,
			   p.tx_roll,
			   p.tx_pitch,
			   0.0,
			   p.txrx_dx,
			   0.0,
			   p.txrx_dz,
			   0.0,
			   0.0,
			   0.0);

    for (int j = 0; j < image.rows; j ++) {
      earth1d.conductivity[j] = exp(image.conductivity[j * image.columns + i]);
    }

    
    for (int k = 0; k < (int)forwardmodel.size(); k ++) {

      cTDEmSystem *f = forwardmodel[k];
      const aemresponse &r = p.responses[k];

      cTDEmResponse response;

      f->forwardmodel(geometry, earth1d, response);

      switch (r.d) {
      case aemresponse::DIRECTION_Z:
	if (r.response.size() != response.SZ.size()) {
	  throw AEMEXCEPTION("Size mismatch\n");
	}

	fprintf(fp_out, "%d ", (int)r.response.size());
	for (int l = 0; l < (int)response.SZ.size(); l ++) {
	  double dx = r.response[l] - response.SZ[l];
	  fprintf(fp_out, "%.9g ", dx);
	}

	if (fp_response != nullptr) {
	  fprintf(fp_response, "%d ", (int)response.SZ.size());
	  for (int l = 0; l < (int)response.SZ.size(); l ++) {
	    fprintf(fp_response, "%.9g ", response.SZ[l]);
	  }
	}

	break;
      default:
	throw AEMEXCEPTION("Unimplemented\n");
      }
    }

    fprintf(fp_out, "\n");
    if (output_response != nullptr) {
      fprintf(fp_response, "\n");
    }
  }

  if (output != nullptr) {
    fclose(fp_out);
  }

  if (output_response != nullptr) {
    fclose(fp_response);
  }
  return 0;

}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i | --input <filename>           Input observations\n"
	  " -I | --input-image <filename>     Raw input image\n"
	  "\n"
	  " -s | --stm <filename>             Stm File(s) for forward model\n"
	  "\n"
	  " -o | --output <filename>          Residuals file\n"
	  " -r | --response <filename>        Output computed response\n"
	  "\n"
	  " -d | --degree-depth <int>         Height of image as power of 2\n"
	  " -D | --depth <float>              Physical depth of model in metres\n"
	  "\n"
	  " -h | --help                       Show usage\n"
	  "\n",
	  pname);
}

static int load_image(const char *filename, int width, int height, double *img)
{
  FILE *fp;
  int i;
  int j;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "load_image: failed to open file\n");
    return -1;
  }

  for (j = 0; j < height; j ++) {
    for (i = 0; i < width; i ++) {

      if (fscanf(fp, "%lf", &(img[j*width + i])) != 1) {
        fprintf(stderr, "load_image: failed to read pixel\n");
        return -1;
      }

    }
  }

  fclose(fp);
  return 0;
}

