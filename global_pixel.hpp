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

#pragma once
#ifndef global_pixel_hpp
#define global_pixel_hpp

#include <vector>
#include <string>

#include <gmp.h>

#include <mpi.h>

#include "rng.hpp"
#include "aemimage.hpp"
#include "aemobservations.hpp"

#include "tdemsystem.h"
#include "general_types.h"

#include "constants.hpp"

#include "chainhistory_pixel.hpp"
#include "hierarchicalmodel.hpp"

class GlobalPixel {
public:

  GlobalPixel(const char *filename,
	      const std::vector<std::string> &stm_files,
	      const char *initial_model,
	      double prior_min,
	      double prior_max,
	      double proposal_stddev,
	      int degreex,
	      int degreey,
	      double depth,
	      int hierarchical,
	      const std::vector<double> &initial_lambda,
	      int seed);
  ~GlobalPixel();

  double likelihood();

  void initialize_mpi(MPI_Comm communicator);
  
  double likelihood_mpi();

  void load_initial_model(const char *filename);

  double depth;

  int degreex;
  int degreey;

  std::vector<cTDEmSystem*> forwardmodel;
  std::vector<double*> forwardmodel_time;
  std::vector<hierarchicalmodel*> lambda;
  double lambda_scale;
  
  aemobservations *observations;
  aemimage *image;

  int width;
  int height;
  int size;

  double current_likelihood;

  Rng random;

  double prior_min;
  double prior_max;

  double proposal_stddev;

  ChainHistoryPixel *chainhistory;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

};


#endif // global_pixel_hpp
