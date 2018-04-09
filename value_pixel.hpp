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
#ifndef value_pixel_hpp
#define value_pixel_hpp

#include "global_pixel.hpp"

class PixelPerturbation;

class ValuePixel {
public:

  ValuePixel(GlobalPixel &global);
  ~ValuePixel();

  int step(PixelPerturbation &pb);

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  GlobalPixel &global;

  int propose;
  int accept;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;
  
  bool primary() const;

  int choose_value_location_and_value(int &value_idx,
				      double &value,
				      int &valid_proposal);

  int communicate_value_location_and_value(int &valid_proposal,
					   int &value_idx,
					   double &value);



  int compute_likelihood(int value_idx, double &proposed_likelihood);

  int compute_acceptance(int value_idx,
			 double value_prior_ratio,
			 double proposed_likelihood,
			 bool &accept_proposal);
  
  int communicate_acceptance(bool &accept_proposal);

};

#endif // value_pixel_hpp
