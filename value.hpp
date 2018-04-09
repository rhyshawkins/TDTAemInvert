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
#ifndef value_hpp
#define value_hpp

#include <mpi.h>

#include "global.hpp"

class Value {
public:

  Value(Global &global);
  ~Value();

  int step();

  std::string write_short_stats();

  std::string write_long_stats();

  void initialize_mpi(MPI_Comm communicator);

  Global &global;

  int propose;
  int accept;

  int *propose_depth;
  int *accept_depth;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;

private:

  bool primary() const;

  int choose_value_location_and_value(int &value_depth,
				      int &value_idx,
				      double &choose_prob,
				      double &value,
				      int &ii,
				      int &ij,
				      double &value_prior_ratio,
				      int &prior_errors,
				      int &valid_proposal);

  int communicate_value_location_and_value(int &valid_proposal,
					   int &value_idx,
					   int &value_depth,
					   double &value);

  int propose_value(int valid_proposal,
		    int value_idx,
		    int value_depth,
		    double value);

  int compute_likelihood(int valid_idx, double &proposed_likelihood, double &proposed_log_normalization);

  int compute_acceptance(int value_idx,
			 double value_prior_ratio,
			 double proposed_likelihood,
			 double proposed_log_normalization,
			 bool &accept_proposal);
  
  int communicate_acceptance(bool &accept_proposal);
};

#endif // value_hpp
