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
#ifndef death_hpp
#define death_hpp

#include "global.hpp"

class Death {
public:

  Death(Global &global);
  ~Death();
  
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

  int choose_death_location(int k,
			    double &ratio,
			    int &death_depth,
			    int &death_idx,
			    double &choose_prob,
			    int &death_valid);

  int communicate_death_location(int &death_valid,
				 int &death_idx,
				 int &death_depth);

  int propose_death(int death_valid,
		    int death_idx,
		    int death_depth,
		    double &death_value);

  int compute_reverse_death_probability(int death_idx,
					int death_depth,
					double death_value,
					int &ii,
					int &ij,
					double &death_parent_coeff,
					double &death_prob,
					double &reverse_prob,
					double &prior_prob);

  int compute_likelihood(int death_idx,
			 double &proposed_likelihood,
			 double &proposed_log_normalization);

  int compute_acceptance(double proposed_likelihood,
			 double proposed_log_normalization,
			 double reverse_prob,
			 double choose_prob,
			 double birth_prob,
			 double ratio,
			 double prior_prob,
			 bool &accept_proposal);

  int communicate_acceptance(bool &accept_proposal);
 
};

#endif // death_hpp
