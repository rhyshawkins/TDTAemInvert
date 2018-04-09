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

#include "value_pixel.hpp"

#include "aemutil.hpp"

ValuePixel::ValuePixel(GlobalPixel &_global) :
  global(_global),
  propose(0),
  accept(0),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1)
{
}

ValuePixel::~ValuePixel()
{
}

int
ValuePixel::step(PixelPerturbation &pb)
{
  propose ++;

  int value_idx;
  double old_value;
  double new_value;
  int valid_proposal = 0;
  double proposed_likelihood;
  bool accept_proposal = false;
  
  if (choose_value_location_and_value(value_idx, new_value, valid_proposal) < 0) {
    throw AEMEXCEPTION("Failed to choose new pixel/value");
  }

  if (communicate_value_location_and_value(valid_proposal,
					   value_idx,
					   new_value) < 0) {
    throw AEMEXCEPTION("Failed to communicate pixel perturbation");
  }

  old_value = global.image->conductivity[value_idx];
  
  pb.accepted = false;
  pb.idx = value_idx;
  pb.oldvalue = old_value;
  pb.newvalue = new_value;

  if (valid_proposal) {
    
    //
    // Perturb the model
    //
    global.image->conductivity[value_idx] = new_value;
    
    if (compute_likelihood(value_idx, proposed_likelihood) < 0) {
      throw AEMEXCEPTION("Failed to compute likelihood");
    }
    
    if (compute_acceptance(value_idx,
			   1.0,
			   proposed_likelihood,
			   accept_proposal) < 0) {
      throw AEMEXCEPTION("Failed to compute acceptance");
    }
    
    if (communicate_acceptance(accept_proposal) < 0) {
      throw AEMEXCEPTION("Failed to communicate acceptance");
    }
    
    
    if (accept_proposal) {
      
      //
      // Accept
      //
      
      accept ++;
      global.current_likelihood = proposed_likelihood;
      pb.accepted = true;
      
      return 1;
      
    } else {
      
      //
      // Reject and restore old value
      //
      global.image->conductivity[value_idx] = old_value;
      return 0;
    }
  }
 
  return 0;
}

std::string
ValuePixel::write_short_stats()
{
  return mkformatstring("ValuePixel %6d/%6d %7.3f", accept, propose, propose == 0 ? 0.0 : 100.0*(double)accept/(double)propose);
}

std::string
ValuePixel::write_long_stats()
{
  return write_short_stats();
}

void
ValuePixel::initialize_mpi(MPI_Comm _communicator)
{
  MPI_Comm_dup(_communicator, &communicator);

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw AEMEXCEPTION("MPI Failure\n");
  }
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw AEMEXCEPTION("MPI Failure\n");
  }
}

bool
ValuePixel::primary() const
{
  return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int
ValuePixel::choose_value_location_and_value(int &value_idx,
					    double &value,
					    int &valid_proposal)
{
  if (primary()) {

    value_idx = global.random.uniform(global.size);
    value = global.image->conductivity[value_idx] + global.random.normal(global.proposal_stddev);

    valid_proposal = 0;
    if (value >= global.prior_min &&
	value <= global.prior_max) {
      valid_proposal = 1;
    }
  }

  return 0;
}


int
ValuePixel::communicate_value_location_and_value(int &valid_proposal,
						 int &value_idx,
						 double &value)
{
  if (communicator != MPI_COMM_NULL) {
    if (MPI_Bcast(&valid_proposal, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw AEMEXCEPTION("Failed to broadcast valid proposal\n");
    }
    
    if (valid_proposal) {
      
      if (MPI_Bcast(&value_idx, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
	throw AEMEXCEPTION("Failed to broadcast index\n");
      }
      if (MPI_Bcast(&value, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
	throw AEMEXCEPTION("Failed to broadcast value\n");
      }
    }
  }

  return 0;
}

int
ValuePixel::compute_likelihood(int value_idx, double &proposed_likelihood)
{
  if (communicator == MPI_COMM_NULL) {
    proposed_likelihood = global.likelihood();
  } else {
    proposed_likelihood = global.likelihood_mpi();
  }

  return 0;
}

int
ValuePixel::compute_acceptance(int value_idx,
			       double value_prior_ratio,
			       double proposed_likelihood,
			       bool &accept_proposal)
{
  if (primary()) {
    
    double u = log(global.random.uniform());
    
    double alpha = (log(value_prior_ratio) +
		    (global.current_likelihood - proposed_likelihood));
    
    accept_proposal = (u < alpha);

  }

  return 0;
}

int
ValuePixel::communicate_acceptance(bool &accept_proposal)
{
  if (communicator != MPI_COMM_NULL) {

    int ta;

    if (mpi_rank == 0) {
      ta = (int)accept_proposal;
    }

    if (MPI_Bcast(&ta, 1, MPI_INT, 0, communicator) != MPI_SUCCESS) {
      throw AEMEXCEPTION("Failed to broadcast acceptted\n");
    }

    if (mpi_rank != 0) {
      accept_proposal = (bool)ta;
    }

  }

  return 0;
}


