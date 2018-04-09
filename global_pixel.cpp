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

#include "aemexception.hpp"
#include "aemobservations.hpp"

#include "global_pixel.hpp"

#include "constants.hpp"

extern "C" {
  #include "slog.h"
};

GlobalPixel::GlobalPixel(const char *filename,
			 const std::vector<std::string> &stm_files,
			 const char *initial_model,
			 double _prior_min,
			 double _prior_max,
			 double _proposal_stddev,
			 int _degreex,
			 int _degreey,
			 double _depth,
			 int hierarchical,
			 const std::vector<double> &initial_lambda,
			 int seed) :
  depth(_depth),
  degreex(_degreex),
  degreey(_degreey),
  lambda_scale(1.0),
  observations(nullptr),
  image(nullptr),
  width(-1),
  height(-1),
  size(-1),
  current_likelihood(-1.0),
  random(seed),
  prior_min(_prior_min),
  prior_max(_prior_max),
  proposal_stddev(_proposal_stddev),
  chainhistory(nullptr),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1)
{
  if (degreex < 0 || degreex >= 16 ||
      degreey < 0 || degreey >= 16) {
    throw AEMEXCEPTION("Degree(s) out of range: %d x %d\n", degreex, degreey);
  }
  
  if (depth <= 0.0) {
    throw AEMEXCEPTION("Depth out of range\n");
  }

  //
  // Load observations
  //
  observations = new aemobservations(filename);
  
  //
  // Load stm files
  //
  int hoffset = 0;
  for (auto &s : stm_files) {
    cTDEmSystem *p = new cTDEmSystem(s);
    
    forwardmodel.push_back(p);

    double *centre_time = new double[p->WinSpec.size()];

    int t = 0;
    for (auto &w : p->WinSpec) {
      centre_time[t] = (w.TimeLow + w.TimeHigh)/2.0;
      t ++;
    }
    
    forwardmodel_time.push_back(centre_time);

    hierarchicalmodel *m;
    switch(hierarchical) {
    case 0:
      INFO("IID Noise");
      m = new independentgaussianhierarchicalmodel();
      break;
      
    case 1:
      INFO("Hyperbolic Noise");
      m = new hyperbolichierarchicalmodel();
      break;
      
    default:
      throw AEMEXCEPTION("Invalid hierarchical model index\n");
    }
    
    if ((hoffset + m->nparameters()) > (int)initial_lambda.size()) {
      throw AEMEXCEPTION("Not enough lambda initialization parameters for hierarchical model(s) %d\n",
			 (int)initial_lambda.size());
    }
    
    for (int h = 0; h < m->nparameters(); h ++, hoffset ++) {
      m->setparameter(h, initial_lambda[hoffset]);
    }
    
    lambda.push_back(m);
  }
  
  if (forwardmodel.size() != observations->points[0].responses.size()) {
    throw AEMEXCEPTION("Mismatch in STM and responses size: %d != %d\n",
		       (int)forwardmodel.size(),
		       (int)observations->points[0].responses.size());
  }

  width = 1 << degreex;
  height = 1 << degreey;
  size = width * height;

  printf("Image: %d x %d\n", width, height);

  if ((int)observations->points.size() != width) {
    throw AEMEXCEPTION("Image size mismatch to observations: %d != %d\n",
		       width,
		       (int)observations->points.size());
  }
  
  image = new aemimage(height, width, depth, DEFAULT_CONDUCTIVITY);

  if (initial_model) {
    load_initial_model(initial_model);
  }

  chainhistory = new ChainHistoryPixel(*image);
  
  printf("Data: %d total points\n", observations->total_response_datapoints());

}

GlobalPixel::~GlobalPixel()
{
}

double
GlobalPixel::likelihood()
{
  cEarth1D earth1d;
  
  //
  // Setup layer thicknesses
  //
  earth1d.conductivity.resize(image->rows);
  earth1d.thickness.resize(image->rows - 1);
  
  for (int i = 0; i < (image->rows - 1); i ++) {
    earth1d.thickness[i] = image->layer_thickness[i];
  }
  
  double sum = 0.0;
  
  for (int i = 0; i < image->columns; i ++) {
    
    aempoint &p = observations->points[i];
    
    //
    // Construct geometry
    //
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
    //
    // Copy image column to earth model, our model is in log of conductivity so here we use exp
    //
    for (int j = 0; j < image->rows; j ++) {
      earth1d.conductivity[j] = exp(image->conductivity[j * image->columns + i]);
    }
    
    for (int k = 0; k < (int)forwardmodel.size(); k ++) {

      cTDEmSystem *f = forwardmodel[k];
      hierarchicalmodel *h = lambda[k];
      double *time = forwardmodel_time[k];
      const aemresponse &r = p.responses[k];
      
      cTDEmResponse response;
	
      f->forwardmodel(geometry,
		      earth1d,
		      response);
      
      switch (r.d) {
      case aemresponse::DIRECTION_X:
	if (r.response.size() != response.SX.size()) {
	  throw AEMEXCEPTION("Size mismatch in X response\n");
	}
	for (int l = 0; l < (int)response.SX.size(); l ++) {
	  
	  double dx = r.response[l] - response.SX[l];
	  double noise = h->noise(r.response[l], time[l], lambda_scale);
	  
	  sum += dx*dx/(2.0 * noise * noise);
	}
	break;
	
      case aemresponse::DIRECTION_Y:
	if (r.response.size() != response.SY.size()) {
	  throw AEMEXCEPTION("Size mismatch in Y response\n");
	}
	for (int l = 0; l < (int)response.SY.size(); l ++) {
	  
	  double dy = r.response[l] - response.SY[l];
	  double noise = h->noise(r.response[l], time[l], lambda_scale);
	  
	  sum += dy*dy/(2.0 * noise * noise);
	}
	break;
	
      case aemresponse::DIRECTION_Z:
	if (r.response.size() != response.SZ.size()) {
	  throw AEMEXCEPTION("Size mismatch in Z response (%d != %d)\n",
			     (int)r.response.size(),
			     (int)response.SZ.size());
	}
	for (int l = 0; l < (int)response.SZ.size(); l ++) {
	  
	  double dz = r.response[l] - response.SZ[l];
	  double noise = h->noise(r.response[l], time[l], lambda_scale);
	  
	  sum += dz*dz/(2.0 * noise * noise);
	}
	break;
	    
      default:
	throw AEMEXCEPTION("Unhandled direction\n");
      }
      
    }

  }
  
  return sum;
}

void
GlobalPixel::initialize_mpi(MPI_Comm _communicator)
{
  communicator = _communicator;

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw AEMEXCEPTION("MPI Failure\n");
  }
  
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw AEMEXCEPTION("MPI Failure\n");
  }
}

double
GlobalPixel::likelihood_mpi()
{
  if (communicator == MPI_COMM_NULL || mpi_rank < 0 || mpi_size < 0) {
    throw AEMEXCEPTION("MPI Parameters unset\n");
  }
  
  cEarth1D earth1d;
  
  //
  // Setup layer thicknesses
  //
  earth1d.conductivity.resize(image->rows);
  earth1d.thickness.resize(image->rows - 1);
  
  for (int i = 0; i < (image->rows - 1); i ++) {
    earth1d.thickness[i] = image->layer_thickness[i];
  }
  
  double sum = 0.0;
  
  for (int i = mpi_rank; i < image->columns; i += mpi_size) {
    
    aempoint &p = observations->points[i];
    
    //
    // Construct geometry
    //
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
    //
    // Copy image column to earth model, our model is in log of conductivity so here we use exp
    //
    for (int j = 0; j < image->rows; j ++) {
      earth1d.conductivity[j] = exp(image->conductivity[j * image->columns + i]);
    }
    
    for (int k = 0; k < (int)forwardmodel.size(); k ++) {

      cTDEmSystem *f = forwardmodel[k];
      hierarchicalmodel *h = lambda[k];
      double *time = forwardmodel_time[k];
      const aemresponse &r = p.responses[k];
      
      cTDEmResponse response;
	
      f->forwardmodel(geometry,
		      earth1d,
		      response);
      
      switch (r.d) {
      case aemresponse::DIRECTION_X:
	if (r.response.size() != response.SX.size()) {
	  throw AEMEXCEPTION("Size mismatch in X response\n");
	}
	for (int l = 0; l < (int)response.SX.size(); l ++) {
	  
	  double dx = r.response[l] - response.SX[l];
	  double noise = h->noise(r.response[l], time[l], lambda_scale);
	  
	  sum += dx*dx/(2.0 * noise * noise);
	}
	break;
	
      case aemresponse::DIRECTION_Y:
	if (r.response.size() != response.SY.size()) {
	  throw AEMEXCEPTION("Size mismatch in Y response\n");
	}
	for (int l = 0; l < (int)response.SY.size(); l ++) {
	  
	  double dy = r.response[l] - response.SY[l];
	  double noise = h->noise(r.response[l], time[l], lambda_scale);
	  
	  sum += dy*dy/(2.0 * noise * noise);
	}
	break;
	
      case aemresponse::DIRECTION_Z:
	if (r.response.size() != response.SZ.size()) {
	  throw AEMEXCEPTION("Size mismatch in Z response (%d != %d)\n",
			     (int)r.response.size(),
			     (int)response.SZ.size());
	}
	for (int l = 0; l < (int)response.SZ.size(); l ++) {
	  
	  double dz = r.response[l] - response.SZ[l];
	  double noise = h->noise(r.response[l], time[l], lambda_scale);
	  
	  sum += dz*dz/(2.0 * noise * noise);
	}
	break;
	    
      default:
	throw AEMEXCEPTION("Unhandled direction\n");
      }
      
    }

  }
  
  double total;
  
  if (MPI_Reduce(&sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
    throw AEMEXCEPTION("Likelihood failed in reducing\n");
  }
  if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
    throw AEMEXCEPTION("Likelihood failed in broadcast\n");
  }
    
  return total;
}

void
GlobalPixel::load_initial_model(const char *filename)
{
  FILE *fp = fopen(filename, "r");

  if (image == nullptr) {
    throw AEMEXCEPTION("Image null\n");
  }
  
  if (fp == NULL) {
    throw AEMEXCEPTION("Failed to load initial model\n");
  }

  for (int j = 0; j < image->rows; j ++) {
    for (int i = 0; i < image->columns; i ++) {

      double v;
      if (fscanf(fp, "%lf ", &v) != 1) {
	throw AEMEXCEPTION("Failed to read image\n");
      }

      image->conductivity[j * image->columns + i] = log(v);
    }
  }

  fclose(fp);
}
