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

#include "global.hpp"

extern "C" {
  #include "hnk_cartesian_nonsquare.h"

  #include "slog.h"
};

#include "constants.hpp"

const int CHAIN_STEPS = 1000000;

int global_coordtoindex(void *user, int i, int j, int k, int depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  return wavetree2d_sub_from_2dindices(wt, i, j);
}

int global_indextocoord(void *user, int index, int *i, int *j, int *k, int *depth)
{
  wavetree2d_sub_t *wt = (wavetree2d_sub_t*)user;

  if (wavetree2d_sub_2dindices(wt, index, i, j) < 0) {
    return -1;
  }

  *k = 0;
  *depth = wavetree2d_sub_depthofindex(wt, index);
  return 0;
}

Global::Global(const char *filename,
	       const std::vector<std::string> &stm_files,
	       const char *initial_model,
	       const char *prior_file,
	       int _degreex,
	       int _degreey,
	       double _depth,
	       const std::vector<std::string> &hierarchical_files,
	       int seed,
	       int _kmax,
	       bool _posteriork,
	       int hwavelet,
	       int vwavelet) :
  kmax(_kmax),
  treemaxdepth(-1),
  depth(_depth),
  wt(nullptr),
  ch(nullptr),
  hnk(nullptr),
  proposal(nullptr),
  degreex(_degreex),
  degreey(_degreey),
  observations(nullptr),
  image(nullptr),
  model(nullptr),
  workspace(nullptr),
  mean_residual_n(0),
  residual(nullptr),
  mean_residual(nullptr),
  last_valid_residual(nullptr),
  residual_normed(nullptr),
  mean_residual_normed(nullptr),
  last_valid_residual_normed(nullptr),
  residuals_valid(false),
  residual_hist_bins(100),
  residual_hist_min(-5.0),
  residual_hist_max(5.0),
  residual_hist(nullptr),
  width(-1),
  height(-1),
  size(-1),
  ncoeff(-1),
  lambda_scale(1.0),
  current_likelihood(-1.0),
  coeff_hist(nullptr),
  random(seed),
  posteriork(_posteriork),
  communicator(MPI_COMM_NULL),
  mpi_size(-1),
  mpi_rank(-1),
  temperature(1.0),
  column_offsets(nullptr),
  column_sizes(nullptr),
  residual_offsets(nullptr),
  residual_sizes(nullptr)
{
  if (degreex < 0 || degreex >= 16 ||
      degreey < 0 || degreey >= 16) {
    throw AEMEXCEPTION("Degree(s) out of range: %d x %d\n", degreex, degreey);
  }
  
  if (depth <= 0.0) {
    throw AEMEXCEPTION("Depth out of range\n");
  }

  if (!posteriork) {
    //
    // Load observations
    //
    observations = new aemobservations(filename);

    //
    // Load stm files
    //
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

      //
      // Create covariance information
      //
      cov_count.push_back(p->WinSpec.size());
      cov_delta.push_back(new double[p->WinSpec.size()]);
      cov_mu.push_back(new double[p->WinSpec.size()]);
      cov_sigma.push_back(new double[p->WinSpec.size() * p->WinSpec.size()]);
      
    }

    for (auto &h : hierarchical_files) {
      hierarchicalmodel *m = hierarchicalmodel::load(h.c_str());
      if (m == nullptr) {
	throw AEMEXCEPTION("Failed to create/load hierarchical model");
      }
      
      lambda.push_back(m);
    }
    
    if (forwardmodel.size() != observations->points[0].responses.size()) {
      throw AEMEXCEPTION("Mismatch in STM and responses size: %d != %d\n",
			 (int)forwardmodel.size(),
			 (int)observations->points[0].responses.size());
    }
  }

  wt = wavetree2d_sub_create(degreex, degreey, 0.0);
  if (wt == NULL) {
    throw AEMEXCEPTION("Failed to create wavetree\n");
  }

  width = wavetree2d_sub_get_width(wt);
  height = wavetree2d_sub_get_height(wt);
  size = wavetree2d_sub_get_size(wt);
  ncoeff = wavetree2d_sub_get_ncoeff(wt);
  treemaxdepth = wavetree2d_sub_maxdepth(wt);

  INFO("Image: %d x %d\n", width, height);

  if (!posteriork) {
    if ((int)observations->points.size() != width) {
      throw AEMEXCEPTION("Image size mismatch to observations: %d != %d\n",
    			 width,
    			 (int)observations->points.size());
    }
    
    image = new aemimage(height, width, depth, DEFAULT_CONDUCTIVITY);

    model = new double[size];
    int workspacesize = width;
    if (height > width) {
      workspacesize = height;
    }
    workspace = new double[workspacesize];

    int ntotal = observations->total_response_datapoints();
    INFO("Data: %d total points\n", ntotal);

    residual_size = ntotal;
    residual = new double[residual_size];
    mean_residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];
    last_valid_residual_normed = new double [residual_size];

    residuals_per_column = residual_size/image->columns;

    residual_hist = new int[residual_size * residual_hist_bins];

    reset_residuals();
  }
  
  if (initial_model == NULL) {
    if (wavetree2d_sub_initialize(wt, log(DEFAULT_CONDUCTIVITY)) < 0) {
      throw AEMEXCEPTION("Failed to initialize wavetree\n");
    }
  } else {

    if (wavetree2d_sub_load_promote(wt, initial_model) < 0) {
      throw AEMEXCEPTION("Failed to load initial model: %s\n", initial_model);
    }

    INFO("Loaded model with %d coefficients\n", wavetree2d_sub_coeff_count(wt));
  }

  //
  // Hnk Ratio
  //
  if (kmax > ncoeff) {
    INFO("Warning: kmax truncated to %d\n", ncoeff);
    kmax = ncoeff;
  }
  
  hnk = hnk_cartesian_nonsquare_2D_create_sub(degreex,
					      degreey,
					      kmax);
  if (hnk == NULL) {
    throw AEMEXCEPTION("Failed to create hnk table\n");
  }

  //
  // Chain History
  //
  ch = chain_history_create(CHAIN_STEPS);
  if (ch == nullptr) {
    throw AEMEXCEPTION("Failed to create chain history\n");
  }

  
  //
  // Initialse coeff histogram
  //
  coeff_hist = coefficient_histogram_create(ncoeff, 100, -1.0, 1.0,
					    global_coordtoindex,
					    global_indextocoord,
					    wt);
  if (coeff_hist == NULL) {
    throw AEMEXCEPTION("Failed to create coefficient histogram\n");
  }

  //
  // Create proposal structure.
  //
  if (prior_file != nullptr) {
    proposal = wavetree_pp_load(prior_file, seed, coeff_hist);
    if (proposal == NULL) {
      throw AEMEXCEPTION("Failed to load proposal file\n");
    }
    
    for (int i = 0; i < ncoeff; i ++) {
      int depth = wavetree2d_sub_depthofindex(wt, i);
      
      int ii, ij;
      if (wavetree2d_sub_2dindices(wt, i, &ii, &ij) < 0) {
	throw AEMEXCEPTION("Failed to get 2d indices\n");
      }
      
      double vmin, vmax;
      if (wavetree_pp_prior_range2d(proposal,
				    ii,
				    ij,
				    depth,
				    treemaxdepth,
				    0.0,
				    &vmin,
				    &vmax) < 0) {
	throw AEMEXCEPTION("Failed to get coefficient range\n");
      }
      
      if (coefficient_histogram_set_range(coeff_hist, 
					  i,
					  vmin,
					  vmax) < 0) {
	throw AEMEXCEPTION("Failed to set coefficient histogram range\n");
      }
    }
  }
    
  hwaveletf = wavelet_inverse_function_from_id(hwavelet);
  if (hwaveletf == nullptr) {
    throw AEMEXCEPTION("Invalid horizontal wavelet %d\n", hwavelet);
  }

  vwaveletf = wavelet_inverse_function_from_id(vwavelet);
  if (vwaveletf == nullptr) {
    throw AEMEXCEPTION("Invalid vertical wavelet %d\n", vwavelet);
  }

}

Global::~Global()
{
}

double
Global::likelihood(double &log_normalization)
{
  if (!posteriork) {
    
    cEarth1D earth1d;

    //
    // Setup layer thicknesses
    //
    earth1d.conductivity.resize(image->rows);
    earth1d.thickness.resize(image->rows - 1);
    
    for (int i = 0; i < (image->rows - 1); i ++) {
      earth1d.thickness[i] = image->layer_thickness[i];
    }
    
    //
    // Get tree model wavelet coefficients
    //
    memset(image->conductivity, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, image->conductivity, size) < 0) {
      throw AEMEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(image->conductivity,
			       width,
			       height,
			       width,
			       workspace,
			       hwaveletf,
			       vwaveletf,
			       1) < 0) {
      throw AEMEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    double sum = 0.0;
    int residual_offset;

    log_normalization = 0.0;
    
    for (int i = 0; i < image->columns; i ++) {

      residual_offset = i * residuals_per_column;
      aempoint &p = observations->points[i];
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
      //
      // Copy image column to earth model, our model is in log of conductivity so here we use exp
      //
      for (int j = 0; j < image->rows; j ++) {
	earth1d.conductivity[j] = exp(image->conductivity[j * image->columns + i]);
      }
      
      double point_sum = 0.0;
      
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
	    throw AEMEXCEPTION("Size mismatch in X response (%d != %d)\n", (int)r.response.size(), (int)response.SX.size());
	  }
	  for (int l = 0; l < (int)response.SX.size(); l ++) {
	    residual[residual_offset + l] = r.response[l] - response.SX[l];
	  }
	  point_sum +=
	    h->nll(r.response,
		   time,
		   residual + residual_offset,
		   lambda_scale,
		   residual_normed + residual_offset,
		   log_normalization);

	  residual_offset += response.SX.size();
	  break;
	  
	case aemresponse::DIRECTION_Y:
	  if (r.response.size() != response.SY.size()) {
	    throw AEMEXCEPTION("Size mismatch in Y response (%d != %d)\n", (int)r.response.size(), (int)response.SY.size());
	  }

	  for (int l = 0; l < (int)response.SY.size(); l ++) {
	    residual[residual_offset + l] = r.response[l] - response.SY[l];
	  }
	  point_sum +=
	    h->nll(r.response,
		   time,
		   residual + residual_offset,
		   lambda_scale,
		   residual_normed + residual_offset,
		   log_normalization);

	  residual_offset += response.SY.size();
	  break;
	  
	case aemresponse::DIRECTION_Z:
	  if (r.response.size() != response.SZ.size()) {
	    throw AEMEXCEPTION("Size mismatch in Z response (%d != %d)\n",
			       (int)r.response.size(),
			       (int)response.SZ.size());
	  }
	  for (int l = 0; l < (int)response.SZ.size(); l ++) {
	    residual[residual_offset + l] = r.response[l] - response.SZ[l];
	  }
	  point_sum +=
	    h->nll(r.response,
		   time,
		   residual + residual_offset,
		   lambda_scale,
		   residual_normed + residual_offset,
		   log_normalization);

	  residual_offset += response.SZ.size();
	  break;
	  
	default:
	  throw AEMEXCEPTION("Unhandled direction\n");
	}
	
	
      }

      sum += point_sum;
    }
    
    return sum;
  } else {
    return 1.0;
  }
}

double
Global::hierarchical_likelihood(double proposed_lambda_scale,
				double &log_normalization)
{
  log_normalization = 0.0;

  if (!posteriork) {

    if (!residuals_valid) {
      double x;
      (void)likelihood(x);
      accept();
    }

    double sum = 0.0;
    int residual_offset;
  
    for (int i = 0; i < image->columns; i ++) {
      
      residual_offset = i * residuals_per_column;
      aempoint &p = observations->points[i];
      
      for (int k = 0; k < (int)forwardmodel.size(); k ++) {
	
	hierarchicalmodel *h = lambda[k];
	double *time = forwardmodel_time[k];
	const aemresponse &r = p.responses[k];
	
	sum += h->nll(r.response,
		      time,
		      last_valid_residual + residual_offset,
		      proposed_lambda_scale,
		      residual_normed + residual_offset,
		      log_normalization);

	residual_offset += r.response.size();
      }
    }

    return sum;

  } else {
    return 1.0;
  }
}

void
Global::initialize_mpi(MPI_Comm _communicator, double _temperature)
{
  communicator = _communicator;

  if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS) {
    throw AEMEXCEPTION("MPI Failure\n");
  }
  
  if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS) {
    throw AEMEXCEPTION("MPI Failure\n");
  }

  column_offsets = new int[mpi_size];
  column_sizes = new int[mpi_size];
  residual_offsets = new int[mpi_size];
  residual_sizes = new int[mpi_size];

  int columns = image->columns;
  int processes = mpi_size;

  //
  // Evenly distribute columns
  //
  for (int i = 0; i < mpi_size; i ++) {
    column_sizes[i] = columns/processes;
    residual_sizes[i] = column_sizes[i] * residuals_per_column;

    columns -= column_sizes[i];
    processes --;
  }

  column_offsets[0] = 0;
  residual_offsets[0] = 0;
  for (int i = 1; i < mpi_size; i ++) {
    column_offsets[i] = column_offsets[i - 1] + column_sizes[i - 1];
    residual_offsets[i] = column_offsets[i] * residuals_per_column;
    INFO("Split: %4d %4d", column_offsets[i], column_sizes[i]);
  }

  if (column_offsets[mpi_size - 1] + column_sizes[mpi_size - 1] != image->columns) {
    throw AEMEXCEPTION("Column sharing intialization failure\n");
  }

  temperature = _temperature;
}

double
Global::likelihood_mpi(double &log_normalization)
{
  if (communicator == MPI_COMM_NULL || mpi_rank < 0 || mpi_size < 0) {
    throw AEMEXCEPTION("MPI Parameters unset\n");
  }
  
  if (!posteriork) {
    
    cEarth1D earth1d;

    //
    // Setup layer thicknesses
    //
    earth1d.conductivity.resize(image->rows);
    earth1d.thickness.resize(image->rows - 1);
    
    for (int i = 0; i < (image->rows - 1); i ++) {
      earth1d.thickness[i] = image->layer_thickness[i];
    }

    //
    // Get tree model wavelet coefficients
    //
    memset(image->conductivity, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, image->conductivity, size) < 0) {
      throw AEMEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(image->conductivity,
			       width,
			       height,
			       width,
			       workspace,
			       hwaveletf,
			       vwaveletf,
			       1) < 0) {
      throw AEMEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    double sum = 0.0;
    double local_log_normalization = 0.0;
    int residual_offset;
    
    for (int mi = 0, i = column_offsets[mpi_rank]; mi < column_sizes[mpi_rank]; mi ++, i ++) {

      residual_offset = i * residuals_per_column;
      
      aempoint &p = observations->points[i];

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
      //
      // Copy image column to earth model, our model is in log of conductivity so here we use exp
      //
      for (int j = 0; j < image->rows; j ++) {
	earth1d.conductivity[j] = exp(image->conductivity[j * image->columns + i]);
      }
      
      double point_sum = 0.0;
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
	    residual[residual_offset + l] = r.response[l] - response.SX[l];
	  }
	  point_sum +=
	    h->nll(r.response,
		   time,
		   residual + residual_offset,
		   lambda_scale,
		   residual_normed + residual_offset,
		   local_log_normalization);

	  residual_offset += response.SX.size();
	  break;
	  
	case aemresponse::DIRECTION_Y:
	  if (r.response.size() != response.SY.size()) {
	    throw AEMEXCEPTION("Size mismatch in Y response\n");
	  }
	  for (int l = 0; l < (int)response.SY.size(); l ++) {
	    residual[residual_offset + l] = r.response[l] - response.SY[l];
	  }
	  point_sum +=
	    h->nll(r.response,
		   time,
		   residual + residual_offset,
		   lambda_scale,
		   residual_normed + residual_offset,
		   local_log_normalization);

	  residual_offset += response.SY.size();
	  break;
	  
	case aemresponse::DIRECTION_Z:
	  if (r.response.size() != response.SZ.size()) {
	    throw AEMEXCEPTION("Size mismatch in Z response (%d != %d)\n", (int)r.response.size(), (int)response.SZ.size());
	  }
	  for (int l = 0; l < (int)response.SY.size(); l ++) {
	    residual[residual_offset + l] = r.response[l] - response.SZ[l];
	  }
	  point_sum +=
	    h->nll(r.response,
		   time,
		   residual + residual_offset,
		   lambda_scale,
		   residual_normed + residual_offset,
		   local_log_normalization);

	  residual_offset += response.SZ.size();
	  break;
	  
	default:
	  throw AEMEXCEPTION("Unhandled direction\n");
	}
	
	
      }
      
      sum += point_sum;
	
    }

    double total;
    
    if (MPI_Reduce(&local_log_normalization, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
      throw AEMEXCEPTION("Likelihood failed in reducing\n");
    }
    if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      throw AEMEXCEPTION("Likelihood failed in broadcast\n");
    }

    log_normalization = total;
    
    if (MPI_Reduce(&sum, &total, 1, MPI_DOUBLE, MPI_SUM, 0, communicator) != MPI_SUCCESS) {
      throw AEMEXCEPTION("Likelihood failed in reducing\n");
    }
    if (MPI_Bcast(&total, 1, MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      throw AEMEXCEPTION("Likelihood failed in broadcast\n");
    }


    MPI_Allgatherv(residual + residual_offsets[mpi_rank],
		   residual_sizes[mpi_rank],
		   MPI_DOUBLE,
		   residual,
		   residual_sizes,
		   residual_offsets,
		   MPI_DOUBLE,
		   communicator);

    MPI_Allgatherv(residual_normed + residual_offsets[mpi_rank],
		   residual_sizes[mpi_rank],
		   MPI_DOUBLE,
		   residual_normed,
		   residual_sizes,
		   residual_offsets,
		   MPI_DOUBLE,
		   communicator);

    return total;
    
  } else {
    log_normalization = 0.0;
    
    return 1.0;
  }
}

double
Global::hierarchical_likelihood_mpi(double proposed_lambda_scale,
				    double &log_normalization)
{
  double like;

  log_normalization = 0.0;

  if (!residuals_valid) {
    double x;
    like = likelihood_mpi(x);
    accept();
  }
  
  //
  // 
  //
  // Since this is just a sum over residuals, just compute on root node and broadcast.
  //
  if (mpi_rank == 0) {
    like = hierarchical_likelihood(proposed_lambda_scale, log_normalization);
  }
  
  MPI_Bcast(&like, 1, MPI_DOUBLE, 0, communicator);
  MPI_Bcast(&log_normalization, 1, MPI_DOUBLE, 0, communicator);

  return like;
}

void
Global::reset_residuals()
{
  mean_residual_n = 0;
  for (int i = 0; i < residual_size; i ++) {
    residual[i] = 0.0;
    mean_residual[i] = 0.0;
    last_valid_residual[i] = 0.0;
      
    residual_normed[i] = 0.0;
    mean_residual_normed[i] = 0.0;
    last_valid_residual_normed[i] = 0.0;

    for (int j = 0; j < residual_hist_bins; j ++) {
      residual_hist[i * residual_hist_bins + j] = 0;
    }
  }

  cov_n = 0;

  for (int i = 0; i < (int)cov_count.size(); i ++) {
    int N = cov_count[i];

    for (int j = 0; j < N; j ++) {
      cov_delta[i][j] = 0.0;
      cov_mu[i][j] = 0.0;
    }

    N = N*N;
    for (int j = 0; j < N; j ++) {
      cov_sigma[i][j] = 0.0;
    }
  }
}

void
Global::invalidate_residuals()
{
  residuals_valid = false;
}

void
Global::accept()
{
  residuals_valid = true;
  if (!posteriork) {
    for (int i = 0; i < residual_size; i ++) {
      last_valid_residual[i] = residual[i];
      last_valid_residual_normed[i] = residual_normed[i];
    }

    update_residual_mean();
    update_residual_covariance();
  }
}

void
Global::accept_hierarchical()
{
}

void
Global::reject()
{
  update_residual_mean();
}

void
Global::reject_hierarchical()
{
}

void
Global::update_residual_mean()
{
  mean_residual_n ++;

  for (int i = 0; i < residual_size; i ++) {

    double delta = last_valid_residual[i] - mean_residual[i];
    mean_residual[i] += delta/(double)(mean_residual_n);

    delta = last_valid_residual_normed[i] - mean_residual_normed[i];
    mean_residual_normed[i] += delta/(double)(mean_residual_n);

    int hi = (int)((last_valid_residual_normed[i] - residual_hist_min)/(residual_hist_max - residual_hist_min) * (double)residual_hist_bins);
    if (hi >= 0 && hi < residual_hist_bins) {
      residual_hist[i * residual_hist_bins + hi] ++;
    }
  }
}

void
Global::update_residual_covariance()
{
  double *p = last_valid_residual;

  for (int k = 0; k < (int)(observations->points.size()); k ++) {

    cov_n ++;
    
    for (int i = 0; i < (int)cov_count.size(); i ++) {

      int N = cov_count[i];
      
      for (int j = 0; j < N; j ++) {
	cov_delta[i][j] = (p[j] - cov_mu[i][j])/(double)(cov_n);
	cov_mu[i][j] += cov_delta[i][j];
      }

      for (int j = 0; j < N; j ++) {
	for (int l = j; l < N; l ++) {

	  cov_sigma[i][j * N + l] +=
	    (double)(cov_n - 1)*cov_delta[i][j]*cov_delta[i][l] -
	    cov_sigma[i][j * N + l]/(double)(cov_n);
	}
      }

      p += N;
    }
  }
}


int
Global::get_residual_size() const
{
  return residual_size;
}

const double *
Global::get_mean_residuals() const
{
  return mean_residual;
}
  
const double *
Global::get_mean_normed_residuals() const
{
  return mean_residual_normed;
}

bool
Global::save_residual_histogram(const char *filename) const
{
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("Failed to create file");
    return false;
  }

  fprintf(fp, "%d %d %f %f\n", residual_size, residual_hist_bins, residual_hist_min, residual_hist_max);
  for (int i = 0; i < residual_size; i ++) {
    for (int j = 0; j < residual_hist_bins; j ++) {
      fprintf(fp, "%d ", residual_hist[i * residual_hist_bins + j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  return true;
}

bool
Global::save_residual_covariance(const char *filename) const
{
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    ERROR("Failed to create file\n");
    return false;
  }

  fprintf(fp, "%d\n", (int)cov_count.size());
  
  for (int i = 0; i < (int)cov_count.size(); i ++) {

    int N = cov_count[i];

    fprintf(fp, "%d\n", N);
    for (int j = 0; j < N; j ++) {
      fprintf(fp, "%.9g ", cov_mu[i][j]);
    }
    fprintf(fp, "\n");

    for (int j = 0; j < N; j ++) {
      for (int k = 0; k < N; k ++) {
	fprintf(fp, "%.9g ", cov_sigma[i][j * N + k]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);

  return true;
}

generic_lift_inverse1d_step_t
Global::wavelet_inverse_function_from_id(int id)
{
  switch (id) {
  case WAVELET_HAAR:
    return haar_lift_inverse1d_haar_step;

  case WAVELET_DAUB4:
    return daub4_dwt_inverse1d_daub4_step;

  case WAVELET_DAUB6:
    return daub6_dwt_inverse1d_daub6_step;

  case WAVELET_DAUB8:
    return daub8_dwt_inverse1d_daub8_step;

  case WAVELET_CDF97:
    return cdf97_lift_inverse1d_cdf97_step;

  case WAVELET_CDF97_PERIODIC:
    return cdf97_lift_periodic_inverse1d_cdf97_step;

  default:
    return nullptr;
  }
}

generic_lift_forward1d_step_t
Global::wavelet_forward_function_from_id(int id)
{
  switch (id) {
  case WAVELET_HAAR:
    return haar_lift_forward1d_haar_step;

  case WAVELET_DAUB4:
    return daub4_dwt_forward1d_daub4_step;

  case WAVELET_DAUB6:
    return daub6_dwt_forward1d_daub6_step;

  case WAVELET_DAUB8:
    return daub8_dwt_forward1d_daub8_step;

  case WAVELET_CDF97:
    return cdf97_lift_forward1d_cdf97_step;

  case WAVELET_CDF97_PERIODIC:
    return cdf97_lift_periodic_forward1d_cdf97_step;

  default:
    return nullptr;
  }
}
