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
#ifndef global_hpp
#define global_hpp

#include <vector>
#include <string>
#include <set>

#include <gmp.h>

#include <mpi.h>

extern "C" {
#include "hnk.h"
#include "wavetree2d_sub.h"
#include "wavetreepp.h"
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
#include "generic_lift.h"
};

#include "rng.hpp"
#include "aemimage.hpp"
#include "aemobservations.hpp"
#include "hierarchicalmodel.hpp"

#include "tdemsystem.h"
#include "general_types.h"

extern const double DEFAULT_CONDUCTIVITY;

class Global {
public:

  enum {
    WAVELET_HAAR = 0,
    WAVELET_DAUB4 = 1,
    WAVELET_DAUB6 = 2,
    WAVELET_DAUB8 = 3,
    WAVELET_CDF97 = 4,
    WAVELET_CDF97_PERIODIC = 5,
    WAVELET_MAX = 5
  };

  Global(const char *filename,
	 const std::vector<std::string> &stm_files,
	 const char *initial_model,
	 const char *prior_file,
	 int degreex,
	 int degreey,
	 double depth,
	 const std::vector<std::string> &hierarchical_files,
	 int seed,
	 int kmax,
	 bool posteriork,
	 int hwavelet,
	 int vwavelet);
  ~Global();

  double likelihood(double &log_normalization);

  double hierarchical_likelihood(double proposed_lambda_scale,
				 double &log_hierarchical_normalization);

  void initialize_mpi(MPI_Comm communicator, double temperature = 1.0);

  double likelihood_mpi(double &log_normalization);

  double hierarchical_likelihood_mpi(double proposed_lambda_scale,
				     double &log_hierarchical_normalization);

  void resample(MPI_Comm temperature_communicator, double resample_temperature);

  void reset_residuals();

  void invalidate_residuals();

  void accept();

  void accept_hierarchical();

  void reject();

  void reject_hierarchical();

  void update_residual_mean();

  void update_residual_covariance();

  int get_residual_size() const;
  
  const double *get_mean_residuals() const;

  const double *get_mean_normed_residuals() const;

  bool save_residual_histogram(const char *filename) const;
  bool save_residual_covariance(const char *filename) const;
  
  static generic_lift_inverse1d_step_t wavelet_inverse_function_from_id(int id);
  static generic_lift_forward1d_step_t wavelet_forward_function_from_id(int id);
  
  int kmax;
  int maxdepth;
  int treemaxdepth;

  double depth;

  wavetree2d_sub_t *wt;
  chain_history_t *ch;

  hnk_t *hnk;
  wavetree_pp_t *proposal;

  int degreex;
  int degreey;

  std::vector<cTDEmSystem*> forwardmodel;
  std::vector<double*> forwardmodel_time;
  
  aemobservations *observations;
  aemimage *image;

  double *model;
  double *workspace;

  int mean_residual_n;
  int residual_size;
  int residuals_per_column;
  
  double *residual;
  double *mean_residual;
  double *last_valid_residual;

  double *residual_normed;
  double *mean_residual_normed;
  double *last_valid_residual_normed;

  bool residuals_valid;

  int residual_hist_bins;
  double residual_hist_min;
  double residual_hist_max;
  int *residual_hist;

  int width;
  int height;
  int size;
  int ncoeff;

  std::vector<hierarchicalmodel*> lambda;
  double lambda_scale;
  
  double current_likelihood;
  double current_log_normalization;

  coefficient_histogram_t *coeff_hist;

  Rng random;
  bool posteriork;

  generic_lift_inverse1d_step_t hwaveletf;
  generic_lift_inverse1d_step_t vwaveletf;

  MPI_Comm communicator;
  int mpi_size;
  int mpi_rank;
  double temperature;

  int *column_offsets;
  int *column_sizes;
  int *residual_offsets;
  int *residual_sizes;

  int cov_n;
  std::vector<int> cov_count;
  std::vector<double*> cov_delta;
  std::vector<double*> cov_mu;
  std::vector<double*> cov_sigma;

};


#endif // global_hpp
