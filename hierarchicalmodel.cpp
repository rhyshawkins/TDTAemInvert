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

#include <math.h>
#include <string.h>

#include "aemexception.hpp"

#include "hierarchicalmodel.hpp"

extern "C" {
  #include "slog.h"
};

std::map<std::string, hierarchicalmodel::reader_function_t> hierarchicalmodel::readers =
  {
    {"IndependentGaussian", independentgaussianhierarchicalmodel::read},
    {"Hyperbolic", hyperbolichierarchicalmodel::read},
    {"Brodie", brodiehierarchicalmodel::read},
    {"Covariance", covariancehierarchicalmodel::read}
  };

hierarchicalmodel::hierarchicalmodel()
{
}

hierarchicalmodel::~hierarchicalmodel()
{
}

hierarchicalmodel *
hierarchicalmodel::load(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    ERROR("Failed to open file: %s", filename);
    return nullptr;
  }

  char modelname[1024];
  if (fgets(modelname, sizeof(modelname) - 1, fp) == NULL) {
    ERROR("Failed to read hierarchical model name");
    return nullptr;
  }

  // Strip trailing whitespace
  int i = strlen(modelname) - 1;
  while (i > 0 && isspace(modelname[i])) {
    modelname[i] = '\0';
    i --;
  }

  std::map<std::string, hierarchicalmodel::reader_function_t>::iterator r = readers.find(modelname);
  if (r == readers.end()) {
    ERROR("Unknown hierarchical model");
    fclose(fp);
    return nullptr;
  } else {
    hierarchicalmodel *m = (r->second)(fp);
    fclose(fp);

    return m;
  }
}

//
// Independent Gaussian
//
independentgaussianhierarchicalmodel::independentgaussianhierarchicalmodel() :
  sigma(1.0)
{
}

independentgaussianhierarchicalmodel::~independentgaussianhierarchicalmodel()
{
}

int
independentgaussianhierarchicalmodel::nparameters() const
{
  return 1;
}

double
independentgaussianhierarchicalmodel::getparameter(int i) const
{
  if (i == 0) {
    return sigma;
  } else {
    throw AEMEXCEPTION("Invalid index\n");
  }
}

void
independentgaussianhierarchicalmodel::setparameter(int i, double v)
{
  if (i == 0) {
    if (v <= 0.0) {
      throw AEMEXCEPTION("Sigma out of range\n");
    }

    sigma = v;
  } else {
    throw AEMEXCEPTION("Invalid index\n");
  }
}
  
double
independentgaussianhierarchicalmodel::noise(double observed_magnitude,
					    double observed_time,
					    double scale)
{
  return sigma * scale;
}

double
independentgaussianhierarchicalmodel::nll(const std::vector<double> &observed_response,
					  const double *time,
					  const double *residuals,
					  double lambda_scale,
					  double *residuals_normed,
					  double &log_normalization)
{
  int i = 0;
  double sum = 0.0;
  
  for (auto &r : observed_response) {

    double n = noise(r, time[i], lambda_scale);

    residuals_normed[i] = residuals[i]/n;

    sum += residuals_normed[i] * residuals_normed[i] * 0.5;

    log_normalization += log(n);
    i ++;
  }

  return sum;
}

hierarchicalmodel *
independentgaussianhierarchicalmodel::read(FILE *fp)
{
  double sigma;

  if (fscanf(fp, "%lf", &sigma) != 1) {
    ERROR("Failed to read sigma");
    return nullptr;
  }

  independentgaussianhierarchicalmodel *m = new independentgaussianhierarchicalmodel();
  m->sigma = sigma;

  return m;
}

hyperbolichierarchicalmodel::hyperbolichierarchicalmodel() :
  A(0.0),
  B(1.0),
  C(1.0)
{
}

hyperbolichierarchicalmodel::~hyperbolichierarchicalmodel()
{
}

int
hyperbolichierarchicalmodel::nparameters() const
{
  return 3;
}

double
hyperbolichierarchicalmodel::getparameter(int i) const
{
  switch (i) {
  case 0:
    return A;
  case 1:
    return B;
  case 2:
    return C;
  default:
    throw AEMEXCEPTION("Invalid index\n");
  }
}

void
hyperbolichierarchicalmodel::setparameter(int i, double v)
{
  switch(i) {
  case 0:
    if (v < 0.0) {
      throw AEMEXCEPTION("Invalid value for A\n");
    }
    A = v;
    break;

  case 1:
    if (v <= 0.0) {
      throw AEMEXCEPTION("Invalid value for B\n");
    }
    B = v;
    break;

  case 2:
    C = v;
    break;
    
  default:
    throw AEMEXCEPTION("Invalid index\n");
  }
}
  
double
hyperbolichierarchicalmodel::noise(double observed_magnitude,
				   double observed_time,
				   double scale)
{
  return A + (B * scale)/pow(observed_time, C);
}

double
hyperbolichierarchicalmodel::nll(const std::vector<double> &observed_response,
				 const double *time,
				 const double *residuals,
				 double lambda_scale,
				 double *residuals_normed,
				 double &log_normalization)
{
  int i = 0;
  double sum = 0.0;
  
  for (auto &r : observed_response) {

    double n = noise(r, time[i], lambda_scale);

    residuals_normed[i] = residuals[i]/n;

    sum += residuals_normed[i] * residuals_normed[i] * 0.5;

    log_normalization += log(n);
    i ++;
  }

  return sum;
}

hierarchicalmodel *
hyperbolichierarchicalmodel::read(FILE *fp)
{
  double A, B, C;

  if (fscanf(fp, "%lf %lf %lf", &A, &B, &C) != 3) {
    ERROR("Failed to read hyperbolic parameters");
    return nullptr;
  }

  hyperbolichierarchicalmodel *m = new hyperbolichierarchicalmodel();
  m->A = A;
  m->B = B;
  m->C = C;

  return m;
}

//
// Ross Brodie's noise model
//

brodiehierarchicalmodel::brodiehierarchicalmodel() :
  ntimes(-1),
  time(nullptr),
  additive(nullptr),
  relative(0.0)
{
}

brodiehierarchicalmodel::~brodiehierarchicalmodel()
{
  delete [] time;
  delete [] additive;
}

int
brodiehierarchicalmodel::nparameters() const
{
  return 1;
}

double
brodiehierarchicalmodel::getparameter(int i) const
{
  if (i == 0) {
    return relative;
  } else {
    throw AEMEXCEPTION("Invalid index");
  }
}

void
brodiehierarchicalmodel::setparameter(int i, double v)
{
  if (i == 0 && v >= 0.0) {
    relative = v;
  } else {
    throw AEMEXCEPTION("Invalid values");
  }
}

double
brodiehierarchicalmodel::noise(double observed_magnitude,
			       double observed_time,
			       double scale)
{
  //
  // Copy/Paste from Ross's email
  //
  // We use a noise model as follows
  //
  // N(i) = sqrt((0.01*rn*d(i))^2 + an(i)^2)
  //
  // Where
  //   d(i) is the observed datum
  //   rn is the relative noise in %
  // an(i) is the additive noise floor in the same units as the data (in this case /A.m^4)
  //
  // the (guesstimated) relative and additive noises are shown below for each moment.
  //
  // Low moment data
  // RelativeNoise   = 3.6%
  // AdditiveNoise   = 5.7761e-011 7.7154e-012 5.7849e-012 3.9164e-012 3.1502e-012 2.5105e-012
  // 2.2912e-012  1.921e-012  1.733e-012  1.529e-012 1.2258e-012 9.6876e-013 9.0323e-013
  // 8.2181e-013 7.4835e-013 6.2648e-013 6.2901e-013 5.7157e-013 5.1475e-013
  //
  // High moment data
  // RelativeNoise    = 3.6%
  //
  // AdditiveNoise    = 2.5545e-013 2.0815e-013 1.9144e-013  1.592e-013 1.4598e-013 1.3402e-013
  // 1.2712e-013 1.0844e-013 1.0214e-013 9.7184e-014 9.0881e-014 8.4579e-014 7.7776e-014 6.9864e-014
  // 6.6747e-014 5.9365e-014   5.33e-014  4.843e-014 4.2199e-014 3.7096e-014  3.571e-014

  double rnoise = scale * relative * observed_magnitude;
  double anoise = additive_noise(observed_time);

  return sqrt(rnoise*rnoise + anoise*anoise);
}

double
brodiehierarchicalmodel::nll(const std::vector<double> &observed_response,
			     const double *time,
			     const double *residuals,
			     double lambda_scale,
			     double *residuals_normed,
			     double &log_normalization)
{
  int i = 0;
  double sum = 0.0;
  
  for (auto &r : observed_response) {

    double n = noise(r, time[i], lambda_scale);

    residuals_normed[i] = residuals[i]/n;

    sum += residuals_normed[i] * residuals_normed[i] * 0.5;

    log_normalization += log(n);
    i ++;
  }

  return sum;
}

double
brodiehierarchicalmodel::additive_noise(double observed_time)
{
  if (ntimes < 0 ||
      time == nullptr ||
      additive == nullptr) {
    throw AEMEXCEPTION("Uninitialized additive noise");
  }

  if (observed_time <= time[0]) {
    return additive[0];
  }

  if (observed_time >= time[ntimes - 1]) {
    return additive[ntimes - 1];
  }

  int i = 0;
  while (time[i + 1] < observed_time) {
    i ++;
  }

  // Return linear interpolation
  double alpha = (observed_time - time[i])/(time[i + 1] - time[i]);
  return (1.0 - alpha)*additive[i] + alpha*additive[i + 1];
}

hierarchicalmodel *
brodiehierarchicalmodel::read(FILE *fp)
{
  int ntimes;
  double relative;

  if (fscanf(fp, "%d %lf", &ntimes, &relative) != 2) {
    ERROR("Failed to parameters");
    return nullptr;
  }

  brodiehierarchicalmodel *m = new brodiehierarchicalmodel();
  m->ntimes = ntimes;
  m->relative = relative;

  m->time = new double[ntimes];
  m->additive = new double[ntimes];

  for (int i = 0; i < ntimes; i ++) {
    if (fscanf(fp, "%lf %lf", &(m->time[i]), &(m->additive[i])) != 2) {
      ERROR("Failed to read additive noise entry");
      return nullptr;
    }
  }

  return m;
}


//
// Covariance model
//

covariancehierarchicalmodel::covariancehierarchicalmodel() :
  size(-1),
  w(nullptr),
  v(nullptr)
{
}

covariancehierarchicalmodel::~covariancehierarchicalmodel()
{
}

int
covariancehierarchicalmodel::nparameters() const
{
  return 0;
}

double
covariancehierarchicalmodel::getparameter(int i) const
{
  throw AEMEXCEPTION("Invalid index\n");
}

void
covariancehierarchicalmodel::setparameter(int i, double v)
{
  throw AEMEXCEPTION("Invalid index\n");
}

double
covariancehierarchicalmodel::noise(double observed_magnitude,
				   double observed_time,
				   double scale)
{
  throw AEMEXCEPTION("Unsupported");
}

double
covariancehierarchicalmodel::nll(const std::vector<double> &observed_response,
				 const double *time,
				 const double *residuals,
				 double lambda_scale,
				 double *residuals_normed,
				 double &log_normalization)
{
  double sum = 0.0;
  double norm = 0.0;

  if (size != (int)observed_response.size()) {
    throw AEMEXCEPTION("Size mismatch");
  }

  for (int i = 0; i < size; i ++) {
    residuals_normed[i] = 0.0;

    for (int j = 0; j < size; j ++) {
      residuals_normed[i] += (residuals[j] * v[j * size + i])/sqrt(lambda_scale * w[i]);
    }

    norm += 0.5 * log(lambda_scale * w[i]);

    sum += residuals_normed[i] * residuals_normed[i] * 0.5;
	   
  }

  log_normalization += norm;

  return sum;
}

hierarchicalmodel *
covariancehierarchicalmodel::read(FILE *fp)
{
  int size;
  double scale;

  if (fscanf(fp, "%d", &size) != 1) {
    ERROR("Failed to parameters");
    return nullptr;
  }

  if (fscanf(fp, "%lf", &scale) != 1) {
    ERROR("Failed to read scale\n");
    return nullptr;
  }

  covariancehierarchicalmodel *m = new covariancehierarchicalmodel();

  m->size = size;
  m->w = new double[size];
  m->v = new double[size * size];

  for (int i = 0; i < size; i ++) {
    if (fscanf(fp, "%lf", &(m->w[i])) != 1) {
      ERROR("Failed to read eigen value");
      return nullptr;
    }

    if (m->w[i] != 1.0) {
      m->w[i] *= scale;
    }
  }

  for (int i = 0; i < (size * size); i ++) {
    if (fscanf(fp, "%lf", &(m->v[i])) != 1) {
      ERROR("Failed to read eigen vector entry");
      return nullptr;
    }
  }

  return m;
}
