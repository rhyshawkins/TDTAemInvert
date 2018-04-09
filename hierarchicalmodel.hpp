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
#ifndef hierarchicalmodel_hpp
#define hierarchicalmodel_hpp

#include <stdio.h>
#include <map>
#include <string>
#include <vector>

class hierarchicalmodel {
public:

  hierarchicalmodel();
  virtual ~hierarchicalmodel();

  virtual int nparameters() const = 0;

  virtual double getparameter(int i) const = 0;

  virtual void setparameter(int i, double v) = 0;
  
  virtual double noise(double observed_magnitude,
		       double observed_time,
		       double scale) = 0;

  virtual double nll(const std::vector<double> &observed_response,
		     const double *time,
		     const double *residuals,
		     double lambda_scale,
		     double *residuals_normed,
		     double &log_normalization) = 0;
  
  static hierarchicalmodel *load(const char *filename);

  typedef hierarchicalmodel* (*reader_function_t)(FILE *fp);

private:
  
  static std::map<std::string, reader_function_t> readers;
  
};

class independentgaussianhierarchicalmodel : public hierarchicalmodel {
public:
  
  independentgaussianhierarchicalmodel();
  virtual ~independentgaussianhierarchicalmodel();

  virtual int nparameters() const;

  virtual double getparameter(int i) const;

  virtual void setparameter(int i, double v);
  
  virtual double noise(double observed_magnitude,
		       double observed_time,
		       double scale);

  virtual double nll(const std::vector<double> &observed_response,
		     const double *time,
		     const double *residuals,
		     double lambda_scale,
		     double *residuals_normed,
		     double &log_normalization);  

  static hierarchicalmodel *read(FILE *fp);

private:

  double sigma;
  
};

class hyperbolichierarchicalmodel : public hierarchicalmodel {
public:
  
  hyperbolichierarchicalmodel();
  virtual ~hyperbolichierarchicalmodel();

  virtual int nparameters() const;

  virtual double getparameter(int i) const;

  virtual void setparameter(int i, double v);
  
  virtual double noise(double observed_magnitude,
		       double observed_time,
		       double scale);

  virtual double nll(const std::vector<double> &observed_response,
		     const double *time,
		     const double *residuals,
		     double lambda_scale,
		     double *residuals_normed,
		     double &log_normalization);  

  static hierarchicalmodel *read(FILE *fp);

private:

  double A, B, C;

};

class brodiehierarchicalmodel : public hierarchicalmodel {
public:

  brodiehierarchicalmodel();
  virtual ~brodiehierarchicalmodel();

  virtual int nparameters() const;

  virtual double getparameter(int i) const;

  virtual void setparameter(int i, double v);

  virtual double noise(double observed_magnitude,
		       double observed_time,
		       double scale);

  virtual double nll(const std::vector<double> &observed_response,
		     const double *time,
		     const double *residuals,
		     double lambda_scale,
		     double *residuals_normed,
		     double &log_normalization);  

  static hierarchicalmodel *read(FILE *fp);

private:

  double additive_noise(double observed_time);
  
  int ntimes;
  double *time;
  double *additive;

  double relative;

};

class covariancehierarchicalmodel : public hierarchicalmodel {
public:

  covariancehierarchicalmodel();
  virtual ~covariancehierarchicalmodel();

  virtual int nparameters() const;

  virtual double getparameter(int i) const;

  virtual void setparameter(int i, double v);

  virtual double noise(double observed_magnitude,
		       double observed_time,
		       double scale);

  virtual double nll(const std::vector<double> &observed_response,
		     const double *time,
		     const double *residuals,
		     double lambda_scale,
		     double *residuals_normed,
		     double &log_normalization);  

  static hierarchicalmodel *read(FILE *fp);

private:

  int size;
  double *w;
  double *v;
  
};
    
#endif // hierarchicalmodel.hpp
  
