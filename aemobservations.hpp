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
#ifndef aemobservations_hpp
#define aemobservations_hpp

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "aemexception.hpp"

class aemresponse {
public:

  typedef enum {
    DIRECTION_X = 0,
    DIRECTION_Y = 1,
    DIRECTION_Z = 2
  } direction_t;

  aemresponse(direction_t _d = DIRECTION_X) :
    d(_d)
  {
  }

  bool write_text(FILE *fp) const
  {
    fprintf(fp, "%d %d ", (int)d, (int)response.size());
    for (auto &r : response) {
      fprintf(fp, "%.9g ", r);
    }
    return true;
  }

  bool read_text(FILE *fp)
  {
    int id, nr;
    if (fscanf(fp, "%d %d", &id, &nr) != 2) {
      return false;
    }
    if (id < DIRECTION_X || id > DIRECTION_Z) {
      return false;
    }
    
    d = (direction_t)id;

    for (int i = 0; i < nr; i ++) {
      double r;
      if (fscanf(fp, "%lf", &r) != 1) {
	fprintf(stderr, "aempoint::read_txt: failed to read value\n");
	return false;
      }

      response.push_back(r);
    }

    return true;
  }

  direction_t d;
  std::vector<double> response;
};

class aempoint {
public:

  aempoint() :
    tx_height(0.0),
    tx_roll(0.0),
    tx_pitch(0.0),
    tx_yaw(0.0),
    txrx_dx(0.0),
    txrx_dy(0.0),
    txrx_dz(0.0),
    rx_roll(0.0),
    rx_pitch(0.0),
    rx_yaw(0.0)
  {
    reset();

    cached_residual = -1.0;
  }

  aempoint(double height,
	   double roll,
	   double pitch,
	   double yaw,
	   double dx,
	   double dy,
	   double dz,
	   double rxroll,
	   double rxpitch,
	   double rxyaw) :
    tx_height(height),
    tx_roll(roll),
    tx_pitch(pitch),
    tx_yaw(yaw),
    txrx_dx(dx),
    txrx_dy(dy),
    txrx_dz(dz),
    rx_roll(rxroll),
    rx_pitch(rxpitch),
    rx_yaw(rxyaw)
  {

    cached_residual = -1.0;
  }

  void reset()
  {
    responses.clear();
  }

  bool write_text(FILE *fp) const
  {
    fprintf(fp,
	    "%15.9f "
	    "%15.9f %15.9f %15.9f "
	    "%15.9f %15.9f %15.9f "
	    "%15.9f %15.9f %15.9f "
	    "%d ",
	    tx_height,
	    tx_roll,
	    tx_pitch,
	    tx_yaw,
	    txrx_dx,
	    txrx_dy,
	    txrx_dz,
	    rx_roll,
	    rx_pitch,
	    rx_yaw,
	    (int)responses.size());

    for (auto &r : responses) {
      if (!r.write_text(fp)) {
	return false;
      }
    }
    
    fprintf(fp, "\n");

    return true;
  }

  bool read_text(FILE *fp)
  {
    int nresponse;

    if (fscanf(fp,
	       "%lf "
	       "%lf %lf %lf "
	       "%lf %lf %lf "
	       "%lf %lf %lf "
	       "%d ",
	       &tx_height,
	       &tx_roll, &tx_pitch, &tx_yaw,
	       &txrx_dx, &txrx_dy, &txrx_dz,
	       &rx_roll, &rx_pitch, &rx_yaw,
	       &nresponse) != 11) {
      return false;
    }

    for (int i = 0; i < nresponse; i ++) {
      aemresponse r;

      if (!r.read_text(fp)) {
	return false;
      }

      responses.push_back(r);
    }

    return true;
  }

  double tx_height;
  double tx_roll;
  double tx_pitch;
  double tx_yaw;
  double txrx_dx;
  double txrx_dy;
  double txrx_dz;
  double rx_roll;
  double rx_pitch;
  double rx_yaw;

  std::vector<aemresponse> responses;

  double cached_residual;

};

    
class aemobservations {
public:

  
  aemobservations()
  {
  }

  aemobservations(const char *filename)
  {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
      throw AEMEXCEPTION("Failed to open %s for reading\n", filename);
    }

    while (true) {
      aempoint p;

      if (!p.read_text(fp)) {
	if (feof(fp)) {
	  break;
	} else {
	  throw AEMEXCEPTION("Failed to read line from file\n");
	}
      }

      points.push_back(p);
    }
  }

  bool save(const char *filename) const
  {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
      return false;
    }

    for (auto &a : points) {
      if (!a.write_text(fp)) {
	return false;
      }
    }

    fclose(fp);
    return true;
  }

  int total_response_datapoints()
  {
    int c = 0;

    for (auto &p : points) {
      for (auto &r : p.responses) {

	c += r.response.size();

      }
    }

    return c;
  }
  
  std::vector<aempoint> points;
};

#endif // aemobservations_hpp
