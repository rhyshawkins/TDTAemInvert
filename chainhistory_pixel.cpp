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

#include <stdio.h>
#include <stdlib.h>

#include "chainhistory_pixel.hpp"

#include "aemimage.hpp"

PixelPerturbation::PixelPerturbation() :
  accepted(false),
  idx(-1),
  oldvalue(0.0),
  newvalue(0.0)
{
}

PixelPerturbation::PixelPerturbation(bool _accepted,
				     int _idx,
				     double _oldvalue,
				     double _newvalue) :
  accepted(_accepted),
  idx(_idx),
  oldvalue(_oldvalue),
  newvalue(_newvalue)
{
}

ChainHistoryPixel::ChainHistoryPixel(int _rows, int _columns) :
  rows(_rows),
  columns(_columns),
  initial_image(new double[rows * columns])
{
}


ChainHistoryPixel::ChainHistoryPixel(aemimage &initial_model) :
  rows(initial_model.rows),
  columns(initial_model.columns),
  initial_image(new double[initial_model.rows * initial_model.columns])
{
  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {
      initial_image[j * columns + i] = initial_model.conductivity[j * columns + i];
    }
  }
}
  
ChainHistoryPixel *
ChainHistoryPixel::load(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    return NULL;
  }

  int rows, columns;
  if (fscanf(fp, "%d %d\n", &rows, &columns) != 2) {
    return NULL;
  }

  ChainHistoryPixel *r = new ChainHistoryPixel(rows, columns);
  
  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {
      double v;
      if (fscanf(fp, "%lf\n", &v) != 1) {
	return NULL;
      }
      r->initial_image[j * columns + i] = v;
    }
  }

  while (!feof(fp)) {
    int accept;
    int idx;
    double oldv;
    double newv;
    
    if (fscanf(fp, "%d %d %lf %lf\n", &accept, &idx, &oldv, &newv) != 4) {
      if (feof(fp)) {
	break;
      } else {
	return NULL;
      }
    }

    r->history.push_back(PixelPerturbation((bool)accept, idx, oldv, newv));
  }

  fclose(fp);
  return r;
}

bool
ChainHistoryPixel::save(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create %s\n", filename);
    return false;
  }

  fprintf(fp, "%d %d\n", rows, columns);
  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      fprintf(fp, "%15.9f ", initial_image[j * columns + i]);
    }
    fprintf(fp, "\n");
  }

  for (auto &pp : history) {

    fprintf(fp, "%d %d %15.9f %15.9f\n", (int)pp.accepted, pp.idx, pp.oldvalue, pp.newvalue);

  }

  fclose(fp);
  return true;
}
  
