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
#include <math.h>

#include "aemimage.hpp"

#include "aemexception.hpp"

#include "logspace.hpp"


aemimage::aemimage() :
  rows(-1),
  columns(-1),
  depth(0.0),
  conductivity(nullptr)
{
}

aemimage::aemimage(int _rows, int _columns, double _depth, double const_conductivity) :
  rows(_rows),
  columns(_columns),
  depth(_depth),
  conductivity(new double[rows * columns])
{
  int s = rows * columns;
  for (int i = 0; i < s; i ++) {
    conductivity[i] = const_conductivity;
  }

  update_layer_thickness();
}

aemimage::~aemimage()
{
  delete [] conductivity;
}

bool
aemimage::load(const char *filename)
{
  int nrows, ncols;
  double ndepth;
  
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    return false;
  }

  if (fscanf(fp, "%d %d %lf\n", &nrows, &ncols, &ndepth) != 3) {
    return false;
  }
  
  if (conductivity != nullptr) {
    delete [] conductivity;
  }

  rows = nrows;
  columns = ncols;
  depth = ndepth;

  conductivity = new double[rows * columns];

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      if (fscanf(fp, "%lf", &conductivity[j * columns + i]) != 1) {
	return false;
      }
    }
  }

  update_layer_thickness();
  
  fclose(fp);
  return true;
}

bool
aemimage::save(const char *filename)
{
  if (conductivity == nullptr) {
    return false;
  }

  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  fprintf(fp, "%d %d %15.9f\n", rows, columns, depth);

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      fprintf(fp, "%15.9f ", conductivity[j * columns + i]);

    }

    fprintf(fp, "\n");
  }

  fclose(fp);

  return true;
}

bool
aemimage::save_image(const char *filename)
{
  if (conductivity == nullptr) {
    return false;
  }

  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  for (int j = 0; j < rows; j ++) {
    for (int i = 0; i < columns; i ++) {

      fprintf(fp, "%15.9f ", conductivity[j * columns + i]);

    }

    fprintf(fp, "\n");
  }

  fclose(fp);

  return true;
}

void aemimage::update_layer_thickness()
{
  logspace(rows, depth, layer_thickness);

  double sum = 0.0;
  for (auto &d: layer_thickness) {
    sum += d;
  }

  if (fabs(depth - sum) > 1.0e-3) {
    //
    // Sanity check to make sure log space working.
    //
    for (auto &d: layer_thickness) {
      printf("  %10.6f\n", d);
    }
    printf("s %10.6f\n", sum);
    
    throw AEMEXCEPTION("Mismatch in layer thickness sum and depth\n");
  }
}
