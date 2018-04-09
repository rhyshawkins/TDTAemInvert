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
#ifndef chainhistory_pixel_hpp
#define chainhistory_pixel_hpp

#include <vector>

class aemimage;

class PixelPerturbation {
public:

  PixelPerturbation();
  PixelPerturbation(bool accepted,
		    int idx,
		    double oldvalue,
		    double newvalue);

  bool accepted;
  int idx;
  double oldvalue;
  double newvalue;
};

class ChainHistoryPixel
{
public:
  ChainHistoryPixel(int rows, int columns);
  ChainHistoryPixel(aemimage &initial_model);
  
  static ChainHistoryPixel *load(const char *filename);

  bool save(const char *filename);
  
  int rows;
  int columns;
  double *initial_image;
  std::vector<PixelPerturbation> history;
};

#endif // chainhistory_pixel.hpp


  
