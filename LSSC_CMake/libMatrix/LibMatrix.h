/*
* Copyright (c) 2013, Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
* All rights reserved.
*
* This program is free software: you can use, modify and/or
* redistribute it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version. You should have received a copy of this license along
* this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef LIBMATRIX_H_INCLUDED
#define LIBMATRIX_H_INCLUDED

#include <vector>

#include "ClassMatrix.h"

using namespace std;


/**
* @brief Compute z = x +/- y.
*
* @param i_x : vector of size m;
* @param i_y : matrix of size m;
* @param o_z : will contain i_x +/- i_y of size m;
* @param p_minus : if true, compute x - y, otherwise compute x + y.
*
* @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
**/
int add(
  vector<float> const& i_x,
  vector<float> const& i_y,
  vector<float> &o_z,
  const bool p_minus = false);


/**
* @brief Compute the dot product z = x . y
*
* @param i_x : vector of size m;
* @param i_y : matrix of size m;
* @param o_z : will contain the float x . y;
* @param p_max : if > -1, the dot product is only done for the p_iMax first values.
*
* @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
**/
int dotProduct(
  vector<float> const& i_x,
  vector<float> const& i_y,
  float &o_z,
  const int p_max = -1);

#endif // LIBMATRIX_H_INCLUDED
