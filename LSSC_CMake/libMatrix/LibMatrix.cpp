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

/**
* @file LibMatrix.cpp
* @brief Small matrix library.
*
* @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
**/

#include "LibMatrix.h"

#include <stdlib.h>
#include <iostream>

using namespace std;


//! Compute z = x +/- y.
int add(
  vector<float> const& i_x,
  vector<float> const& i_y,
  vector<float> &o_z,
  const bool p_minus) {

    //! Initialization
    const unsigned int m = i_x.size();
    const float sign     = (p_minus ? -1.f : 1.f);
    if (o_z.size() != m) {
      o_z.resize(m);
    }

    //! Check size
    if (i_y.size() != m) {
      cerr << "addxy - error : vector sizes aren't consistent" << endl;
      return EXIT_FAILURE;
    }

    //! Compute the addition
    for (unsigned int i = 0; i < m; i++) {
      o_z[i] = i_x[i] + sign * i_y[i];
    }

    return EXIT_SUCCESS;
}


//! Compute the dot product z = x . y
int dotProduct(
  vector<float> const& i_x,
  vector<float> const& i_y,
  float &o_z,
  const int p_max) {

  //! Initialization
  const unsigned int m = (p_max > -1 ? p_max : i_x.size());
  o_z = 0;

  //! Check size
  if ((p_max == -1 && i_x.size() != i_y.size()) || p_max > (int) i_x.size() || p_max > (int) i_y.size()) {
    cerr << "dotProduct - error : vector sizes aren't consistent" << endl;
    return EXIT_FAILURE;
  }

  //! Compute the dot product
  for (unsigned int i = 0; i < m; i++) {
    o_z += i_x[i] * i_y[i];
  }

  return EXIT_SUCCESS;
}












