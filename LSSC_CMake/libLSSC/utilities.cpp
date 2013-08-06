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
 * @file utilities.cpp
 * @brief Side functions (matrices operations...) for LSSC algorithm
 *
 * @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "utilities.h"

using namespace std;

void display(const char* msg,  const Parameters &params, bool endline){
  if(params.verbose){
    cout << msg;
    if(endline){
      cout << endl;
    }
  }
}

void display(const vector<float>& vec){
  for(unsigned int i = 0; i< vec.size(); ++i){
    cout << vec[i] << "/";
  }
  cout << endl;
}

void display(const Matrix2 & mat, const int iMax){
  unsigned int m, n;
  mat.getSize(m, n);
  
  m = (iMax == -1)? m: iMax;
  n = (iMax == -1)? n: iMax;

  for(int i = 0; i< m; ++i){
    for(int j = 0; j< n; ++j){
      cout << mat(i, j) << "/";
    }
    cout << endl;
  }
}