/*
 * Copyright (c) 2013, Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
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
 * @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include "utilities.h"

using namespace std;

int add_xyT(vector<float> &A, vector<float> &x, vector<float> &y){
	unsigned k = x.size();
	unsigned m = y.size();
	
	// size issue
	if(k*m != A.size()){
		return 0;
	}

	for(unsigned i = 0; i < k; ++i){
		for(unsigned j = 0; j < m; ++j){
			A[i*m + j] += x[i] * y[j];
		}
	}

	return 1;
}
    