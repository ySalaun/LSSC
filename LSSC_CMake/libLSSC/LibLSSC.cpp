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
 * @file LibLSSC.cpp
 * @brief Library for LSSC algorithm
 *
 * @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include "LibLSSC.h"

using namespace std;

void trainL1(vector<float> &io_dict, vector<float> &i_noisy, unsigned p_nPatch, Parameters &params){
	// matrices used for the dictionnary update
	vector<float> A(params.k*params.k, 0.f);
	vector<float> B(params.m*params.k, 0.f);

	// sparse coefficients of the dictionnary
	vector<float> alpha(params.k, 0.f);

	// generation of p_nPatch iid patches 
	vector<float> iidPatches(p_nPatch, 0.f);
	// TODO

	// initialize patch vector
	vector<float> patch(params.sPatch, 0.f);

	for(unsigned i = 0; i < p_nPatch; ++i){
		// compute patch vector
		// TODO
		
		// LARS
		// TODO

		add_xyT(A, alpha, alpha);
		add_xyT(B, patch, alpha);

		// Update
		// TODO
	}
}

    