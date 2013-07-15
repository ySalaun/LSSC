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

// train the algorithm with the L1 norm
void trainL1(Matrix &io_dict, vector<float> &i_noisy, unsigned p_nPatch, Parameters &params){
	// matrices used for the dictionnary update
	Matrix A(params.k, params.k);
	Matrix B(params.m, params.k);

	// sparse coefficients of the dictionnary
	vector<float> alpha(params.k, 0.f);

	// generation of p_nPatch iid patches 
	vector<int> iidPatches = randPatches(p_nPatch, params.nPatch);

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
		updateDictionary(io_dict, A, B, params);
	}

	// free memory
	A.~Matrix();
	B.~Matrix();
}

// generation of random patches
vector<int> randPatches(int nPatch, int nPatchMax){
	// list with all possibilities
	vector<int> list(nPatchMax, 0);	
	for(unsigned i = 0; i < nPatchMax; ++i){
		list[i] = i;
	}

	// list of random patches
	vector<int> randPatches(nPatch, 0);

	// seed the random function
	srand (time(NULL));

	// find the nPatch iid patches
	for(unsigned i = 0; i < nPatch; ++i){
		unsigned index = rand() % list.size();
		randPatches[i] = list[index];
		list.erase(list.begin() + index);
	}

	return randPatches;
}

// lars algorithm that minimizes ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
vector<float> lars(const Matrix &p_dict, const vector<float> &p_patch, Parameters &params){
	// code
	vector<float> alpha(params.k);

	// noisy picture norm
	// TODO: generate norm computation
	float norm = 0;

	// most correlated element
	float correlation = 0.f;
	int currentIndex = 0;
	// TODO: compute max of correlation


	// sort the coefficients inside alpha
	return alpha;
}

// dictionary update algorithm
void updateDictionary(Matrix &D, Matrix &A, Matrix &B, Parameters &params){
	// initialize dictionnary column vector
	Matrix u(params.m, params.k);
	
	float norm2;

	for(unsigned iter = 0; iter < params.update_iteration; ++iter){
		// u = DA
		product_AB(D, A, u);
		
		// u = B - u
		add(B, u, u, true);
		
		// uj = uj/ajj for each column j of u
		for(unsigned j = 0; j < params.k; ++j){
			for(unsigned i = 0; i < params.m; ++i){
				u.matrix[i*params.k + j] /= A.matrix[j*params.k + j];
			}
		}

		// u = u + D
		add(u, D, u);
		
		// uj = uj / max(||u||_2, 1) for each column j
		for(unsigned j = 0; j < params.k; ++j){
			norm2 = 0.f;
			for(unsigned i = 0; i < params.m; ++i){
				norm2 += u.matrix[i*params.k + j] * u.matrix[i*params.k + j];
			}
			if(norm2 > 1){
				norm2 = sqrt(norm2);
				for(unsigned i = 0; i < params.m; ++i){
					u.matrix[i*params.k + j] /= norm2;
				}
			}
		}

		// D = copy(u)
		D.setMatrix(params.m, params.k, u.matrix);
	}
	// free memory
	u.~Matrix();
}