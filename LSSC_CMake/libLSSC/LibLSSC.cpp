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

#define INFINITY	std::numeric_limits<float>::max()
#define EPSILON		1e-15

// train the algorithm with the L1 norm
void trainL1(Matrix &io_dict, const vector<float> &i_noisy, unsigned p_nPatch, const Parameters &params){

	// matrices used for the dictionnary update
	Matrix A(params.k, params.k);
	Matrix B(params.m, params.k);

	// sparse coefficients of the dictionnary
	vector<float> alpha(params.k, 0.f);

	// generation of p_nPatch iid patches 
	vector<int> iidPatches = randPatches(p_nPatch, params.nPatch);

	// initialize patch vector
	vector<float> patch(params.m, 0.f);

	for(unsigned i = 0; i < p_nPatch; ++i){
    // compute patch vector
		int I = iidPatches[i] / params.nRowPatches * params.sPatch;
		int J = iidPatches[i] % params.nRowPatches * params.sPatch;
		for(unsigned ii = 0; ii < params.sPatch; ++ii){
			for(unsigned jj = 0; jj < params.sPatch; ++jj){
				patch[ii * params.sPatch + jj] = i_noisy[(I + ii) * params.w + (J + jj)];
			}
		}
		
		// LARS
		lars(io_dict, patch, params);

		add_xyT(A, alpha, alpha);
		add_xyT(B, patch, alpha);

		// Update
		updateDictionary(io_dict, A, B, params);

		cout << i << " done" << endl;
	}
}

// generation of random patches
vector<int> randPatches(const int nPatch, const int nPatchMax){
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
vector<float> lars(const Matrix &p_dict, const vector<float> &p_patch, const Parameters &params){
	// code
	vector<float> alpha(params.k, 0);
	
	// active set indexes
	vector<int> activeIndexes(0);

	// noisy picture norm
	float norm = dotProduct(p_patch, p_patch);
	// if the norm is lower than the regularization parameter, stop the algorithm
	if(norm > params.reg){
		return alpha;
	}

	// most correlated element
	vector<float> correlation = product_Ax(p_dict, p_patch, true);
	float cMax = correlation[0];
	int currentIndex = 0;
	for(unsigned i=1; i<correlation.size(); ++i){
		if(correlation[i] > cMax){
			cMax = correlation[i];
			currentIndex = i;
		}
	}

	// begin by adding a new atom in the code
	bool newAtom = true;

	// matrices parameters
	Matrix G(params.k, params.k);
  product_AB(p_dict, p_dict, G, true);
	Matrix Ga(params.k, params.k);
	Matrix Gs(params.k, params.k);
	Matrix invGs(params.k, params.k);

	// add the new most correlated element at each iteration
	for(unsigned i = 0; i < params.k; ++i){
		/*------------UPDATE------------*/
		if(newAtom){
			activeIndexes.push_back(currentIndex);
			Ga.copyRow(G, currentIndex, i);
			Gs.copyCol(Ga, currentIndex, i);
			Gs.symmetrizeUpperPart();
      // TODO: update Gs^-1
		}

		/*-------VARIABLES UPDATE-------*/
		// compute sign vector
		// sgn = sgn(c)
		vector<float> sgn(i+1, 1.f);
		for(unsigned j = 0; j <= i; ++j){
			if(correlation[j] < 0){
				sgn[j] = -1.f;
			}
		}
		// compute direction vector
		// Ua = invGs * sgn
		// TODO change sth cause the size of Gs and invGs do not correspond with sgn one
		vector<float> Ua = product_Ax(invGs, sgn);
				
		/*-------------STEP-------------*/
		float step = INFINITY;
		float value;
		// current maximum of correlation
		cMax = correlation[0];
		Ga.nRow = i;
		vector<float> tGaUa = product_Ax(Ga, Ua, true);
		Ga.nRow = params.k;
		for(unsigned j = 0; j <= i; ++j){
			tGaUa[activeIndexes[j]] = INFINITY;
		}
		for(unsigned j = 0; j < params.k; ++j){
			if(tGaUa[j] == INFINITY){
				continue;
			}

			if(tGaUa[j] > -1){
				value = (cMax + correlation[j]) / (1 + tGaUa[j]);
				if(value < step){
					step = value;
					currentIndex = j;
				}
			}

			if(tGaUa[j] < 1){
				value = (cMax - correlation[j]) / (1 - tGaUa[j]);
				if(value < step){
					step = value;
					currentIndex = j;
				}
			}
		}

		/*---------DOWNDATE STEP--------*/
		// if this step is reached, downdate
		float ratio, stepDowndate = INFINITY;
		int indexDowndate;
		for(unsigned j = 0; j <= i; ++j){
			float ratio = -alpha[i]/Ua[i];
			if( ratio < stepDowndate && ratio >= 0 ){
				stepDowndate = ratio;
				indexDowndate = j;
			}
		}
		if(stepDowndate < 0){
			cout << "step downdate sign issue" << endl;
		}

		/*---------STOPPING STEP--------*/
		float a = dotProduct(sgn,Ua);
		float b = dotProduct(correlation, Ua);
		float c = norm - params.reg;
		float delta = b*b - a*c;
		float stepStop = (delta > 0)? (b - sqrt(delta))/a: INFINITY;
		stepStop = min(stepStop, cMax);

		/*---------TAKE THE STEP--------*/
		float finalStep = min(step, min(stepDowndate, stepStop));
		for(unsigned j = 0; j <= i; ++j){
			alpha[activeIndexes[j]] += finalStep*Ua[j];
		}
		for(unsigned j = 0; j < params.k; ++j){
			correlation[j] -= finalStep * tGaUa[j];
		}
		norm += (a * finalStep - 2 * b) * finalStep;

		/*------STOPPING CONDITIONS-----*/
		if( abs(finalStep) < EPSILON	||
			finalStep == stepStop		||
			norm < EPSILON				||
			norm - params.reg < EPSILON	){
			break;
		}

		/*-----------DOWNDATE-----------*/
		if( finalStep == stepDowndate){
			// TODO
			i -= 2;
			newAtom = false;
		}
		else{
			newAtom = true;
		}
	}	

	// return final code
	return alpha;
}

void updateGram(Matrix &invGs, Matrix &Gs, const unsigned iter, const Parameters &params){
  // case when the matrix is 1x1
  if(iter == 0){
    invGs.matrix[0] = 1/Gs.matrix[0];
  }
  else{
    // iRow = Gs[i-th row]
    vector<float> iRow = Gs.row(iter);
    
    // u = Gs^-1 Gs[iter-th row]
    vector<float> u = product_Ax(invGs, iRow);

    // sigma = 1/(Gs_ii - u. Gs[i-th row])
    float sigma = 1/((*Gs(iter, iter)) - dotProduct(u, iRow));

    //
    //invGs(iter, iter) = sigma; 
   
  }
}

// dictionary update algorithm
void updateDictionary(Matrix &D, const Matrix &A, const Matrix &B, const Parameters &params){
	time_t t1, t2, t3, t4;
	
	float norm2;
  Matrix u(params.m, params.k);
  //D = Matrix(D.nRow, D.nCol);
	for(unsigned iter = 0; iter < params.update_iteration; ++iter){
		// u = DA
		product_AB(D, A, u);
		
		// u = B - u
		add(B, u, u, true);
		
		// uj = uj/ajj for each column j of u
		for(unsigned j = 0; j < params.k; ++j){
			for(unsigned i = 0; i < params.m; ++i){
        // TODO Beware, this condition should really be there ?
        if(A.matrix[j*params.k + j] != 0){
          u.matrix[i*params.k + j] /= A.matrix[j*params.k + j];
        }
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

    // D = u
    D = u;
	}
}