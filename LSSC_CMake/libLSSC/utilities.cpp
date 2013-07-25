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

int add_xyT(Matrix &A, vector<float> &x, vector<float> &y){
	unsigned k = x.size();
	unsigned m = y.size();
	
	// size issue
	if(k*m != A.nCol*A.nRow){
		cout << "size issue in function add_xyT" << endl;
		return 0;
	}

	for(unsigned i = 0; i < k; ++i){
		for(unsigned j = 0; j < m; ++j){
			A.matrix[i*m + j] += x[i] * y[j];
		}
	}

	return 1;
}

int product_AB(const Matrix &A, const Matrix &B, Matrix &AB, const bool transpose){
	// TODO: the transpose part
  unsigned K = B.nRow;

	// size issue
	if(A.nCol != K){
		cout << "size issue in function product_AB" << endl;
		return 0;
	}

	unsigned m = A.nRow;
	unsigned n = B.nCol;
	float sum;

	for(unsigned i = 0; i < m; ++i){
		for(unsigned j = 0; j < n; ++j){
			sum = 0.f;
			for(unsigned k = 0; k < K; ++k){
				sum += A.matrix[i*K + k]*B.matrix[k*n + j];
			}
			AB.matrix[i*n + j] = sum;
		}
		
	}

	/*std::vector<float> q0(K);
	float * const pq0 = &q0[0];

	for (unsigned i = 0; i < n; i++) 
	{
		const float *pB = &B.matrix[i];
		for (unsigned k = 0; k < K; k++, pB += n)
		{
			pq0[k] = *pB;
		}

		float * pO = &AB.matrix[i];
		const float * pA = &A.matrix[0];
		for (unsigned j = 0; j < m; j++, pO += n, pA += K)
		{
			float z = 0.f;
			for (unsigned k = 0; k < K; k++)
			{
				z += pA[k] * pq0[k];
			}

			*pO = z;
		}
	}*/

	return 1;
}

vector<float> product_Ax(const Matrix &A, const vector<float> &x, const bool transpose){
	unsigned n = x.size();
	
	// size issue
	if((A.nCol != n && !transpose) || (A.nRow != n && transpose)){
		cout << "size issue in function product_Ax" << endl;
		return vector<float>(n, 0.f);
	}

	unsigned m = (transpose)? A.nCol : A.nRow;
	vector<float> Ax(m, 0.f);
	float sum;

	for(unsigned i = 0; i < m; ++i){
		sum = 0;
		for(unsigned j = 0; j < n; ++j){
			if(transpose){
				sum += A.matrix[j*m + i]*x[j];
			}
			else{
				sum += A.matrix[i*n + j]*x[j];
			}
		}
		Ax[i] = sum;
	}

	return Ax;
}

vector<float> add(const vector<float> &A, const vector<float> &B, bool minus){
	unsigned n = A.size();
	int sign = (minus)? -1:1;

	// size issue
	if(B.size() != n){
		cout << "size issue in function add" << endl;
		return vector<float>(n, 0.f);
	}

	vector<float> C(n, 0.f);

	for(unsigned i = 0; i < n; ++i){
		C[i] = A[i] + sign * B[i];
	}

	return C;
}

int add(const Matrix &A, const Matrix &B, Matrix &C, const bool minus){
	// size issue
	if(B.nRow != A.nRow || B.nCol != A.nCol){
		cout << "size issue in function add" << endl;
		return 0;
	}

  C.matrix = add(A.matrix, B.matrix, minus);

	return 1;
}

float dotProduct(const vector<float> x, const vector<float> y){
	unsigned n = x.size();
	
	// size issue
	if(n != y.size()){
		cout << "size issue in function dot product" << endl;
		return 0.f;
	}

	float res = 0.f;
	for(unsigned i = 0; i < n; ++i){
		res += x[i]*y[i];
	}
	return res;
}
    