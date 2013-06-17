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

#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

// small class for Matrices
// elements are accessed this way:
// M(i,j) = M.matrix[i * M.nCol + j];
// where i is the row and j the column
class Matrix{
public:
	int nRow, nCol;
	float *matrix;

	Matrix(int m, int n){
		nRow = m;
		nCol = n;
		matrix = new float[n*m];
		for(unsigned i = 0; i < n*m; ++i){
			matrix[i] = 0.f;
		}
	}

	void setMatrix(int m, int n, float *mat){
		delete[] matrix;
		nRow = m;
		nCol = n;
		matrix = new float[m*n];
		for(unsigned i = 0; i < m*n; ++i){
			matrix[i] = mat[i];
		}
	}

	~Matrix(){
		delete[] matrix;
	}
};

/**
 * @brief Compute A = A + xy' where y' is the transposed version of y
 *
 * @param A : matrix of size k x m;
 * @param x : vector of size k;
 * @param y : vector of size m.
 *
 * @return 0 if size issue and 1 else
 **/
int add_xyT(Matrix &A, vector<float> &x, vector<float> &y);

int product_AB(Matrix &A, Matrix &B, Matrix &AB);
vector<float> product_Ax(Matrix &A, vector<float> &x);

/**
 * @brief Compute C = A + B
 *
 * @param A : matrix/vector of size n;
 * @param B : matrix/vector of size n;
 *
 * @return 0 if size issue and A + B else 1
 **/
vector<float> add(vector<float> &A, vector<float> &B, bool minus = false);
int add(Matrix &A, Matrix &B, Matrix &C, bool minus = false);

#endif // UTILITIES_H_INCLUDED
