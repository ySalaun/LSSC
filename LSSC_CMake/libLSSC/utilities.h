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

vector<float> product_ABj(Matrix &A, Matrix &B, int col);
vector<float> product_Ax(Matrix &A, vector<float> &x);

/**
 * @brief Compute C = A + B
 *
 * @param A : matrix/vector of size n;
 * @param B : matrix/vector of size n;
 *
 * @return 0 if size issue and A + B else
 **/
vector<float> add(vector<float> &A, vector<float> &B);
Matrix add(Matrix &A, Matrix &B);

#endif // UTILITIES_H_INCLUDED
