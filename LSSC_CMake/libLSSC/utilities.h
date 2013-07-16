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

	// set a Matrix as a Gram Matrix version of the input
	void setGram(Matrix D){
		nRow = D.nCol;
		nCol = D.nCol;
		delete[] matrix;
		matrix = new float[nRow*nRow];
		for(unsigned i = 0; i < nCol; ++i){
			for(unsigned j = 0; j < nCol; ++j){
				float sum = 0;
				for(unsigned k = 0; k < D.nCol; ++k){
					sum += D.matrix[k * nCol + i] * D.matrix[k * nCol + j];
				}
				matrix[i * nCol + j] = sum;
			}
		}
	}

	void copyRow(Matrix M, int i_from, int i_to){
		if(M.nCol != nCol){
			cout << "error in row size" << endl;
		}
		else{
			for(unsigned j=0; j<nCol; ++j){
				matrix[i_to * nCol + j] = M.matrix[i_from * nCol + j];
			}
		}
	}

	void copyCol(Matrix M, int j_from, int j_to){
		if(M.nRow != nRow){
			cout << "error in column size" << endl;
		}
		else{
			for(unsigned i=0; i<nRow; ++i){
				matrix[i * nCol + j_to] = M.matrix[i * nCol + j_from];
			}
		}
	}

	void symmetrizeUpperPart(){
		if(nRow != nCol){
			cout << "it is not a square matrix" << endl;
		}
		for(unsigned i=0; i<nRow; ++i){
			for(unsigned j=i+1; j<nCol; ++j){
				matrix[j * nCol + i] = matrix[i * nCol + j];
			}
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
vector<float> product_Ax(const Matrix &A, const vector<float> &x, const bool transpose = false);

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

float dotProduct(const vector<float> x, const vector<float> y);

#endif // UTILITIES_H_INCLUDED
