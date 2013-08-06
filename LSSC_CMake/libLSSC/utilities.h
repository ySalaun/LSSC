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

#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <stdlib.h>
#include <iostream>
#include <vector>

#ifdef __linux__
    #include "../Main/params.h"
    #include "../libMatrix/LibMatrix.h"
#else
    #include "Main/params.h"
    #include "libMatrix/LibMatrix.h"
#endif

using namespace std;

// small class for Matrices
// elements are accessed this way:
// M(i,j) = M.matrix[i * M.nCol + j];
// where i is the row and j the column
class Matrix{
public:
	unsigned int nRow, nCol; // no need to be integer
	vector<float> matrix;

	Matrix(int m, int n):
    nRow(m),
		nCol(n),
    matrix(n*m, 0.f)
  {}

  Matrix(Matrix &m):
    nRow(m.nRow),
    nCol(m.nCol),
    matrix(m.matrix)
  {}

	Matrix& operator=(Matrix& mat){
    nRow = mat.nRow;
    nCol = mat.nCol;
    matrix = mat.matrix;
    return *this;
	}

	float& operator()(const unsigned i, const unsigned j){
    return matrix[i * nCol + j];
  }

  // set a Matrix as a Gram Matrix version of the input
	void setGram(const Matrix &D){
		nRow = D.nCol;
		nCol = D.nCol;
		matrix.assign(nRow*nRow, 0.f);
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

	void copyRow(const Matrix &M, int i_from, int i_to){
		if(M.nCol != nCol){
			cout << "error in row size" << endl;
		}
		else{
			for(unsigned j=0; j<nCol; ++j){
				matrix[i_to * nCol + j] = M.matrix[i_from * nCol + j];
			}
		}
	}

  vector<float> row(const unsigned iRow, const unsigned length){
    vector<float> row(length);
    for(unsigned j=0; j<length; ++j){
      row[j] = matrix[iRow * nCol + j];
    }
    return row;
  }

	void copyCol(const Matrix &M, int j_from, int j_to){
		if(M.nRow != nRow){
			cout << "error in column size" << endl;
		}
		else{
			for(unsigned i=0; i<nRow; ++i){
				matrix[i * nCol + j_to] = M.matrix[i * nCol + j_from];
			}
		}
	}

	void symmetrizeUpperPart(unsigned int iMax){ // no need to be an integer
		if(nRow != nCol){
			cout << "it is not a square matrix" << endl;
		}
		for(unsigned i=0; i<iMax; ++i){
			for(unsigned j=i+1; j<iMax; ++j){
				matrix[j * nCol + i] = matrix[i * nCol + j];
			}
		}
	}

	~Matrix(){
	}
};

void display(const char* msg, const Parameters &params, bool endline = true);
void display(const vector<float>& vec);
void display(const Matrix & mat, const int iMax = -1);
/**
 * @brief Compute A = A + xy' where y' is the transposed version of y
 *
 * @param A : matrix of size k x m;
 * @param x : vector of size k;
 * @param y : vector of size m.
 *
 * @return 0 if size issue and 1 else
 **/
int add_xyT(Matrix &A, vector<float> &x, vector<float> &y, const int iMax = -1);

int product_AB(const Matrix &A, const Matrix &B, Matrix &AB, const bool transpose = false);
vector<float> product_Ax(const Matrix &A, const vector<float> &x, const bool transpose = false, const int iMax = -1);

/**
 * @brief Compute C = A + B
 *
 * @param A : matrix/vector of size n;
 * @param B : matrix/vector of size n;
 *
 * @return 0 if size issue and A + B else 1
 **/
vector<float> add(const vector<float> &A, const vector<float> &B, bool minus = false);
int add(const Matrix &A, const Matrix &B, Matrix &C, const bool minus = false);

float dotProduct(const vector<float> x, const vector<float> y, const int iMax = -1);

#endif // UTILITIES_H_INCLUDED
