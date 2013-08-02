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

void display(const Matrix & mat, const int iMax){
  int m = (iMax == -1)? mat.nRow: iMax;
  int n = (iMax == -1)? mat.nCol: iMax;
  for(int i = 0; i< m; ++i){
    for(int j = 0; j< n; ++j){
      cout << mat.matrix[i*mat.nCol + j] << "/";
    }
    cout << endl;
  }
}

int add_xyT(Matrix &A, vector<float> &x, vector<float> &y, const int iMax){
	unsigned int k = x.size();
	unsigned int m = y.size();

	// size issue
	if( (iMax == -1 && (k != A.nRow || m != A.nCol))
    ||(iMax >=  0 && (k < (unsigned int) iMax || m < (unsigned int) iMax))){
		cout << "size issue in function add_xyT" << endl;
		return 0;
	}

  if(iMax > -1){
    k = iMax;
    m = iMax;
  }

	for(unsigned int i = 0; i < k; ++i){
		for(unsigned int j = 0; j < m; ++j){
      A.matrix[i*A.nCol + j] += x[i] * y[j];
		}
	}

	return 1;
}

int product_AB(const Matrix &A, const Matrix &B, Matrix &AB, const bool transpose){
  unsigned K = (transpose)? A.nRow:A.nCol;

	// size issue
	if(B.nRow != K){
		cout << "size issue in function product_AB" << endl;
		return 0;
	}

	unsigned m = (transpose)? A.nCol:A.nRow;
	unsigned n = B.nCol;
	float sum;

	for(unsigned i = 0; i < m; ++i){
		for(unsigned j = 0; j < n; ++j){
			sum = 0.f;
			for(unsigned k = 0; k < K; ++k){
        if(transpose){
          sum += A.matrix[k*A.nCol + i]*B.matrix[k*B.nCol + j];
        }
        else{
          sum += A.matrix[i*A.nCol + k]*B.matrix[k*B.nCol + j];
        }
			}
			AB.matrix[i*n + j] = sum;
		}

	}

  // TODO keep it or change it code from Marc
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

vector<float> product_Ax(const Matrix &A, const vector<float> &x, const bool transpose, const int iMax){
	unsigned int n = x.size();

	// size issue
	if( (iMax == -1 && ((A.nCol != n && !transpose) || (A.nRow != n && transpose)))
    ||(iMax >=  0 && ((A.nCol < (unsigned int) iMax && !transpose)
    || (A.nRow < (unsigned int) iMax && transpose) || (n < (unsigned int) iMax)))
    ){
		cout << "size issue in function product_Ax" << endl;
		return vector<float>(n, 0.f);
	}

	unsigned int m = (transpose)? A.nCol : A.nRow;
    m = (iMax > -1 && m > (unsigned int) iMax)? (unsigned int) iMax : m;
	vector<float> Ax(m, 0.f);
	float sum;

	for(unsigned i = 0; i < m; ++i){
		sum = 0;
		for(unsigned j = 0; j < n; ++j){
			if(transpose){
        sum += A.matrix[j*A.nCol + i]*x[j];
			}
			else{
        sum += A.matrix[i*A.nCol + j]*x[j];
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

float dotProduct(const vector<float> x, const vector<float> y, const int iMax){
	unsigned n = x.size();

	// size issue
	if( (iMax == -1 && n != y.size())
    ||(iMax >=  0 && (n < (unsigned int) iMax || y.size() < (unsigned int) iMax) )
    ){
		cout << "size issue in function dot product" << endl;
		return 0.f;
	}

    n = (iMax > -1) ? iMax : n;

	float res = 0.f;
	for(unsigned i = 0; i < n; ++i){
		res += x[i]*y[i];
	}
	return res;
}
