/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
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
 * @file LibLars.cpp
 * @brief LARS functions, based on Julien Mairal's SPAM toolbox.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "LibLars.h"

using namespace std;

/**
 * @brief Auxiliary function for lasso
 *      solve min_{ alpha } | | alpha | | _1 s . t . | | x-Dalpha | | _2^2 <= lambda
 *
 * @param io_DtR : D'x
 * @param i_G : Gram Matrix = D'D
 * @param Gs, Ga, invGs, u, coeffs, ind, work, normX seems to be references that are filled in the programm below
 * @param pos : positivity constraint
 *
 * @return none.
 **/
void coreLARS2new(
    std::vector<float>& io_DtR
,   std::vector<float> const& i_G
,   float& io_normX
,   const float p_constraint
,   const bool p_pos
){
    //! Size of matrices
    const unsigned K  = sqrt(i_G.size()); 	// size of G
    const unsigned LL = K;					// size of Gs ==> number max of iteration, NOT necessarily K !!!
    const unsigned L  = getMin(LL, K);

    //! Initializations
	// maybe better to fill them with zeros
    vector<float> Gs    (LL * LL);
    vector<float> Ga    (K  * K );
    vector<float> invGs (LL * LL);
    vector<float> work  (K  * K );
    vector<float> u     (LL);
    vector<float> coeffs(L);
    vector<int  > ind   (L);
    for (unsigned k = 0; k < L; k++) {
        coeffs[k] = 0.f;
        ind   [k] = -1;
    }

	//! Find the most correlated element ==> c_j in LARS article
	int currentInd = findMax(io_DtR, true);

	//! why issue if normX < constraint ?
	if (io_normX < p_constraint) {
	    return;
	}

	//! begin by adding a new atom
	bool newAtom = true;

	//! loop
	//! add one atom at each iteration except when gone too far and need to come back
	//! stop when the criterion is no more checked
	for (int i = 0; i < (int) L; i++) { // check that an int is really necessary for i

		//! new atom, update Ga, Gs, invGs
		if (newAtom) {
			//! add new index for the current most correlated element
			ind[i] = currentInd;
			
			//! Ga[i-th line] = G[ind[i]-th column]
			//! seems to be used for computation purpose with u
			for (unsigned k = 0; k < K; k++) {
                Ga[i * K + k] = i_G[k * K + ind[i]]; // To improve
			}
			//! Gs = triangular_inf(D_A^T D_A) where D_A = selected index of D (LARS article)
			for (int j = 0; j <= i; j++) {
				Gs[i * LL + j] = Ga[i * K + ind[j]]; // To improve
			}

			//! Update inverse of Gs : invGs
			if (i == 0) {
				//! Gs^(-1)(0,0) = 1/Gs(0,0)
				invGs[0] = 1.f / Gs[0];
			}
			else {
				//! u = Gs^(-1)*Gs[i-th line]
				//! WARNING: only the lower part of Gs^(-1) is correctly filled ....
                for (int p = 0; p < i; p++) {
                    float value = 0.f;
                    for (int q = 0; q < i; q++) {
						if(q > p){	//! use the symmetry property of the matrix
							value += invGs[q * LL + p] * Gs[i * LL + q];
						}
						else{
							value += invGs[p * LL + q] * Gs[i * LL + q];
						}
                    }
                    u[p] = value;
                }

				//! schur = 1/(Gs(i,i) -  u.Gs[i-th line])
				float dotProduct = 0.f;
				for (int k = 0; k < i; k++) {
                    dotProduct += u[k] * Gs[i * LL + k];
				}
				const float schur = 1.f / (Gs[i * LL + i] - dotProduct);

				//! Gs^(-1)(i,i) = 1/(Gs(i,i) -  u.Gs[i-th line])
				invGs[i * LL + i] = schur;

				//! Gs^(-1)[i-th line] = -schur * u
                for (int k = 0; k < i; k++) {
                    invGs[i * LL + k] = -schur * u[k];
                }

				//! Gs^(-1) = schur * uu' + Gs^(-1)
				//! WARNING: only the triangular inf part of Gs^(-1) is correct (the other is never used....)
                for (int p = 0; p <= i; p++) {
                    for (int q = 0; q <= p; q++) {
                        invGs[p * LL + q] += schur * u[p] * u[q];
                    }
                }
			}
		}

		//! Compute the path direction
		//! work[j] = sgn(c_j)
		for (int j = 0; j <= i; j++) {
			work[j] = (io_DtR[ind[j]] > 0 ? 1.f : -1.f);
		}

		//! u = Gs^(-1)*work
		//! WARNING: only the lower part of Gs^(-1) is correctly filled ....
        for (int p = 0; p < i + 1; p++) {
            float value = 0.f;
            for (int q = 0; q < i + 1; q++) {
                if(q > p){	//! use the symmetry property of the matrix
					value += invGs[q * LL + p] * work[q];
				}
				else{
					value += invGs[p * LL + q] * work[q];
				}
            }
            u[p] = value;
        }

		//! Compute the step on the path
		float step_max = INFINITY;
		int first_zero = -1;
		for (int j = 0; j <= i; j++) {
			const float ratio = -coeffs[j] / u[j];
			if (ratio > 0 && ratio <= step_max) {
				step_max   = ratio;
				first_zero = j;
			}
		}

		//! correl = abs(io_DtR(ind[0]))
		float current_correlation = fabs(io_DtR[ind[0]]);
		

		//! work[3rd row] = work[2nd row] = work[1st row] = Ga*u
        for (unsigned p = 0; p < K; p++) {
            float value = 0.f;
            for (int q = 0; q < i + 1; q++) {
                value += Ga[q * K + p] * u[q];
            }
            work[p + K * 2] = value;
            work[p + K    ] = value;
            work[p        ] = value;
        }

        //! work[1st/2nd row active indexes] = inf
		for (int j = 0; j <= i; j++) {
			work[ind[j]    ] = INFINITY;
			work[ind[j] + K] = INFINITY;
		}

		//! fill inactive indexes coeff of 1st col
		//! when positivity constraint, eliminate 1st col of work
		for (unsigned j = 0; j < K; j++) {
			work[j] = (((work[j] < INFINITY) && (work[j] > -1.f) && !p_pos)
                    ? (current_correlation + io_DtR[j]) / (1.f + work[j])
                    : INFINITY);
		}

		//! fill inactive indexes coeff of 2nd col
		for (unsigned j = 0; j < K; j++) {
			work[j + K] = (((work[j + K] < INFINITY) && (work[j + K] < 1.f))
                        ? (current_correlation - io_DtR[j]) / (1.f - work[j + K])
                        : INFINITY);
		}
		
		//! index of the minimum absolute value of work (1st and 2nd col only)
		//! new most correlated element ?
		int index  = findMin(work, 2 * K);
		float step = work[index];

		//! Choose next element
		currentInd = index % K;

		//! compute the coefficients of the polynome representing normX^2
		//! coeff1 = sum_{j <= i} sgn(io_DtR[ind[j]])*u[j]
		float coeff1 = 0;
		for (int j = 0; j <= i; j++) {
			coeff1 += (io_DtR[ind[j]] > 0 ? u[j] : -u[j]);
		}

		//! coeff2 = sum_{j <= i} io_DtR[ind[j]]*u[j]
		float coeff2 = 0;
		for (int j = 0; j <= i; j++){
			coeff2 += io_DtR[ind[j]] * u[j];
		}
		float coeff3 = io_normX - p_constraint;

		//! what is the meaning of delta ??
		//! delta = (sum_{j <= i} io_DtR[ind[j]]*u[j])^2 - (sum_{j <= i} sgn(io_DtR[ind[j]])*u[j])*(normX-constraint)
		const float delta = coeff2 * coeff2 - coeff1 * coeff3;

		//! what are the meanings of the step ?
		float step_max2;
		step_max2 = (delta < 0 ? INFINITY : (coeff2 - sqrt(delta)) / coeff1);
		step_max2 = getMin(current_correlation, step_max2);

		step = getMin(getMin(step, step_max2), step_max);

		//! stop the path
		if (step == INFINITY) {
		    break;
        }

		//! Update coefficients
		//! coeffs = step * u + coeffs
        for (int k = 0; k < i + 1; k++) {
            coeffs[k] += step * u[k];
        }

		//! if positivity constraint ==> cancel the negative coeffs
		if (p_pos) {
			for (int j = 0; j < i + 1; j++) {
				if (coeffs[j] < 0) {
                    coeffs[j] = 0;
				}
			}
		}

		//! Update correlations
        for (unsigned k = 0; k < K; k++) {
            io_DtR[k] -= step * work[2 * K + k];
        }

		//! Update normX
		//! normX = normX + (sum_{j <= i} sgn(io_DtR[ind[j]])*u[j]) * step^2
		//!         - 2 * (sum_{j <= i} io_DtR[ind[j]]*u[j]) * step
		io_normX += coeff1 * step * step - 2 * coeff2 * step;

		//! Choose next action
		//! case 1: when need to go back on the path cause of sign issues <== check
		//! case 2: next index, new atom
		if (step == step_max) {

            //! Downdate, remove first_zero
            //! Downdate Ga, Gs, invGs, ind, coeffs
			//! translate rows of Ga backward
			for (int j = first_zero; j < i; j++) {
                for (unsigned k = 0; k < K; k++) {
                    Ga[j * K + k] = Ga[(j + 1) * K + k];
                }
				ind   [j] = ind   [j + 1];
				coeffs[j] = coeffs[j + 1];
			}
			ind   [i] = -1;
			coeffs[i] = 0.f;

			//! translate rows of Gs
			for (int j = first_zero; j < i; j++) {
				for (int k = 0; k < first_zero; k++) {
                    Gs[j * LL + k] = Gs[(j + 1) * LL + k];
				}
				for (int k = first_zero; k < i; k++) {
                    Gs[j * LL + k] = Gs[(j + 1) * LL + k + 1];
				}
			}

			//! schur = Gs^(-1)[first_zero, first_zero]
			const float invSchur = 1.f / invGs[first_zero * LL + first_zero];

            //! Update u from Gs^(-1)
			for (int k = 0; k < first_zero; k++) {
                u[k] = invGs[first_zero * LL + k];
			}
			for (int k = first_zero; k < i; k++) {
                u[k] = invGs[(first_zero + 1) * LL + k];
			}

            //! Update Gs^(-1)
			for (int j = first_zero; j < i; j++) {
			    for (int k = 0; k < first_zero; k++) {
                    invGs[j * LL + k] = invGs[(j + 1) * LL + k];
			    }
			    for (int k = first_zero; k < i; k++) {
                    invGs[j * LL + k] = invGs[(j + 1) * LL + k + 1];
			    }
			}
			for (int p = 0; p < i; p++) {
                for (int q = 0; q <= p; q++) {
                    invGs[p * LL + q] -= invSchur * u[p] * u[q];
                }
            }
			newAtom = false;

			//! why i-2 ? i-1 ??
			i = i - 2;
		}
		else {
			newAtom = true;
		}

		if (abs(step) < 1e-15 ||
            step == step_max2 ||
            (io_normX < 1e-15) ||
            //(mode == L2ERROR && normX - constraint < 1e-15) ||
            (i == (L-1))) {
			break;
		}
	}
}

/**
 * @brief Find the index of the (absolute) maximum value of a vector.
 *
 * @param i_vec : input vector;
 * @param p_abs : if true, look for the absolute maximum, otherwise the maximum.
 *
 * @return the index of the maximum value of i_vec.
 **/
int findMax(
    std::vector<float> const& i_vec
,   const bool p_abs
){
    float max = i_vec[0];
    int ind = 0;
    for (unsigned k = 0; k < i_vec.size(); k++) {
        const float value = (p_abs ? fabs(i_vec[k]) : i_vec[k]);
        if (max < value) {
            max = value;
            ind = k;
        }
    }

    return ind;
}

/**
 * @brief Find the index of the minimum absolute value between the N first values of a vector.
 *
 * @param i_vec : input vector;
 * @param i_N : number of the first values on which we are looking for.
 *
 * @return the index of the minimum absolute value of the i_N first values of i_vec.
 **/
int findMin(
    std::vector<float> const& i_vec
,   const unsigned i_N
){
    float min = fabs(i_vec[0]);
    int ind = 0;
    for (unsigned k = 0; k < i_N; k++) {
        if (min > fabs(i_vec[k])) {
            min = fabs(i_vec[k]);
            ind = k;
        }
    }

    return ind;
}

/**
 * @brief Find the minimum of two values
 *
 * @param i_a : first value;
 * @param i_b : second value.
 *
 * @return min(i_a, i_b).
 **/
int getMin(
    const int i_a
,   const int i_b
){
    return (i_a < i_b ? i_a : i_b);
}

/**
 * @brief Find the minimum of two values
 *
 * @param i_a : first value;
 * @param i_b : second value.
 *
 * @return min(i_a, i_b).
 **/
unsigned getMin(
    const unsigned i_a
,   const unsigned i_b
){
    return (i_a < i_b ? i_a : i_b);
}

/**
 * @brief Find the minimum of two values
 *
 * @param i_a : first value;
 * @param i_b : second value.
 *
 * @return min(i_a, i_b).
 **/
float getMin(
    const float i_a
,   const float i_b
){
    return (i_a < i_b ? i_a : i_b);
}

/**
 * @brief Print a matrix.
 *
 * @param i_mat : matrix to print on screen;
 * @param p_row : number of rows;
 * @param p_col : number of cols;
 * @param p_name : name of the matrix.
 *
 * @return none.
 **/
void printMat(
    std::vector<float> const& i_mat
,   const unsigned p_row
,   const unsigned p_col
,   const char* p_name
){
    cout << p_name << " = " << endl;
    for (unsigned i = 0; i < p_row; i++) {
        for (unsigned j = 0; j < p_col; j++) {
            cout << i_mat[i * p_col + j] << " ";
        }
        cout << endl;
    }
}








