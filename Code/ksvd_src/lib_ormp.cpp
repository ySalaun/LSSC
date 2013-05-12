/*
 * Copyright (c) 2011, Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
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
 * @file lib_ormp.cpp
 * @brief ORMP process, based on the SPAMS toolbox by Julien Mairal
 *
 * @author Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lib_ormp.h"

using namespace std;

#define MIN(a,b) (((a) > (b)) ? (b) : (a))

/**
 * @brief compute an Orthogonal Recursive Matching Pursuit
 *
 * @param X : (Np x n) X[i] is the vector to be represented in D
 * @param D : (k x n)	 dictionary (D[j] is the j-th atom)
 * @param A : (k x Np) sparse coordinates : A[j][i] is the
 *            coordinate of X[i] on D[j]
 **/

//! The matrix are accessed in column order i.e. X[i] is the i-th column of X
//! CONVENTION : we will denote by A_B the variable containing the
//! matrix transpose(A)*B (which contains the scalar products between the columns of A
//! and the columns of B.
void ormp_process(matD_t      &X,
                  matD_t      &D,
                  matU_t      &ind_v,
                  matD_t      &val_v,
                  unsigned     L,
                  const double eps)
{
    //! Declarations
    const unsigned n = X[0].size();
    const unsigned Np = X.size();
    const unsigned k = D.size();

    //! Initializations
    if (L <= 0)
        return;

    L = MIN(n, MIN(L,k));

    #pragma omp parallel shared(ind_v, val_v, X, D)
    {
        //! Declarations
        vecD_t norm(k), scores(k), x_T(k);
        matD_t A(L, vecD_t(L)), D_ELj(k, vecD_t(L)), D_D(k, vecD_t(k)), D_DLj(k, vecD_t(L));

        //! Compute the scalar products between the atoms of the dictionary
        for (unsigned j = 0; j < k; j++)
        {
            iterD_t it_D_D = D_D[j].begin();
            for (unsigned i = 0; i < k; i++, it_D_D++)
            {
                iterD_t D_i = D[i].begin();
                iterD_t D_j = D[j].begin();
                long double val = 0;
                for (unsigned s = 0; s < n; s++, D_i++, D_j++)
                    val += (long double) (*D_i) * (long double) (*D_j);
                (*it_D_D) = (double) val;
            }
        }

        #pragma omp for schedule(dynamic) nowait
            for (unsigned i = 0; i < Np; ++i)
            {
                //! Initialization
                ind_v[i].clear();
                val_v[i].clear();
                iterD_t X_i = X[i].begin();

                //! Compute the norm of X[i]
                long double normX = 0.0l;
                for (iterD_t it_x = X_i; it_x < X[i].end(); it_x++)
                    normX += (long double) (*it_x) * (long double) (*it_x);

                //! Compute the scalar products between X[i] and the elements
                //! of the dictionary
                for (unsigned j = 0; j < k; j++)
                {
                    iterD_t it_d = D[j].begin();
                    iterD_t it_x = X_i;
                    long double val = 0.0l;
                    for (unsigned s = 0; s < n; s++, it_d++, it_x++)
                        val += (long double) (*it_d) * (long double) (*it_x);
                    x_T[j] = (double)val;
                }

                coreORMP(D, D_D, scores, norm, A, D_ELj, D_DLj, x_T,
                                            ind_v[i], val_v[i], eps, (double) normX);

            } //! End of I-lop

    } //! End of parallel section
}

/**
 * @brief Sub function of ormp_process
 *
 * @param D : dictionary;
 * @param D_D : precomputed matrix D * D';
 * @param scores[i] = <x, e_i^j>^2
 * @param norm[i] = ||t_i^j||^2 = 1 - \sum_{p=0}^{j_1} <d_i, e_l_p>^2
 * @param A : sparse coordinates
 * @param D_ELj[i][j] : equals <d_i, e_l_j> at the end of the j-loop
 * @param D_DLj[i][s] = <d_i, d_l_s>
 * @param x_T[i] = <x, t_i^j> = <x, d_i> - \sum_{p=0}^{j-1} <d_i, e_l_p> <x, e_l_p>
 * @param ind : will contain index of atoms of D at the end of the j-loop
 * @param coord[q] = \alpha_l_q = \sum_{p=q}^j <x, e_l_p> a_{pq} : coordinate
 *                 of x on d_l_q at the end of the j-loop
 * @param eps : break condition
 * @param normr = ||x||^2 - \sum_{p=0}^j <x, e_l_p>^2
 *
 * @return none.
 **/
void coreORMP(matD_t      &D,
              matD_t      &D_D,
              vecD_t      &scores,
              vecD_t      &norm,
              matD_t      &A,
              matD_t      &D_ELj,
              matD_t      &D_DLj,
              vecD_t      &x_T,
              vecU_t      &ind,
              vecD_t      &coord,
              const double eps,
              double       normr)
{
    //! Declarations
    const unsigned L = A.size();
    const unsigned p = D.size();

    if (normr <= eps || L == 0)
        return;

    //! Initializations
    scores = x_T;
    for (iterD_t it = norm.begin(); it < norm.end(); it++)
        *it = 1.0l;

    for (matD_t::iterator it_A = A.begin(); it_A < A.end(); it_A++)
        for (iterD_t it = it_A->begin(); it < it_A->end(); it++)
            *it = 0.0l;

    //! Declarations
    iterD_t A_j, it_A_j, it_D_ELj_lj, it_gs, x_T_i, norm_i, scores_i, it_A_j_tmp;
    iterD_t D_lj, it_D_DLj;
    matD_t::iterator it_D_ELj, it_A;
    vecD_t x_el;

    //! Loop over j
    unsigned j;
    for (j = 0; j < L; j++)
    {
        //! Initialization
        const unsigned lj = ind_fmax(scores);
        A_j = A[j].begin();

        //! Stop if we cannot inverse norm[lj]
        if (norm[lj] < 1e-6)
            break;

        const long double invNorm = 1.0l / (long double) sqrtl(norm[lj]);
        const long double x_elj = (long double) x_T[lj] * invNorm;
        const long double delta = x_elj * x_elj;

        x_el.push_back(x_elj);
        coord.push_back(x_elj); //! The coordinate of x on the last chosen vector
        normr  = (double) ((long double) normr - delta);//! Update of the residual
        ind.push_back(lj); //! Memorize the chosen index

        //! Gram-Schmidt Algorithm, Update of A
        it_A_j = A_j;
        it_D_ELj_lj = D_ELj[lj].begin();
        for (unsigned i = 0; i < j; i++, it_A_j++, it_D_ELj_lj++)
            (*it_A_j) = (*it_D_ELj_lj);

        it_A_j = A_j;
        for (unsigned i = 0; i < j; i++, it_A_j++)
        {
            long double sum = 0.0l;
            it_A = A.begin() + i;
            it_A_j_tmp = it_A_j;

            for (unsigned s = 0; s < j - i; s++, it_A++, it_A_j_tmp++)
                sum -= (long double) (*it_A)[i] * (long double) (*it_A_j_tmp);

            (*it_A_j) = (double) (sum * invNorm);
        }
        (*it_A_j) = (double) invNorm;

        if (j == L-1 || (normr <= eps))
        {
            j++;
            break;
        }

        //! Update of D_DLj
        it_D_DLj = D_D[lj].begin();
        for (unsigned i = 0; i < p; i++, it_D_DLj++)
            D_DLj[i][j] = (*it_D_DLj);

        //! Compute the scalar product D[j]_D[lj] and memorize it now
        long double val = 0.0l;
        D_lj = D[lj].begin();
        for (iterD_t D_j = D[j].begin(); D_j < D[j].end(); D_j++, D_lj++)
            val += (long double) (*D_j) * (long double) (*D_lj);
        D_DLj[j][j] = (double) val;

        //! Update of D_ELj, x_T, norm, and scores.
        x_T_i = x_T.begin();
        norm_i = norm.begin();
        scores_i = scores.begin();
        it_D_ELj = D_ELj.begin();
        for (unsigned i = 0; i < p; i++, x_T_i++, norm_i++, scores_i++, it_D_ELj++)
        {
            long double val = 0;
            it_A_j = A_j;
            it_gs = D_DLj[i].begin();
            for (unsigned s = 0; s < j + 1; s++, it_A_j++, it_gs++)
                val += (long double) (*it_A_j) * (long double) (*it_gs);

            (*it_D_ELj)[j] = (double) val;
            (*x_T_i)      = (double) ((long double) (*x_T_i) - x_elj * val);
            (*norm_i)     = (double) ((long double) (*norm_i) -  val * val);
            (*scores_i)   = (double) ((long double) (*x_T_i) * (long double) (*x_T_i) \
                                 / (long double) (*norm_i));
        }

        for (iterU_t ind_s = ind.begin(); ind_s <= ind.begin() + j; ind_s++)
            scores[*ind_s] = double();
    }

    //! Compute the final coordinates of x on the chosen atoms of the dictionary
    iterD_t it_coord = coord.begin();
    for (unsigned i = 0; i < j; i++, it_coord++)
    {
        long double sum = 0;
        matD_t::iterator it_a = A.begin() + i;
        iterD_t it_x = x_el.begin() + i;
        for (unsigned s = 0; s < j - i; s++, it_a++, it_x++)
            sum += (long double) (*it_a)[i] * (long double) (*it_x);
        (*it_coord) = (double) sum;
    }
}

/**
 * @brief Find the index of the maximum absolute value
 *        of a table
 *
 * @param V : table of values.
 *
 * @return : return the index of the maximum absolute
 *           value of V.
 **/
unsigned ind_fmax(vecD_t &V)
{
    double   val = 0.0f;
    unsigned ind = 0;
    iterD_t it_v;
    unsigned j = 0;

    //! Find the index of the maximum of V
    for (it_v = V.begin(); it_v < V.end(); it_v++, j++)
        if (val < fabs(*it_v))
        {
            val = fabs(*it_v);
            ind = j;
        }

    return ind;
}



