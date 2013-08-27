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
 * @file LibOrmp.cpp
 * @brief ORMP process, based on the SPAMS toolbox by Julien Mairal
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "LibOrmp.h"



using namespace std;


//! Compute the Orthogonal Recursive Matching Pursuit over i_X.
void computeORMP(
  std::vector<double> const& i_X,
  std::vector<double> const& i_D,
  std::vector<std::vector<unsigned int> > &o_indV,
  std::vector<std::vector<double> > &o_valV,
  const double p_eps,
  const unsigned int p_L,
  const unsigned int p_M,
  const unsigned int p_N){

  //! Declarations
  const unsigned int K = std::min(p_N, p_M);
  o_indV.resize(p_L);
  o_valV.resize(p_L);

#ifdef _OPENMP
#pragma omp parallel shared(o_indV, o_valV, i_X, i_D)
#endif
  {
    //! Declarations
    vector<double> Xt(p_M);
    vector<double> DDt(p_M * p_M);

    //! Compute the scalar products between the atoms of the dictionary
    for (unsigned int j = 0; j < p_M; j++) {
      double* iDDt = &DDt[j * p_M];
      const double* iDj = &i_D[j * p_N];

      for (unsigned int i = 0; i < p_M; i++) {
        const double* iDi = &i_D[i * p_N];
        long double val = 0;

        for (unsigned int k = 0; k < p_N; k++) {
          val += (long double) iDi[k] * (long double) iDj[k];
        }

        iDDt[i] = (double) val;
      }
    }

#ifdef _OPENMP
  #pragma omp for schedule(dynamic) nowait
#endif
    for (unsigned int i = 0; i < p_L; i++) {

      //! Initialization
      const double* iX = &i_X[i * p_N];

      //! Compute the norm of the i-th row of i_X
      long double normX = 0.l;
      for (unsigned int j = 0; j < p_N; j++) {
        normX += (long double) iX[j] * (long double) iX[j];
      }

      //! Compute the scalar products between the i-th row of i_X
      //! and the elements of the dictionary
      for (unsigned int j = 0; j < p_M; j++) {
        const double* iD = &i_D[j * p_N];
        long double value = 0.l;

        for (unsigned k = 0; k < p_N; k++) {
          value += (long double) iD[k] * (long double) iX[k];
        }
        Xt[j] = (double) value;
      }

      coreORMP(i_D, DDt, Xt, o_indV[i], o_valV[i], p_eps, (double) normX, p_M, p_N, K);

    } //! End of I-lop

  } //! End of parallel section
}


//! Sub routine of computeORMP.
void coreORMP(
  std::vector<double> const& i_D,
  std::vector<double> const& i_DDt,
  std::vector<double> &io_Xt,
  std::vector<unsigned int> &o_ind,
  std::vector<double> &o_val,
  const double p_eps,
  const double p_normr,
  const unsigned int p_M,
  const unsigned int p_N,
  const unsigned int p_K) {

  //! Break condition
  double normr = p_normr;
  if (p_normr <= p_eps || p_K == 0) {
    return;
  }

  //! Initializations
  vector<double> scores = io_Xt;
  vector<double> norm(p_M, 1.l);
  unsigned int j;
  vector<double> Xel;
  o_ind.clear();
  o_val.clear();

  //! D_ELj[i][j] : equals <d_i, e_l_j> at the end of the j-loop
  vector<double> DELj(p_M * p_K, 0.l);

  //!D_DLj[i][s] = <d_i, d_l_s>
  vector<double> DDLj(p_M * p_K, 0.l);

  //! sparse coordinates : spCoord[j][i] is the coordinate of X[i] on D[j]
  vector<double> spCoord(p_K * p_K, 0.l);

  //! Loop over j
  for (j = 0; j < p_K; j++) {

    //! Initialization
    const unsigned int lj = findMax(scores);
    double* iSc = &spCoord[j * p_K];

    //! Stop if we cannot inverse norm[lj]
    if (norm[lj] < 1e-6) {
      break;
    }

    const long double invNorm = 1.l / (long double) sqrtl(norm[lj]);
    const long double Xelj    = (long double) io_Xt[lj] * invNorm;
    const long double delta   = Xelj * Xelj;

    Xel  .push_back(Xelj);
    o_val.push_back(Xelj); //! The coordinate of x on the last chosen vector
    o_ind.push_back(lj); //! Memorize the chosen index
    normr = (double) ((long double) normr - delta);//! Update of the residual

    //! Gram-Schmidt Algorithm, Update of spCoord
    double* iDELj = &DELj[lj * p_K];
    for (unsigned int i = 0; i < j; i++) {
      iSc[i] = iDELj[i];
    }

    for (unsigned int i = 0; i < j; i++) {
      long double sum = 0.l;
      const double* iScj = &spCoord[j * p_K + i];
      const double* iSci = &spCoord[i * p_K + i];

      for (unsigned int s = 0; s < j - i; s++) {
        sum -= (long double) iSci[s * p_K] * (long double) iScj[s];
      }

      iSc[i] = (double) (sum * invNorm);
    }
    iSc[j] = (double) invNorm;

    //! Break condition
    if (j == p_K - 1 || (normr <= p_eps)) {
      j++;
      break;
    }

    //! Update of DDLj
    const double* iDDt = &i_DDt[lj * p_M];
    double* iDDLj = &DDLj[j];
    for (unsigned int i = 0; i < p_M; i++) {
      iDDLj[i * p_K] = iDDt[i];
    }

    //! Compute the scalar product D[j]_D[lj] and memorize it now
    long double sum    = 0.0l;
    const double* iDlj = &i_D[lj * p_N];
    const double* iDj  = &i_D[j  * p_N];
    for (unsigned int i = 0; i < p_N; i++) {
      sum += (long double) iDlj[i] * (long double) iDj[i];
    }
    DDLj[j * p_K + j] = (double) sum;

    //! Update of io_DELj, io_Xt, norm, and scores.
    for (unsigned int i = 0; i < p_K; i++) {
      long double val = 0;
      const double* iDDLj = &DDLj[i];

      for (unsigned s = 0; s < j + 1; s++) {
        val += (long double) iSc[s] * (long double) iDDLj[s];
      }

      DELj[i * p_K + j] = (double) val;
      io_Xt[i]          = (double) ((long double) io_Xt[i] - Xelj * val);
      norm[i]           = (double) ((long double) norm[i] -  val * val);
      scores[i]         = (double) ((long double) io_Xt[i] * (long double) io_Xt[i] \
                                 / (long double) norm[i]);
    }

    for (unsigned int s = 0; s <= j; s++) {
      scores[o_ind[s]] = double();
    }
  }

  //! Compute the final coordinates of x on the chosen atoms of the dictionary
  for (unsigned int i = 0; i < j; i++) {
    long double sum = 0;
    const double* iXel = &Xel[i];
    const double* iSc  = &spCoord[i * p_K + i];

    for (unsigned s = 0; s < j - i; s++) {
      sum += (long double) iSc[s * p_K] * (long double) iXel[s];
    }

    o_val[i] = (double) sum;
  }
}


//! Find the index of the maximum absolute value of a table
unsigned int findMax(
  std::vector<double> const& i_vec){
  double val       = 0.l;
  unsigned int ind = 0;
  const double* iV = &i_vec[0];

  //! Find the index of the absolute maximum value of i_vec
  for (unsigned int k = 0; k < i_vec.size(); k++) {
    if (val < fabs(iV[k])) {
      val = fabs(iV[k]);
      ind = k;
    }
  }

  return ind;
}



