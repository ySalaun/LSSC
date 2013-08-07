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

#ifdef __linux__
#include <limits>
#include <math.h>
#else
#define INFINITY	numeric_limits<float>::max()
#endif
#define EPSILON		1e-15

//! train the algorithm with the L1 norm
void trainL1(
  Matrix2 &io_dict,
  const vector<float> &i_imNoisy,
  const unsigned int p_nPatch,
  const Parameters &params){

    //! For convenience
    const unsigned int sP = params.sPatch;
    const unsigned int w  = params.w;
    const unsigned int k  = params.k;
    const unsigned int m  = params.m;

    //! matrices used for the dictionnary update
    Matrix2 A(k, k);
    Matrix2 B(m, k);

    //! sparse coefficients of the dictionnary
    vector<float> alpha(k, 0.f);

    //! generation of p_nPatch iid patches
    display("- generation of iid patches", params, false);
    vector<unsigned int> iidPatches(p_nPatch, 0);
    getRandList(p_nPatch, params.nPatch, iidPatches);
    display("...done", params);

    for (unsigned int n = 0; n < p_nPatch; n++) {
      vector<float> patch(m, 0.f);
      const unsigned int i = iidPatches[n] / params.nRowPatches * sP;
      const unsigned int j = iidPatches[n] % params.nRowPatches * sP;
      const float* iI      = &i_imNoisy[i * w + j];
      float* iP            = &patch[0];

      //! initialize patch vector
      for (unsigned int p = 0; p < sP; p++) {
        for (unsigned int q = 0; q < sP; q++) {
          iP[p * sP + q] = iI[p * w + q];
        }
      }

      //! LARS
      display("- LARS", params);
      vector<float> alpha;
      computeLars(io_dict, patch, params, alpha);
      // TODO debug mode
      for(int n = 0; n < k; n++){
        cout << alpha[n] << " " ;
      }
      cout << endl;
      addXYt(A, alpha, alpha);
      addXYt(B, patch, alpha);

      //! Update
      display("- dictionary update", params, false);
      updateDictionary(io_dict, A, B, params);
      display("...done", params);
    }
}

//! generation of a list of random number
void getRandList(
  const unsigned int p_nb,
  const unsigned int p_maxRange,
  vector<unsigned int> &o_randList){
    //! list with all possibilities
    vector<unsigned int> list(p_maxRange);
    for(unsigned int n = 0; n < p_maxRange; n++) {
      list[n] = n;
    }

    //! seed the random function
    srand(unsigned int(time(time_t(NULL))));

    //! find the p_nb random numbers
    for(unsigned int n = 0; n < p_nb; n++){
      const unsigned int index = rand() % list.size();
      o_randList[n] = list[index];
      list.erase(list.begin() + index);
    }
}

//! lars algorithm that minimizes ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
void computeLars(
  const Matrix2 &p_dict,
  const vector<float> &p_patch,
  const Parameters &p_params,
  vector<float> &o_alpha){

    //! For convenience
    const unsigned int nb = p_params.k;
    const float reg       = p_params.reg;

    //! Initialization
    display("-- initialization", p_params);
    o_alpha.assign(nb, 0);

    //! active set indexes
    vector<int> activeIndexes(0);

    //! noisy picture norm
    float norm;
    dotProduct(p_patch, p_patch, norm);

    //! if the norm is lower than the regularization parameter, stop the algorithm
    if (norm > reg) {
      return;
    }

    //! most correlated element
    vector<float> correlation;
    productAtx(p_dict, p_patch, correlation);
    float cMax = abs(correlation[0]);
    unsigned int currentIndex = 0;
    for (unsigned int n = 1; n < correlation.size(); n++) {
      if (fabs(correlation[n]) > cMax) {
        cMax = fabs(correlation[n]);
        currentIndex = n;
      }
    }

    //! begin by adding a new atom in the code
    bool newAtom = true;

    //! matrices parameters
    Matrix2 G    (nb, nb);
    Matrix2 Ga   (nb, nb);
    Matrix2 Gs   (nb, nb);
    Matrix2 invGs(nb, nb);
    productAtB(p_dict, p_dict, G);

    //! add EPSILON to diagonal elements in G
    //! TODO: dunno why but inside Mairal's code
    for(unsigned int k = 0; k < nb; k++){
      G(k, k) += 1e-10;
    }

    //! add the new most correlated element at each iteration
    for (unsigned int iter = 0; iter < nb; iter++) {

      //! UPDATE
      if (newAtom) {
        display("-- update", p_params);
        activeIndexes.push_back(currentIndex);
        Ga.copyRow(G , currentIndex, iter);
        Gs.copyCol(Ga, currentIndex, iter);
        Gs.symmetrizeUpperPart(iter);
        updateGram(invGs, Gs, iter);
      }

      //! VARIABLES UPDATE
      //! compute sign vector
      //! sgn = sgn(c)
      vector<float> sgn(iter + 1, 1.f);
      for (unsigned int j = 0; j <= iter; j++) {
        if (correlation[j] < 0) {
          sgn[j] = -1.f;
        }
      }

      //! compute direction vector
      //! Ua = invGs * sgn
      vector<float> Ua;
      productAx(invGs, sgn, Ua, iter + 1);

      //! STEP
      display("-- compute step", p_params);
      float step = INFINITY;

      //! current maximum of correlation
      cMax = correlation[0];
      vector<float> tGaUa;
      productAtx(Ga, Ua, tGaUa, iter + 1);

      for (unsigned int j = 0; j <= iter; j++) {
        tGaUa[activeIndexes[j]] = INFINITY;
      }
      for (unsigned int j = 0; j < nb; j++) {
        if (tGaUa[j] == INFINITY) {
          continue;
        }

        if (tGaUa[j] > -1) {
          const float value = (cMax + correlation[j]) / (1 + tGaUa[j]);
          if (value < step) {
            step = value;
            currentIndex = j;
          }
        }

        if (tGaUa[j] < 1) {
          const float value = (cMax - correlation[j]) / (1 - tGaUa[j]);
          if (value < step) {
            step = value;
            currentIndex = j;
          }
        }
      }

      //! DOWNDATE STEP
      display("-- compute downdate step", p_params);
      //! if this step is reached, downdate
      float stepDowndate = INFINITY;
      //! Marc: WARNING : indexDowndate not used after that !!!
      //! Yohann:the function is not defined for now
      int indexDowndate;
      for (unsigned int j = 0; j <= iter; j++) {
        const float ratio = -o_alpha[iter] / Ua[iter];
        if (ratio < stepDowndate && ratio >= 0) {
          stepDowndate = ratio;
          indexDowndate = j;
        }
      }

      //! STOPPING STEP
      display("-- compute stopping step", p_params);
      float a, b, c;
      dotProduct(sgn, Ua, a);
      dotProduct(correlation, Ua, b, iter + 1);
      c = norm - reg;

      const float delta = b * b - a * c;
      const float stepStop = min((delta > 0)? (b - sqrtf(delta)) / a : INFINITY, cMax);

      //! TAKE THE STEP
      display("-- take the step ", p_params, false);
      float finalStep = min(step, min(stepDowndate, stepStop));
      cout << finalStep << endl;
      for (unsigned int j = 0; j <= iter; j++) {
        o_alpha[activeIndexes[j]] += finalStep * Ua[j];
      }
      for (unsigned int j = 0; j < nb; j++) {
        correlation[j] -= finalStep * tGaUa[j];
      }
      norm += (a * finalStep - 2 * b) * finalStep;

      //! STOPPING CONDITIONS
      if (fabs(finalStep) < EPSILON ||
        finalStep == stepStop		||
        norm < EPSILON				||
        norm - reg < EPSILON	){
          break;
      }

      //! DOWNDATE
      if (false && finalStep == stepDowndate) {
        display("-- downdate", p_params);
        downdateGram(invGs, Gs, Ga, activeIndexes, o_alpha, iter, indexDowndate);
        iter -= 2;
        newAtom = false;
      }
      else {
        newAtom = true;
      }
    }
}

//! Update the inverse of the Gram matrix
void updateGram(
  Matrix2 &io_invGs,
  Matrix2 const& i_Gs,
  const unsigned int p_iter){

    //! case when the matrix is 1x1
    if (p_iter == 0) {
      // TODO fix the div by zero issue
      // should never happen, only in degenerated cases ==> but to check
      if(i_Gs(0, 0) == 0){
        cerr << "ERROR: the Gs matrix of size 1x1 is 0 and cannot be inverted ==> see updateGram function" << endl;
      }
      io_invGs(0, 0) = 1.f / i_Gs(0, 0);
    }
    //! General case
    else {
      //! iRow = Gs[i-th row]
      vector<float> iRow(p_iter, 0);
      i_Gs.getRow(iRow, p_iter);

      //! u = Gs^-1 Gs[i-th row]
      vector<float> u;
      productAx(io_invGs, iRow, u, p_iter);

      //! sigma = 1/(Gs_ii - u. Gs[i-th row])
      float uGs;
      dotProduct(u, iRow, uGs);
      const float sigma = 1.f / (i_Gs(p_iter, p_iter) - uGs);

      //! Gs^-1_iter,iter = sigma;
      io_invGs(p_iter, p_iter) = sigma;

      //! Gs^-1[i-th row] = - sigma u
      for (unsigned int j = 0; j < p_iter; j++) {
        io_invGs(p_iter, j) = -sigma * u[j];
      }

      //! Gs^-1 = Gs^-1 + sigma u u^t
      addXYt(io_invGs, u, u, p_iter);
    }
}

//! Downdate the inverse of the Gram matrix
void downdateGram(
  Matrix2 &io_invGs,
  Matrix2 &io_Gs,
  Matrix2 &io_Ga,
  vector<int> &io_activeIndexes,
  vector<float> &io_alpha,
  const unsigned int p_iter,
  const unsigned int p_critIndex){

    const float sigma = 1.f / io_invGs(p_critIndex, p_critIndex);

    //! u = Gs^-1[critInd-th row]\{Gs^-1_critInd,critInd}
    vector<float> u(p_iter+1, 0.f);
    io_invGs.getRow(u, p_critIndex);
    u.erase(u.begin() + p_critIndex);

    //! delete critInd-th row and critInd-th column (only for Gs and Gs^-1)
    for(unsigned int i = p_critIndex; i < p_iter; i++){
      for(unsigned int j = 0; j < p_critIndex; j++){
        io_Ga(i,j) = io_Ga(i, j+1);
        io_Gs(i,j) = io_Gs(i, j+1);
        io_invGs(i,j) = io_invGs(i, j+1);
      }
      for(unsigned int j = p_critIndex; j < p_iter; j++){
        io_Ga(i,j) = io_Ga(i, j+1);
        io_Gs(i,j) = io_Gs(i+1, j+1);
        io_invGs(i,j) = io_invGs(i+1, j+1);
      }
    }
    io_Gs.symmetrizeUpperPart(p_iter-1);
    io_invGs.symmetrizeUpperPart(p_iter-1);

    //! delete the critical index coefficients from the output code
    io_alpha[io_activeIndexes[p_critIndex]] = 0;

    //! delete the critical index from the active indexes set
    io_activeIndexes.erase(io_activeIndexes.begin() + p_critIndex);

    //! Gs^-1 = Gs^-1 - sigma u u^t
    addXYt(io_invGs, u, u, p_iter-1);
}

//! dictionary update algorithm
// TODO Marc : travailler sur u^t tout du long, afin d'optimiser les boucles !!
void updateDictionary(
  Matrix2 &io_D,
  const Matrix2 &i_A,
  const Matrix2 &i_B,
  const Parameters &params){

    //! D = Matrix(D.nRow, D.nCol);
    for(unsigned int iter = 0; iter < params.updateIteration; iter++) {

      //! u = DA
      Matrix2 u(params.m, params.k);
      productAB(io_D, i_A, u);

      //! u = B - u
      add(i_B, u, u, true);

      // loop on the column uj of u
      for (unsigned int j = 0; j < params.k; j++) {
        //! if the A coefficient can be inverted, update the column
        // TODO: in Mairal's code it is 1e-6 and the epsilon used for LARS is 1e-15.... 
        if (i_A(j, j) > 1e-6) {
          const float ajjInv = 1.f / i_A(j ,j);
          
          //! uj = uj/ajj + Dj
          for (unsigned int i = 0; i < params.m; i++) {
            u(i, j) = ajjInv * u(i, j) + io_D(i, j);
          }
          
          //! compute ||uj||
          float norm2 = 0.f;
          for (unsigned int i = 0; i < params.m; i++) {
            norm2 += u(i, j) * u(i, j);
          }

          //! uj = uj / max(||uj||_2, 1) 
          if (norm2 > 1.f) {
            const float norm2Inv = 1.f / sqrtf(norm2);
            for (unsigned int i = 0; i < params.m; i++) {
              u(i, j) *= norm2Inv;
            }
          }
        }
        //! else, fill it with 0
        else{
          for (unsigned int i = 0; i < params.m; i++) {
            u(i, j) = 0;
          }
        }
      }

      //! D = u
      io_D = u;
    }
}