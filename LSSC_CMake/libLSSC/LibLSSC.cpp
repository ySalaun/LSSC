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
#include "utilities.h"

using namespace std;

#ifdef __linux__
#include <math.h>
#endif

//! train the algorithm with the L1 norm
void trainL1(
  Matrix &io_dict,
  const vector<float> &i_imNoisy,
  const unsigned int p_nPatch,
  const Parameters &p_params){

  //! For convenience
  const unsigned int sP = p_params.sPatch;
  const unsigned int w  = p_params.w;
  const unsigned int k  = p_params.k;
  const unsigned int m  = p_params.m;

  //! matrices used for the dictionnary update
  Matrix A(k, k);
  Matrix B(m, k);

  //! sparse coefficients of the dictionnary
  vector<float> alpha(k, 0.f);

  //! generation of p_nPatch iid patches
  display("- generation of iid patches", p_params, false);
  vector<unsigned int> iidPatches(p_nPatch, 0);
  getRandList(p_nPatch, p_params.nPatch, iidPatches);
  display("...done", p_params);

  for (unsigned int n = 0; n < p_nPatch; n++) {
    vector<float> patch(m, 0.f);
    const unsigned int i = iidPatches[n] / p_params.nRowPatches * sP;
    const unsigned int j = iidPatches[n] % p_params.nRowPatches * sP;
    const float* iI      = &i_imNoisy[i * w + j];
    float* iP            = &patch[0];

    //! initialize patch vector
    for (unsigned int p = 0; p < sP; p++) {
      for (unsigned int q = 0; q < sP; q++) {
        iP[p * sP + q] = iI[p * w + q];
      }
    }

    //! LARS
    display("- LARS", p_params);
    vector<float> alpha;
    computeLars(io_dict, patch, p_params, alpha);
    // TODO debug mode
    if (p_params.debug) {
      for(unsigned int n = 0; n < k; n++){
        cout << alpha[n] << " " ;
      }
      cout << endl;
      cout << "-------------------------------------" << endl;
      for(unsigned int n = 0; sP < k; n++){
        cout << patch[n] << " " ;
      }
      cout << endl;
    }

    A.addXYt(alpha, alpha);
    B.addXYt(patch, alpha);

    //! Update
    display("- dictionary update", p_params, false);
    updateDictionary(io_dict, A, B, p_params);
    display("...done", p_params);
  }
}


//! generation of a list of random number
// TODO Marc : validé
void getRandList(
  const unsigned int p_nb,
  const unsigned int p_maxRange,
  vector<unsigned int> &o_randList){

  //! Initializations
  if (o_randList.size() < p_nb) {
    o_randList.resize(p_nb);
  }

  //! list with all possibilities
  vector<unsigned int> list(p_maxRange);
  for(unsigned int n = 0; n < p_maxRange; n++) {
    list[n] = n;
  }

  //! seed the random function
  srand((unsigned int)(time(time_t(NULL))));

  //! find the p_nb random numbers
  for(unsigned int n = 0; n < p_nb; n++){
    const unsigned int index = rand() % list.size();
    o_randList[n] = list[index];
    list.erase(list.begin() + index);
  }
}


//! lars algorithm that minimizes ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
// TODO Marc : à valider en test réel (algo correspondant)
void computeLars(
  const Matrix &p_dict,
  const vector<float> &p_patch,
  const Parameters &p_params,
  vector<float> &o_alpha){

  //! For convenience
  const unsigned int k  = p_params.k;
  const float lambda    = p_params.reg;

  //! Initialization
  display("-- initialization", p_params);

  //! matrices parameters
  Matrix G    (k, k);
  Matrix Ga   (k, k);
  Matrix Gs   (k, k);
  Matrix invGs(k, k);
  G.productAtB(p_dict, p_dict);

  //! add EPSILON to diagonal elements in G
  // TODO Marc : voir si on peut le supprimer
  for(unsigned int i = 0; i < k; i++){
    G(i, i) += p_params.epsilon;
  }

  //! noisy picture norm
  float normPatch;
  dotProduct(p_patch, p_patch, normPatch);

  //! Code
  o_alpha.assign(k, 0);

  //! active indexes mapping
  vector<int> A(k, -1);

  //! most correlated element
  vector<float> correlation;
  p_dict.productAtx(p_patch, correlation);

  float C = fabs(correlation[0]);
  unsigned int currentIndex = 0;
  for (unsigned int n = 1; n < correlation.size(); n++) {
    if (fabs(correlation[n]) > C) {
      C = fabs(correlation[n]);
      currentIndex = n;
    }
  }

  //! begin by adding a new atom in the code
  bool newAtom = true;

  //! if the norm is lower than the regularization parameter, stop the algorithm
  if (normPatch < lambda) {
    return;
  }

  //! add the new most correlated element at each iteration
  for (unsigned int i = 0; i < k;) {

    // TODO: debug mode
    if(p_params.verbose && p_params.debug) {
      cout << "-- " << i + 1 << "-th iteration,";
      if(currentIndex == -1){
        cout << " there is no current index" << endl;
      }
      else{
        cout << " the current index is: " << currentIndex << endl;
      }
       
    }

    //! NEW ATOM
    if (newAtom) {
      display("-- new atom", p_params);
      A[i] = currentIndex;
      Ga.copyCol(G , currentIndex, i);
      Gs.copyRow(Ga, currentIndex, i);
      Gs.symmetrizeUpperPart(i+1);
      updateGram(invGs, Gs, i);
    }

    //! VARIABLES UPDATES
    //! compute sign vector
    //! sgn = sgn(c)
    vector<float> sgn(i + 1, 1.f);
    for (unsigned int j = 0; j <= i; j++) {
      if (correlation[A[j]] < 0) {
        sgn[j] = -1.f;
      }
    }

    //! compute direction vector
    //! u = invGs * sgn
    vector<float> u;
    invGs.productAx(sgn, u, i + 1, i + 1);

    //! STEP
    display("-- compute step: ", p_params, p_params.debug);
    float gamma = p_params.infinity;

    //! current maximum of correlation
    C = fabs(correlation[A[0]]);


    // TODO: debug to test if c(A[0]) is really the max (should be)
    //  and if all the active indexes have the same correlation
    if (p_params.debug) {
      cout << "current maximum of correlation is: " << C << endl;
      for(unsigned int j = 0; j < k; j++){
        // inactive indexes correlation should be lower than the max of correlation
        if(o_alpha[j] == 0 && j != currentIndex && fabs(correlation[j]) > C){
          cout << "ERROR: wrong correlation for non active index, |c(" << j << ")| = " << fabs(correlation[j]) << " > " << C << endl;
        }
        // active indexes should reach the max of correlation
        else if((o_alpha[j] != 0 || j == currentIndex) &&  fabs(fabs(correlation[j]) - C) > 0.001){ //TODO issue of equality....
          cout << "ERROR: wrong correlation for active index, |c(A[" << j << "])| = " << fabs(correlation[j]) << " != " << C << endl;
        }
        // TODO: what to do when 2 coefficients have the same initial correlation ??
        else if(o_alpha[j] == 0 && j != currentIndex && fabs(correlation[j]) == C){
          cout << "BEWARE: the inactive index " << j << " has reached the maximum of correlation" << endl;
        }
      }
    }

    //! Compute Gamma
    vector<float> GaU;
    Ga.productAx(u, GaU, i + 1, -1); 
    
    for (unsigned int j = 0; j < k; j++) {
      if (o_alpha[j] != 0 || j == currentIndex) {
        continue;
      }

      if (1 + GaU[j] > 0) {
        const float value = (C + correlation[j]) / (1 + GaU[j]);
        if (value < gamma && value > 0) {
          gamma        = value;
          currentIndex = j;
        }
      }

      if (1 - GaU[j] > 0) {
        const float value = (C - correlation[j]) / (1 - GaU[j]);
        if (value < gamma && value > 0) {
          gamma        = value;
          currentIndex = j;
        }
      }
    }

    // TODO: this condition should not appear
    if(gamma < 0){
      cout << "ERROR: gamma should be positive without this condition" << endl;
    }

    //! Compute stepDownDate
    float stepDownDate          = p_params.infinity;
    unsigned int downDateIndex  = 0;
    for (unsigned int j = 0; j <= i; j++) {
      const float ratio = -o_alpha[A[j]] / u[j];
      if (ratio > 0.f && ratio < stepDownDate) {
        stepDownDate  = ratio;
        downDateIndex = j;
      }
    }

    // TODO: this condition should not appear
    if(stepDownDate <= 0){
      cout << "ERROR: stepDownDate should be positive without this condition" << endl;
    }

    //! POLYNOMIAL RESOLUTION
    float a = 0.f;
    float b = 0.f;
    for (unsigned int j = 0; j <= i; j++) {
      const float value = correlation[A[j]];
      a += sgn[j] * u[j];
      b += value  * u[j];
    }
    const float c     = normPatch - lambda;
    const float delta = b * b - a * c;
    float stepMax     = (delta > 0.f ? std::min((b - sqrtf(delta)) / a, C) : p_params.infinity);

    //! FINAL STEP and BREAK
    // TODO: debug mode
    if (p_params.debug) {
      cout << "gamma: " << gamma << " downdate step: " << stepDownDate << " step max: " << stepMax << endl;
      cout << "new current index: " << currentIndex << " with correlation: " << correlation[currentIndex] << " downdate index: " << downDateIndex << endl;
    }
    gamma = std::min(std::min(gamma, stepDownDate), stepMax);
    // break if the step is too high
    if(gamma >= p_params.infinity / 2){
      display(" the step is too high", p_params);
      break;
    }

    for (unsigned int j = 0; j <= i ; j++) {
      o_alpha[A[j]]  += gamma * u  [j];
    }

    for(unsigned int j = 0; j < k; j++){
      correlation[j] -= gamma * GaU[j];
    }

    normPatch += a * gamma * gamma - 2 * b * gamma; // TODO Marc : il faudra peut être passer normPatch, gamma, a, b, et c en double
    if (fabs(gamma)              < p_params.epsilon || 
        fabs(gamma - stepMax)    < p_params.epsilon ||
        normPatch                < p_params.epsilon ||
        fabs(normPatch - lambda) < p_params.epsilon) {
      break;
    }
    if (fabs(gamma - stepDownDate) < p_params.epsilon) {
      display("-- downdate Gram matrices", p_params);
      if(p_params.verbose){
        cout << "for index " << downDateIndex << endl;
      }
      downdateGram(invGs, Gs, Ga, i, downDateIndex);
      A      [downDateIndex]    = -1;
      o_alpha[A[downDateIndex]] = 0;
      newAtom                   = false;
      currentIndex              = -1;
      i--;
    }
    else {
      newAtom                = true;
      i++;
    }
  }
}


//! Update the inverse of the Gram matrix
// TODO Marc : à valider en test réel (algo correspondant)
void updateGram(
  Matrix &io_invGs,
  Matrix const& i_Gs,
  const unsigned int p_iter){
  //! case when the matrix is 1x1
  if (p_iter == 0) {
    // TODO fix the div by zero issue
    // should never happen, only in degenerated cases ==> but to check
    if(i_Gs(0, 0) == 0){
      cerr << "ERROR: the Gs matrix of size 1x1 is 0 and cannot be inverted ==> see updateGram function" << endl;
      return;
    }
    io_invGs(0, 0) = 1.f / i_Gs(0, 0);
  }
  //! General case
  else {
    //! iRow = Gs[i-th row] without the i-th coefficient
    vector<float> iRow;
    i_Gs.getRow(iRow, p_iter, p_iter);

    //! u = Gs^-1 Gs[i-th row]
    vector<float> u;
    io_invGs.productAx(iRow, u, p_iter, p_iter);

    //! sigma = 1/(Gs_ii - u. Gs[i-th row])
    float uGs;
    dotProduct(u, iRow, uGs);
    const float sigma = 1.f / (i_Gs(p_iter, p_iter) - uGs);

    //! Gs^-1_iter,iter = sigma;
    io_invGs(p_iter, p_iter) = sigma;

    //! v = sigma * u
    vector<float> v(p_iter);
    for (unsigned int k = 0; k < p_iter; k++) {
      v[k] = sigma * u[k];
    }

    //! Gs^-1[i-th row] = - sigma u
    io_invGs.setRow(v, p_iter, true, p_iter);

    //! Gs^-1[i-th row] = - sigma u
    io_invGs.setCol(v, p_iter, true, p_iter);

    //! Gs^-1 = Gs^-1 + sigma u u^t
    io_invGs.addXYt(u, v, p_iter);
  }
}

//! Downdate the inverse of the Gram matrix
// TODO Marc : à valider en test réel (algo correspondant)
// TODO issue ==> to debug, after downdate, the correlation is not good !!
void downdateGram(
  Matrix &io_invGs,
  Matrix &io_Gs,
  Matrix &io_Ga,
  const unsigned int p_iter,
  const unsigned int p_critIndex){

  //! Sigma = 1 / Gs^{-1}(critInd, critInd)
  const float sigma = 1.f / io_invGs(p_critIndex, p_critIndex);

  //! u = Gs^-1[critInd-th row]\{Gs^-1_critInd,critInd}
  vector<float> u;
  io_invGs.getRow(u, p_critIndex, p_iter+1);
  u.erase(u.begin() + p_critIndex);

  //! delete critInd-th row and critInd-th column
  io_Gs   .removeRowCol(p_critIndex, p_critIndex, p_iter+1, p_iter+1);
  io_invGs.removeRowCol(p_critIndex, p_critIndex, p_iter+1, p_iter+1);

  //! delete the critInd-th column
  io_Ga   .removeRowCol(-1 , p_critIndex, -1, p_iter+1);

  //! Gs^-1 = Gs^-1 - sigma u u^t
  vector<float> v(p_iter - 1);
  for (unsigned int i = 0; i < p_iter - 1; i++) {
    v[i] = -sigma * u[i];
  }

  io_invGs.addXYt(u, v, p_iter - 1);
}


//! dictionary update algorithm
// TODO Validé !
void updateDictionary(
  Matrix &io_D,
  const Matrix &i_A,
  const Matrix &i_B,
  const Parameters &p_params){

  //! D = Matrix(D.nRow, D.nCol);
  for(unsigned int iter = 0; iter < p_params.updateIteration; iter++) {

    //! u = DA
    Matrix u;
    u.productAB(io_D, i_A);

    //! u = B - u
    u.add(i_B, u, true);

    //! loop on the column uj of u
    for (unsigned int j = 0; j < p_params.k; j++) {

      //! Get the j-th column of u
      vector<float> uj;
      u.getCol(uj, j);

      //! if the A coefficient can be inverted, update the column
      // TODO: in Mairal's code it is 1e-6 and the epsilon used for LARS is 1e-15....
      // TODO Marc: c'est normal, là on travaille en float, et j'imagine que pour le LARS il travaille en double
      // d'ailleurs dans L'ORMP on travaille aussi en double
      if (i_A(j, j) > p_params.epsilon) {
        const float ajjInv = 1.f / i_A(j, j);

        //! uj = uj/ajj + Dj
        for (unsigned int i = 0; i < p_params.m; i++) {
          uj[i] = ajjInv * uj[i] + io_D(i, j);
        }

        //! compute ||uj||
        float norm2 = 0.f;
        for (unsigned int i = 0; i < p_params.m; i++) {
          norm2 += uj[i] * uj[i];
        }

        //! uj = uj / max(||uj||_2, 1)
        if (norm2 > 1.f) {
          const float norm2Inv = 1.f / sqrtf(norm2);
          for (unsigned int i = 0; i < p_params.m; i++) {
            uj[i] *= norm2Inv;
          }
        }
      }
      //! else, fill it with 0
      else{
        for (unsigned int i = 0; i < p_params.m; i++) {
          uj[i] = 0;
        }
      }

      //! D[j] = uj
      io_D.setCol(uj, j);
    }
  }
}





















