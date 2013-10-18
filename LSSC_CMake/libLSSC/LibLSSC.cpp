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

#ifdef __linux__
#include <math.h>
#include <unistd.h>
#endif

using namespace std;

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
    unsigned int mD, pD;
    io_dict.getSize(mD, pD);
    if (mD != m || pD != k) {
      cout << "io_dict has not the good size." << endl;
      return;
    }

    //! matrices used for the dictionnary update
    Matrix A(k, k);
    Matrix B(m, k);

    //! generation of p_nPatch iid patches
    display("- generation of iid patches", p_params, false);
    vector<unsigned int> iidPatches(p_nPatch, 0);
    getRandList(p_nPatch, p_params.nPatch, iidPatches);
    display("...done", p_params);

    for (unsigned int n = 0; n < p_nPatch; n++) {
      vector<float> patch(m, 0.f);
      const unsigned int i = iidPatches[n] / p_params.nColPatches;
      const unsigned int j = iidPatches[n] % p_params.nColPatches;
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

      //! sparse coefficients of the dictionnary
      vector<float> alpha(k, 0.f);

      Matrix patchMatrix(m, 1);
      for(unsigned int i = 0; i < m; i++){
        //patchMatrix(i, 1) = patch[i]; //! NON !!!!!! on commence à 0, pas 1
        patchMatrix(i, 0) = patch[i];
      }
      Matrix alphaMatrix(k, 1);

      computeLarsMairal(patchMatrix, io_dict, alphaMatrix, 1, p_params);

      for(unsigned int i = 0; i < k; i++){
        //alpha[i] = alphaMatrix(i, 1); //! NON !!!!!!!!
        alpha[i] = alphaMatrix(i, 0);
      }
      //computeLars(io_dict, patch, p_params, alpha);
      // TODO debug mode
      /*if (p_params.debug) {
      for(unsigned int n = 0; n < k; n++){
      cout << alpha[n] << " " ;
      }
      cout << endl;
      cout << "-------------------------------------" << endl;
      for(unsigned int n = 0; sP < k; n++){
      cout << patch[n] << " " ;
      }
      cout << endl;
      }*/

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
    int currentIndex = 0;
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
      // TODO: issue seems to come from the fact that the precision is too low....
      if (p_params.debug) {
        bool error = false;
        cout << "current maximum of correlation is: " << C << endl;
        for(int j = 0; j < (int) k; j++){
          // inactive indexes correlation should be lower than the max of correlation
          if(o_alpha[j] == 0 && j != currentIndex && fabs(correlation[j]) > C && fabs(fabs(correlation[j]) - C) > p_params.epsilon){
            error = true;
            cout << "ERROR: wrong correlation for non active index, |c(" << j << ")| = " << fabs(correlation[j]) << " > " << C << endl;
          }
          // active indexes should reach the max of correlation
          else if((o_alpha[j] != 0 || j == currentIndex) &&  fabs(fabs(correlation[j]) - C) > 0.001){ //TODO issue of equality....
            error = true;
            cout << "ERROR: wrong correlation for active index, |c(A[" << j << "])| = " << fabs(correlation[j]) << " != " << C << endl;
          }
          // TODO: what to do when 2 coefficients have the same initial correlation ??
          else if(o_alpha[j] == 0 && j != currentIndex && fabs(correlation[j]) == C){
            cout << "BEWARE: the inactive index " << j << " has reached the maximum of correlation" << endl;
          }
        }
        if(error && false){
        sleep(360);
        }
      }

      //! Compute Gamma
      vector<float> GaU;
      Ga.productAx(u, GaU, i + 1, -1);

      for (unsigned int j = 0; j < k; j++) {
        if (o_alpha[j] != 0 || j == (unsigned int) currentIndex) {
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
      float stepMax     = (delta > 0.f ? min((b - sqrtf(delta)) / a, C) : p_params.infinity);

      //! FINAL STEP and BREAK
      // TODO: debug mode
      if (p_params.debug) {
        cout << "gamma: " << gamma << " downdate step: " << stepDownDate << " step max: " << stepMax << endl;
        cout << "new current index: " << currentIndex << " with correlation: " << correlation[currentIndex] << " downdate index: " << downDateIndex << endl;
        cout << "current code is (in active index order): ";
        for (unsigned int j = 0; j <= i; j++) {
          cout << o_alpha[A[j]] << " ";
        }
        cout << endl;
      }
      gamma = min(min(gamma, stepDownDate), stepMax);
      // break if the step is too high
      if(gamma >= p_params.infinity / 2){
        display(" the step is too high", p_params);
        cout << "gamma: " << gamma << " downdate step: " << stepDownDate << " step max: " << stepMax << endl;
        break;
      }

      for (unsigned int j = 0; j <= i ; j++) {
        o_alpha[A[j]]  += gamma * u  [j];
      }

      if (p_params.debug) {
        cout << "after step, current code is (in active index order): ";
        for (unsigned int j = 0; j <= i; j++) {
          cout << o_alpha[A[j]] << " ";
        }
        cout << endl;
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
          cout << "for index " << A[downDateIndex] << endl;
        }
        //o_alpha[A[downDateIndex]] = 0; // TODO: seems to be non-necessary ==> it is this that causes the downdate
        downdateGram(invGs, Gs, Ga, A, i, downDateIndex);
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
// TODO Marc : validé !!
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
// TODO Marc : validé !!!
// TODO issue ==> to debug, after downdate, the correlation is not good !!
// TODO add change to A to the latex file
void downdateGram(
  Matrix &io_invGs,
  Matrix &io_Gs,
  Matrix &io_Ga,
  vector<int> & io_A,
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

    //! delete the critInd-th element of the active indexes
    // TODO keep or not ?
    for (unsigned int i = p_critIndex; i < p_iter; i++) {
      io_A[i] = io_A[i+1];
    }
    io_A[p_iter] = -1;

    //! Gs^-1 = Gs^-1 - sigma u u^t
    vector<float> v(p_iter);
    for (unsigned int i = 0; i < p_iter; i++) {
      v[i] = -sigma * u[i];
    }

    io_invGs.addXYt(u, v, p_iter);
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



/**
* FROM HERE, IT'S THE CODE of JULIEN MAIRAL REWRITTEN
**/


//! lars algorithm that minimizes ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
void computeLarsMairal(
  const Matrix &i_patch,
  const Matrix &p_dict,
  Matrix &o_alpha,
  const unsigned int p_n,
  const Parameters &p_params){

  // TODO: are these the correct parameters
  //const unsigned int p_m = 0; // useles ???? why here ?
  const unsigned int p_p = p_params.k;
  //const unsigned int p_L = p_params.k; // not really sure ....
  const unsigned int p_L = p_p * p_n;
  //const unsigned int p_lambda = p_params.reg; //! NON !!!! p_lambda est un float, pas un int !
  const float p_lambda = p_params.reg;

  //! Declarations
  Matrix G;
  Matrix DtX;

  //! D'D
  G  .productAtB(p_dict, p_dict);
  DtX.productAtB(p_dict, i_patch);

  //! Add diagonal
  for (unsigned int k = 0; k < p_p; k++) {
    G(k, k) += 1e-10f;
  }

  //! Initializations
  o_alpha.setSize(p_p, p_n);
  for (unsigned int i = 0; i < p_p; i++) {
    for (unsigned int j = 0; j < p_n; j++) {
      o_alpha(i, j) = 0.f;
    }
  }

  //! Compute the norm of patches
  vector<float> norms;
  i_patch.norm2sqCols(norms);

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
  for (unsigned int i = 0; i < p_n; i++) {

    //! Initializations
    vector<int  > ind(p_L);
    vector<float> coeffs(p_L);
    vector<float> DtR(p_p);
    DtX.getCol(DtR, i);
    //cout << "i=" << i << endl;
    //! Call the core LARS function
    //cout << "norms = " << norms[i] << endl;
    coreLarsMairal(DtR, G, coeffs, ind, norms[i], p_lambda, p_L, p_p);
    //cout << "i=" << i << endl;
    //! Get results
    for (unsigned int j = 0; j < p_L; j++) {
      if (ind[j] >= 0) {
        o_alpha(ind[j], i) = coeffs[j];
      }
    }
  }
}


//! Core function for the LARS algorithm.
void coreLarsMairal(
  vector<float> &io_DtR,
  Matrix const& i_G,
  vector<float> &o_coeffs,
  vector<int> &o_ind,
  const float i_normX,
  const float p_lambda,
  const unsigned int p_L,
  const unsigned int p_K){

  //! Initializations
  const unsigned int LL = p_L;
  const unsigned int L = min(LL, p_K);
  const int lengthPath = 4 * L;
  vector<float> u(p_K, 0.f);
  Matrix Gs     (p_L, p_L);
  Matrix Ga     (p_L, p_K);
  Matrix invGs  (p_L, p_L);
  vector<float> values(p_K);
  float normX = i_normX;
  const float infinity = numeric_limits<float>::max();
  for (unsigned int k = 0; k < p_L; k++) {
    o_coeffs[k] = 0.f;
    o_ind   [k] = -1;
  }

  //! Find the most correlated element
  int currentInd = findMax(io_DtR, p_K);
  if (i_normX < p_lambda) {
    return;
  }
  bool newAtom = true;

  int i;
  int iter   = 0;
  float thrs = 0.f;

  for (i = 0; i < (int) L; ++i) {
    iter++;
    if (newAtom) {
      o_ind[i] = currentInd;
      cout << "test a" << endl;
      Ga.copyRow(i_G, o_ind[i], i);
      cout << "test b" << endl;
      const float* iGa = Ga.getRef(i, 0);
      float* iGs       = Gs.getRef(i, 0);
      for (int j = 0; j <= i; j++) {
        iGs[j] = iGa[o_ind[j]];
      }

      //! Update inverse of Gs
      updateGramMairal(Gs, invGs, i);
    }
    cout << "test 1" << endl;
    //! Compute the path direction
    for (int j = 0; j <= i; j++) {
      values[j] = (io_DtR[o_ind[j]] > 0.f ? 1.f : -1.f);
    }
    for (int p = 0; p < i + 1; p++) {
      float val = 0.f;
      const float* iInvGsq = invGs.getRef(0, p);
      const float* iInvGsp = invGs.getRef(p, 0);
      for (int q = 0; q < i + 1; q++) {
        if (q >= p) {
          val += iInvGsq[q * p_L] * values[q];
        }
        else { // Because matrices aren't stored only upper part (they're symmetric)
          val += iInvGsp[q] * values[q];
        }
      }
      u[p] = val;
    }
    cout << "test 2" << endl;
    //! Compute the step on the path
    float step_max = infinity;
    int first_zero = -1;
    for (int j = 0; j <= i; j++) {
      const float ratio = -o_coeffs[j] / u[j];
      if (ratio > 0.f && ratio <= step_max) {
        step_max   = ratio;
        first_zero = j;
      }
    }
    cout << "test 3" << endl;
    const float current_correlation = fabs(io_DtR[o_ind[0]]);
    float step = infinity;
    vector<int> index (p_K, 1);
    for (int p = 0; p <= i; p++) {
      index[o_ind[p]] = -1;
    }
    for (unsigned int p = 0; p < p_K; p++) {
      float val = 0.f;
      const float* iGa = Ga.getRef(0, p);
      for (int q = 0; q < i + 1; q++) {
        val += iGa[q * p_K] * u[q];
      }

      if (index[p] > 0) {
        if (val > -1.f) {
          const float frac = fabs((current_correlation + io_DtR[p]) / (1.f + val));
          if (frac < step) {
            step       = frac;
            currentInd = p;
          }
        }
        if (val < 1.f) {
          const float frac = fabs((current_correlation - io_DtR[p]) / (1.f - val));
          if (frac < step) {
            step       = frac;
            currentInd = p;
          }
        }
      }
      values[p] = val;
    }
    cout << "test 4" << endl;
    //! compute the coefficients of the polynome representing normX^2
    float coeff1 = 0.f;
    for (int j = 0; j <= i; j++) {
      coeff1 += (io_DtR[o_ind[j]] > 0.f ? u[j] : -u[j]);
    }
    float coeff2 = 0.f;
    for (int j = 0; j <= i; j++) {
      coeff2 += io_DtR[o_ind[j]] * u[j];
    }
    float coeff3 = normX - p_lambda;
    cout << "test 5" << endl;
    float step_max2;
    /// L2ERROR
    const float delta = coeff2 * coeff2 - coeff1 * coeff3;
    step_max2 = delta < 0.f ? infinity : (coeff2 - sqrt(delta)) / coeff1;
    step_max2 = min(current_correlation, step_max2);
    step      = min(min(step, step_max2), step_max);

    //! Stop the path
    if (step == infinity) {
      break;
    }

    //! Update coefficients
    for (int p = 0; p < i + 1; p++) {
      o_coeffs[p] += step * u[p];
    }

    //! Update correlations
    for (unsigned int p = 0; p < p_K; p++) {
      io_DtR[p] -= step * values[p];
    }

    //! Update normX
    normX += coeff1 * step * step - 2.f * coeff2 * step;

    //! Update norm1
    thrs += step * coeff1;
    cout << "test 6" << endl;
    //! Choose next action
    if (step == step_max) {
      /// Downdate, remove first_zero
      downdateGramMairal(Gs, invGs, Ga, o_ind, o_coeffs, i, first_zero);

      newAtom = false;
      i = i - 2;
    }
    else {
      newAtom = true;
    }
    cout << "test 7" << endl;
    if (iter >= lengthPath - 1       ||
      fabs(step) < 1e-15             ||
      fabs(step - step_max2) < 1e-15 ||
      normX < 1e-15                  ||
      i == (L - 1)                   ||
      normX - p_lambda < 1e-15       ) {
        break;
    }
  }
}


//! Find the maximum magnitude of a vector.
unsigned int findMax(
  vector<float> const& i_vec,
  const unsigned int p_N){

    //! Initializations
    const float* iV  = &i_vec[0];
    float maxVal     = fabs(iV[0]);
    unsigned int ind = 0;

    for (unsigned int n = 1; n < p_N; n++) {
      if (maxVal < fabs(iV[n])) {
        maxVal = fabs(iV[n]);
        ind = n;
      }
    }

    return ind;
}


//! Update the inverse of the gram matrix Gs = invGs.
void updateGramMairal(
  Matrix const& i_Gs,
  Matrix &io_invGs,
  const unsigned int p_iter){
    //! Trivial case
    if (p_iter == 0) {
      io_invGs(0, 0) = 1.f / i_Gs(0, 0);
    }
    //! Other cases
    else {
      //! Initializations
      vector<float> u(p_iter);
      const float* iGs = ((Matrix) i_Gs).getRef(p_iter, 0);
      cout << "test 1" << endl;
      cout << p_iter << endl;
      for (unsigned int p = 0; p < p_iter; p++) {
        cout << i_Gs(1,1) << endl;
        cout << p << endl;
        cout << "a" << endl;
        unsigned mG, nG;
        i_Gs.getSize(mG, nG);
        cout << "size(i_Gs) = (" << mG << ", " << nG << ")" << endl;
        cout << "p_iter = " << p_iter << ", p = " << p << endl;
        cout << i_Gs(p_iter, p) << endl;
        const float value = iGs[p];
        cout << "b" << endl;
        float* iInvGs = io_invGs.getRef(p, 0);
        cout << "c" << endl;
        for (unsigned int q = 0; q < p_iter; q++) {
          if (p > q) {
            cout << "test <-" << endl;
            u[q] += iInvGs[q] * value;
          }
          else {
            cout << "test ->" << endl;
            //cout << io_invGs(q, p) << endl;
            cout << value << endl;
            u[q] += io_invGs(q, p) * value;
          }
        }
      }
      cout << "test 2" << endl;
      float scalTmp = 0.f;
      for (unsigned int p = 0; p < p_iter; p++) {
        scalTmp += u[p] * iGs[p];
      }
      const float schur = 1.f / (iGs[p_iter] - scalTmp);
      cout << "test 3" << endl;
      float* iInvGs = io_invGs.getRef(p_iter, 0);
      iInvGs[p_iter] = schur;
      for (unsigned int p = 0; p < p_iter; p++) {
        iInvGs[p] = -u[p] * schur;
        float* iInvGsp = io_invGs.getRef(p, 0);

        for (unsigned int q = 0; q <= p; q++) {
          iInvGsp[q] += u[p] * u[q] * schur;
        }
      }
    }
}


//! Downdate the inverse of the Gram matrix.
void downdateGramMairal(
  Matrix &io_Gs,
  Matrix &io_invGs,
  Matrix &io_Ga,
  vector<int> &io_ind,
  vector<float> &io_coeffs,
  const unsigned int p_iter,
  const unsigned int p_firstZero){

    /// Downdate Ga, ind, coeffs
    for (unsigned int p = p_firstZero; p < p_iter; p++) {
      io_Ga.copyRow(io_Ga, p + 1, p);
      io_ind   [p] = io_ind   [p + 1];
      io_coeffs[p] = io_coeffs[p + 1];
    }

    io_ind   [p_iter] = -1;
    io_coeffs[p_iter] = 0.f;

    //! Downdate Gs
    io_Gs.removeRowCol(p_firstZero, p_firstZero, p_iter, p_iter);

    //! Downdate invGs
    vector<float> u(p_iter);
    const float schur = io_invGs(p_firstZero, p_firstZero);
    float* iInvGs = io_invGs.getRef(p_firstZero, 0);
    for (unsigned int p = 0; p < p_firstZero; p++) {
      u[p] = iInvGs[p];
    }

    float* iu = &u[p_firstZero];
    iInvGs = io_invGs.getRef(p_firstZero + 1, p_firstZero);
    for (unsigned int p = 0; p < p_iter - p_firstZero; p++) {
      iu[p] = iInvGs[p];
    }

    for (unsigned int p = p_firstZero; p < p_iter; p++) {
      float* iInvGs1 = io_invGs.getRef(p + 1, 0);
      iInvGs         = io_invGs.getRef(p    , 0);
      for (unsigned int q = 0; q < p_firstZero; q++) {
        iInvGs[q] = iInvGs1[q];
      }
      iInvGs1 = io_invGs.getRef(p + 1, p_firstZero + 1);
      iInvGs  = io_invGs.getRef(p    , p_firstZero    );
      for (unsigned int q = 0; q < p_iter - p_firstZero; q++) {
        iInvGs[q] = iInvGs1[q];
      }
    }

    for (unsigned int p = 0; p < p_iter; p++) {
      iInvGs = io_invGs.getRef(p, 0);
      for (unsigned int q = 0; q <= p; q++) {
        iInvGs[q] -= u[p] * u[q] * (1.f / schur);
      }
    }
}
