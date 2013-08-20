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
#include <limits>
#include <math.h>
#else
#define INFINITY	numeric_limits<float>::max()
#endif
#define EPSILON		1e-15

//! train the algorithm with the L1 norm
void trainL1(
  Matrix &io_dict,
  const vector<float> &i_imNoisy,
  const unsigned int p_nPatch,
  const Parameters &params){

  //! For convenience
  const unsigned int sP = params.sPatch;
  const unsigned int w  = params.w;
  const unsigned int k  = params.k;
  const unsigned int m  = params.m;

  //! matrices used for the dictionnary update
  Matrix A(k, k);
  Matrix B(m, k);

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
    for(unsigned int n = 0; n < k; n++){
      cout << alpha[n] << " " ;
    }
    cout << endl;
    A.addXYt(alpha, alpha);
    B.addXYt(patch, alpha);

    //! Update
    display("- dictionary update", params, false);
    updateDictionary(io_dict, A, B, params);
    display("...done", params);
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
void computeLars(
  const Matrix &p_dict,
  const vector<float> &p_patch,
  const Parameters &p_params,
  vector<float> &o_alpha){

  //! For convenience
  const unsigned int nb = p_params.k;
  const float lambda    = p_params.reg; // TODO : l'appeler lambda comme dans l'article

  //! Initialization
  display("-- initialization", p_params);
  o_alpha.assign(nb, 0);

  // TODO Marc : il n'y a pas d'initialisation d'une variable coeffs de taille nb avec des 0
  vector<float> coeffs(nb, 0.f);

  //! active set indexes
  vector<int> activeIndexes(0);
  // avec, du coup autant initialiser la taille à nb

  //! noisy picture norm
  float norm; // TODO Marc : l'appeler normX en correspondance avec l'article
  dotProduct(p_patch, p_patch, norm);

  //! if the norm is lower than the regularization parameter, stop the algorithm
  if (norm > lambda) { // TODO Marc : à mon avis c'est norm < reg
    return;
  }

  //! most correlated element
  vector<float> correlation;
  p_dict.productAtx(p_patch, correlation);
  float cMax = fabs(correlation[0]);
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
  Matrix G    (nb, nb);
  Matrix Ga   (nb, nb);
  Matrix Gs   (nb, nb);
  Matrix invGs(nb, nb);
  G.productAtB(p_dict, p_dict);

  //! add EPSILON to diagonal elements in G
  // TODO: dunno why but inside Mairal's code
  for(unsigned int k = 0; k < nb; k++){
    G(k, k) += 1e-10; // TODO Marc : ce serait pas += EPSILON ?
  }

  //! add the new most correlated element at each iteration
  for (unsigned int iter = 0; iter < nb; iter++) { // TODO Marc : vu une des conditions de la fin,
    // il vaut mieux mettre iter en int plutôt que unsigned int
    cout << currentIndex << endl;

    //! UPDATE
    if (newAtom) {
      display("-- update", p_params);
      activeIndexes.push_back(currentIndex);
      Ga.copyRow(G , currentIndex, iter);
      Gs.copyCol(Ga, currentIndex, iter); // TODO Marc : est-ce que ça correspond à Gs = D_A * D_A^t? et qu'est-ce que D_A ?
      Gs.symmetrizeUpperPart(iter); // TODO Marc : à quoi ça sert ?
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
    vector<float> Ua; // TODO Marc : dans Algo 3, c'est dénoté u, donc mettre u, ou mettre Ua dans le latex
    invGs.productAx(sgn, Ua, iter + 1);

    //! STEP
    display("-- compute step", p_params);
    float step = INFINITY;

    //! current maximum of correlation
    cMax = correlation[0];
    vector<float> tGaUa;
    Ga.productAtx(Ua, tGaUa, iter + 1); // TODO Marc : ce n'est pas ce qui est décrit dans l'algo 3

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
          currentIndex = j; // TODO Marc : ce n'est pas currentIndex mais criticalIndex d'après Algo 3.
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

    // TODO à la place du code précédent pour la mise à jour de stepMax, criticalInd et gamma,
    // je propose de mettre pour être plus cohérent avec l'article :
    float stepMax = INFINITY;
    unsigned int criticalIndex = 0;
    for (unsigned i = 0; i < iter; i++) {
      const float ratio = -coeffs[i] / Ua[i];
      if (ratio > 0.f && ratio < stepMax) {
        stepMax = ratio;
        criticalIndex = i;
      }
    }
    const float C = correlation[0]; // TODO Marc : pourquoi le 0-ème terme ?
    vector<float> GaUa;
    Ga.productAx(Ua, GaUa, iter + 1);
    float gamma = INFINITY;
    vector<bool> activeIndexesComp(iter + 1);
    for (unsigned int j = 0; j < activeIndexes.size(); j++) {
      activeIndexesComp[activeIndexes[j]] = false;
    }
    for (unsigned int j = 0; j <= iter; j++){
      if (activeIndexesComp[j]) {
        const float val1 = (C + correlation[j]) / (1 + GaUa[j]);
        const float val2 = (C - correlation[j]) / (1 - GaUa[j]);
        if (val1 > 0) {
          if (val2 > 0) {
            gamma = std::min(std::min(val1, val2), gamma);
          }
          else {
            gamma = std::min(val1, gamma);
          }
        }
        else {
          if (val2 > 0) {
            gamma = std::min(val2, gamma);
          }
        }
      }
    }

    // TODO Marc : code pour la partie polynomiale resolution
    //! POLYNOMIAL RESOLUTION
    float aa = 0.f; // TODO Marc : je double le nom des variables pour ne pas entrer en collision avec le code qui suit
    float bb = 0.f;
    for (unsigned int j = 0; j <= iter; j++) {
      const float value = correlation[activeIndexes[j]];
      aa += (value > 0.f ? Ua[j] : -Ua[j]);
      bb += value * Ua[j];
    }
    const float cc = norm - lambda;
    const float Delta = bb * bb - aa * cc;
    const float stepMax2 = std::min((bb - sqrtf(Delta)) / aa, C); // TODO Marc : qu'est-ce qui se passe si Delta < 0 ?

    // TODO Marc : code pour la partie FINAL STEP and BREAK
    //! FINAL STEP and BREAK
    gamma = std::min(std::min(gamma, stepMax), stepMax2);
    for (unsigned int j = 0; j <= iter ; j++) {
      coeffs[j] += gamma * Ua[j];
      correlation[j] -= gamma * GaUa[j];
    }
    norm += aa * gamma * gamma - 2 * bb * gamma; // TODO Marc : il faudra peut être passer norm, gamma, a, b, et c en double
    if (fabs(gamma) < EPSILON || fabs(gamma - stepMax2) < EPSILON || norm < EPSILON || fabs(norm - lambda) < EPSILON) {
      break;
    }
    if (fabs(gamma - stepMax) < EPSILON) {
      downdateGram(invGs, Gs, Ga, activeIndexes, o_alpha, iter, criticalIndex);
      newAtom = false;
      iter -= 2; // TODO Marc : il ne vaudrait mieux pas faire iter = std::max(0, iter - 2) ?
    }
    else {
      newAtom = true;
    }
    // TODO Marc : la suite n'a pas été vérifiée

    //! DOWNDATE STEP
    display("-- compute downdate step", p_params);
    //! if this step is reached, downdate
    float stepDowndate = INFINITY;
    // TODO Marc : WARNING : indexDowndate not used after that !!!
    // TODO Yohann : the function is not defined for now
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
    c = norm - lambda;

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
      norm - lambda < EPSILON	){
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

  // TODO Marc : il reste à trier les coefficients
  //! SORT COEFFICIENTS
  for (unsigned int j = 0; j < activeIndexes.size(); j++) {
    // TODO Marc : je ne suis pas sûr du tout qu'il faille faire ça, je ne comprends pas
    // ce que tu entends par sort(coeffs, A) dans l'article
    o_alpha[activeIndexes[j]] = coeffs[activeIndexes[j]];
  }
}


//! Update the inverse of the Gram matrix
// TODO Marc : à valider
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
    //! iRow = Gs[i-th row]
    vector<float> iRow;
    i_Gs.getRow(iRow, p_iter, p_iter);

    //! u = Gs^-1 Gs[i-th row]
    vector<float> u;
    io_invGs.productAx(iRow, u, p_iter);

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

    //! Gs^-1 = Gs^-1 + sigma u u^t
    io_invGs.addXYt(u, v, p_iter);
  }
}


//! Downdate the inverse of the Gram matrix
void downdateGram(
  Matrix &io_invGs,
  Matrix &io_Gs,
  Matrix &io_Ga,
  vector<int> &io_activeIndexes,
  vector<float> &io_alpha,
  const unsigned int p_iter,
  const unsigned int p_critIndex){

  // TODO Marc: sigma n'est pas utilisé
  //! Sigma = 1 / Gs^{-1}(critInd, critInd)
//  const float sigma = 1.f / io_invGs(p_critIndex, p_critIndex);

  //! u = Gs^-1[critInd-th row]\{Gs^-1_critInd,critInd}
  vector<float> u;
  io_invGs.getRow(u, p_critIndex, p_iter + 1);
  u.erase(u.begin() + p_critIndex);

  //! delete critInd-th row and critInd-th column (only for Gs and Gs^-1)
  for(unsigned int i = p_critIndex; i < p_iter; i++){
    for(unsigned int j = 0; j < p_critIndex; j++){
      io_Ga   (i,j) = io_Ga   (i, j+1);
      io_Gs   (i,j) = io_Gs   (i, j+1);
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
  io_invGs.addXYt(u, u, p_iter-1); // TODO Marc : il manque sigma quelque part non ?
}


//! Downdate the inverse of the Gram matrix
// TODO Marc : à regarder si cette fonction est la bonne
void downdateGramBis(
  Matrix &io_invGs,
  Matrix &io_Gs,
  Matrix &io_Ga,
  vector<int> &io_activeIndexes,
  vector<float> &io_alpha,
  const unsigned int p_iter,
  const unsigned int p_critIndex){

  //! Sigma = 1 / Gs^{-1}(critInd, critInd)
  const float sigma = 1.f / io_invGs(p_critIndex, p_critIndex);

  //! u = Gs^-1[critInd-th row]\{Gs^-1_critInd,critInd}
  vector<float> u;
  io_invGs.getRow(u, p_critIndex, p_iter + 1);
  u.erase(u.begin() + p_critIndex);

  //! delete critInd-th row and critInd-th column
  io_Gs   .removeRowCol(p_critIndex, p_critIndex, p_iter, p_iter);
  io_invGs.removeRowCol(p_critIndex, p_critIndex, p_iter, p_iter);
  io_Ga   .removeRowCol(p_critIndex, p_critIndex, p_iter, p_iter);

  //! delete the critical index coefficients from the output code
  io_alpha[io_activeIndexes[p_critIndex]] = 0; // TODO Marc: why not io_alpha[p_critIndex] = 0 ?

  //! delete the critical index from the active indexes set
  io_activeIndexes.erase(io_activeIndexes.begin() + p_critIndex);

  //! Gs^-1 = Gs^-1 - sigma u u^t
  vector<float> v(p_iter - 1);
  for (unsigned int i = 0; i < p_iter - 1; i++) {
    v[i] = -sigma * u[i];
  }
  io_invGs.addXYt(u, v, p_iter-1);
}


//! dictionary update algorithm
// TODO Validé !
void updateDictionary(
  Matrix &io_D,
  const Matrix &i_A,
  const Matrix &i_B,
  const Parameters &params){

  //! D = Matrix(D.nRow, D.nCol);
  for(unsigned int iter = 0; iter < params.updateIteration; iter++) {

    //! u = DA
    Matrix u;
    u.productAB(io_D, i_A);

    //! u = B - u
    u.add(i_B, u, true);

    //! loop on the column uj of u
    for (unsigned int j = 0; j < params.k; j++) {

      //! Get the j-th column of u
      vector<float> uj;
      u.getCol(uj, j);

      //! if the A coefficient can be inverted, update the column
      // TODO: in Mairal's code it is 1e-6 and the epsilon used for LARS is 1e-15....
      // TODO Marc: c'est normal, là on travaille en float, et j'imagine que pour le LARS il travaille en double
      // d'ailleurs dans L'ORMP on travaille aussi en double
      if (i_A(j, j) > 1e-6) {
        const float ajjInv = 1.f / i_A(j, j);

        //! uj = uj/ajj + Dj
        for (unsigned int i = 0; i < params.m; i++) {
          uj[i] = ajjInv * uj[i] + io_D(i, j);
        }

        //! compute ||uj||
        float norm2 = 0.f;
        for (unsigned int i = 0; i < params.m; i++) {
          norm2 += uj[i] * uj[i];
        }

        //! uj = uj / max(||uj||_2, 1)
        if (norm2 > 1.f) {
          const float norm2Inv = 1.f / sqrtf(norm2);
          for (unsigned int i = 0; i < params.m; i++) {
            uj[i] *= norm2Inv;
          }
        }
      }
      //! else, fill it with 0
      else{
        for (unsigned int i = 0; i < params.m; i++) {
          uj[i] = 0;
        }
      }

      //! D[j] = uj
      io_D.setCol(uj, j);
    }
  }
}





















