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

#ifndef LIB_LSSC_H_INCLUDED
#define LIB_LSSC_H_INCLUDED

#include <stdlib.h>
#include <vector>
#include <time.h>

#include "utilities.h"
#ifdef __linux__
#include "../libMatrix/LibMatrix.h"
#else
#include "libMatrix\LibMatrix.h"
#endif

/**
* @brief Update the dictionnary with the l1 norm
*
* @param io_dict : table that contains the dictionary coefficients;
* @param i_noisy : noisy picture;
* @param p_nPatch : number of iid patches used for the update;
* @param p_params : global parameters.
*
* @return none.
**/
void trainL1(
  Matrix &io_dict,
  const vector<float> &i_noisy,
  const unsigned int p_nPatch,
  const Parameters &p_params);


/**
* @brief Get p_nb random number in a p_maxRange range. p_maxRange must be
*  superior to p_nb.
*
* @param p_nb : number of random number wanted;
* @param p_maxRange : range;
* @param o_randList : will contains the p_nb random numbers.
*
* @return none.
**/
void getRandList(
  const unsigned int p_nb,
  const unsigned int p_maxRange,
  vector<unsigned int> &o_randList);


/**
* @brief Compute the LARS algorithm that minimizes
*      ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
*
* @param p_dict : dictionary;
* @param p_patch : current patch;
* @param p_params : see Parameters;
* @param o_alpha : will contain the minimized coefficients.
*
* @return none.
**/
void computeLars(
  const Matrix &p_dict,
  const vector<float> &p_patch,
  const Parameters &p_params,
  vector<float> &o_alpha);


/**
* @brief Update the inverse of the Gram matrix.
*
* @param io_invGs : Gs^{-1} of size p_iter \times p_iter. Will be updated
*       of size p_iter+1 \times p_iter+1.
* @param i_Gs : Gram matrix of size p_iter+1 \times p_iter+1
* @param p_iter : current index.
*
* @return none.
**/
void updateGram(
  Matrix &io_invGs,
  const Matrix & i_Gs,
  const unsigned int p_iter);

/**
* @brief  Downdate the inverse of the Gram matrix
*
* @param io_invGs : Gs^{-1}. Will be downdated. Gs of size p_iter+1 \times p_iter+1, will
*       be of size p_iter \times p_iter;
* @param io_Gs : Gram matrix. Will be downdated;
* @param io_Ga : pseudo-Gram matrix. Will be downdated;
* @param io_A  : active index list. Will be downdated;
* @param p_iter : current index;
* @param p_critIndex : critical index that shows where the downdate has to be done.
*
* @return none.
**/
void downdateGram(
  Matrix &io_invGs,
  Matrix &io_Gs,
  Matrix &io_Ga,
  vector<int> & io_A,
  const unsigned int p_iter,
  const unsigned int p_critIndex);


/**
* @brief Update the dictionary.
*
* @param io_D : dictionary which will be updated;
* @param i_A : ?
* @param i_B : ?
* @param p_params : see Parameters.
*
* @return none.
**/
void updateDictionary(
  Matrix &io_D,
  const Matrix &i_A,
  const Matrix &i_B,
  const Parameters &params);


/**
 * FROM HERE, IT'S THE CODE of JULIEN MAIRAL REWRITTEN
 **/



/**
 * @brief Compute the LARS algorithm that minimizes
 *      ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
 *
 * @param p_patch : (m x n) matrix containing n patches of size m;
 * @param p_dict : (m x p) matrix (dictionary) containing p elements of length m;
 * @param o_alpha : (p x n) will contain the minimized coefficients;
 * @param p_m : size of a patch;
 * @param p_n : number of patches;
 * @param p_p : number of elements of the dictionary;
 * @param p_L : maximum number of coefficients;
 * @param p_lambda : constraint;

 *
 * @return none.
 **/
void computeLarsMairal(
  const Matrix &i_patch,
  const Matrix &p_dict,
  Matrix &o_alpha,
  const unsigned int p_m,
  const unsigned int p_n,
  const unsigned int p_p,
  const unsigned int p_L,
  const float p_lambda);


/**
 * @brief Core function for the LARS algorithm.
 *
 * @param io_DtR :
 * @param i_G :
 * @param io_Gs :
 * @param io_Ga :
 * @param io_invGs :
 * @param io_u :
 * @param o_coeffs :
 * @param o_ind :
 * @param io_work :
 * @param io_normX :
 * @param p_lambda :
 * @param p_L :
 * @param p_K :
 *
 * @return none.
 **/
void coreLarsMairal(
  std::vector<float> &io_DtR,
  Matrix const& i_G,
  std::vector<float> &o_coeffs,
  std::vector<int> &o_ind,
  const float i_normX,
  const float p_lambda,
  const unsigned int p_L,
  const unsigned int p_K);


/**
 * @brief Find the maximum magnitude of a vector.
 *
 * @param i_vec : input values;
 * @param p_N : number of values.
 *
 * @return the corresponding index to the biggest magnitude.
 **/
unsigned int findMax(
  std::vector<float> const& i_vec,
  const unsigned int p_N);


/**
 * @brief Update the inverse of the gram matrix Gs = invGs.
 *
 * @param i_Gs : gram matrix;
 * @param io_invGs : inverse to update according to i_Gs;
 * @param p_iter : current index.
 *
 * @return none.
 **/
void updateGramMairal(
  Matrix const& i_Gs,
  Matrix &io_invGs,
  const unsigned int p_iter);


/**
 * @brief Downdate the inverse of the Gram matrix.
 *
 * @param io_Gs : Gram matrix. Will be downdated;
 * @param io_invGs : inverse of the Gram matrix. Will be downdated;
 * @param io_Ga : pseudo-Gram matrix. Will be downdated;
 * @param io_ind : index of coefficients. Will be downdated;
 * @param io_coeffs : coefficients. Will be downdated;
 * @param p_iter : current index;
 * @param p_firstZero : critical index where the downdate has to be done.
 *
 * @return none.
 **/
void downdateGramMairal(
  Matrix &io_Gs,
  Matrix &io_invGs,
  Matrix &io_Ga,
  std::vector<int> &io_ind,
  std::vector<float> &io_coeffs,
  const unsigned int p_iter,
  const unsigned int p_firstZero);


#endif // LIB_LSSC_H_INCLUDED
