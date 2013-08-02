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

/**
 * @brief Update the dictionnary with the l1 norm
 *
 * @param io_dict : table that contains the dictionary coefficients;
 * @param i_noisy : noisy picture;
 * @param p_nPatch : number of iid patches used for the update;
 * @param params : global parameters.
 *
 * @return none.
 **/
void trainL1(
    Matrix &io_dict,
    const std::vector<float> &i_noisy,
    const unsigned int p_nPatch,
    const Parameters &params);


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
    std::vector<unsigned int> &o_randList);


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
    const std::vector<float> &p_patch,
    const Parameters &p_params,
    std::vector<float> &o_alpha);


/**
 * @brief Update the inverse of the Gram matrix.
 *
 * @param io_invGs : Gs^{-1}. Will be updated;
 * @param i_Gs : Gram matrix;
 * @param p_iter : current index.
 *
 * @return none.
 **/
void updateGram(
    Matrix &io_invGs,
    Matrix const& i_Gs,
    const unsigned int p_iter);


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


#endif // LIB_LSSC_H_INCLUDED
