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

#include "Main\params.h"
#include "utilities.h"

/**
 * @brief Update the dictionnary with the l1 norm
 *
 * @param io_dict : table that contains the dictionary coefficients;
 * @param i_noisy : noisy picture;
 * @param p_nPatch : number of iid patches used for the update;
 * @param params : global parameters.
 *
 * @return
 **/
void trainL1(Matrix &io_dict, const vector<float> &i_noisy, unsigned p_nPatch, const Parameters &params);
vector<int> randPatches(const int nPatch, const int nPatchMax);
void updateDictionary(Matrix &D, const Matrix &A, const Matrix &B, const Parameters &params);
vector<float> lars(const Matrix &p_dict, const vector<float> &p_patch, const Parameters &params);

#endif // LIB_LSSC_H_INCLUDED
