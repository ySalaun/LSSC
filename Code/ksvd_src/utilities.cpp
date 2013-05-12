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
 * @file utilities.cpp
 * @brief Utilities functions
 *
 * @author Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include "utilities.h"
#include "io_png.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;
/**
 * @brief Compute PSNR and RMSE between img_1 and img_2
 *
 * @param img_1 : pointer to an allocated array of pixels.
 * @param img_2 : pointer to an allocated array of pixels.
 * @param psnr  : will contain the PSNR
 * @param rmse  : will contain the RMSE
 * @param size  : size of both images
 *
 * @return none.
 **/
void psnr_rmse(const float *  img_1,
               const float *  img_2,
               float *        psnr,
               float *        rmse,
               const unsigned size)
{
    float tmp = 0.0f;
    for (unsigned k = 0; k < size; k++)
        tmp += (img_1[k] - img_2[k]) * (img_1[k] - img_2[k]);

    (*rmse) = sqrtf(tmp / (float) size);
    (*psnr) = 20.0f * log10f(255.0f / (*rmse));
}
