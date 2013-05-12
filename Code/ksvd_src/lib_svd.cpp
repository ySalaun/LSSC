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
 * @file lib_svd.cpp
 * @brief Process a truncated SVD by the power method
 *
 * @author Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include "lib_svd.h"

#include <math.h>
#include <string.h>

/**
 * @brief Process a truncated svd (only the largest singular
 *        value is processed).
 *        /!\ Warning: tX is transposed for convenience
 *
 * @param tX : (n, m) (so X is a m x n matrix);
 * @param U (m) : will contain the new coefficients;
 * @param V (n) : will contain the largest principal vectors.
 *
 * @return S, the largest singular value
 **/
double svd_trunc(mat_t &tX,
                 vec_t &U,
                 vec_t &V)
{
    //! Declarations
    const double epsilon = 10e-6;
    const unsigned max_iter = 100;
    unsigned iter = 0;
    bool go_on = true;
    double S_old = 0;
    double S = 0;
    const unsigned m = U.size();
    const unsigned n = V.size();
    iter_t it_v, it_u, it_x;

    long double norm = 0.0l;
    it_v = V.begin();
    for (unsigned j = 0; j < n; j++, it_v++)
    {
        long double val = 0.0l;
        it_x = tX[j].begin();
        for (unsigned i = 0; i < m; i++, it_x++)
            val += (long double) fabsl(*it_x);
        (*it_v) = (double) val;
        norm += val * val;
    }

    long double s_inv = -1.0l / sqrt(norm);

    for (it_v = V.begin(); it_v < V.end(); it_v++)
        (*it_v) *= (double) s_inv;

    while(iter < max_iter && go_on)
    {
        S_old = S;
        it_u = U.begin();
        norm = 0.0l;
        for (unsigned i = 0; i < m; i++, it_u++)
        {
            long double value = 0.0l;
            it_v = V.begin();
            for (unsigned j = 0; j < n; j++, it_v++)
                value += (long double) tX[j][i] * (long double) (*it_v);
            (*it_u) = (double) value;
            norm += value * value;
        }
        s_inv = 1.0l / sqrt(norm);

        for (it_u = U.begin(); it_u < U.end(); it_u++)
            (*it_u) *= s_inv;

        for (it_v = V.begin(); it_v < V.end(); it_v++)
            (*it_v) = 0.0l;

        it_v = V.begin();
        for (unsigned j = 0; j < n; j++, it_v++)
        {
            it_u = U.begin();
            it_x = tX[j].begin();
            for (unsigned i = 0; i < m; i++, it_u++, it_x++)
                (*it_v) += (*it_x) * (*it_u);
        }


        norm = 0.0l;
        for (it_v = V.begin(); it_v < V.end(); it_v++)
            norm += (*it_v) * (*it_v);

        S = sqrt(norm);
        s_inv = 1.0l / (long double) S;

        for (it_v = V.begin(); it_v < V.end(); it_v++)
            (*it_v) *= s_inv;

        iter++;
        go_on = fabsl(S - S_old) > epsilon * S;
    }

    return S;
}
