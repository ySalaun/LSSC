/*
 * Copyright (c) 2013, Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIBMATRIX_H_INCLUDED
#define LIBMATRIX_H_INCLUDED

#include <vector>

#include "ClassMatrix.h"

/**
 * @brief Compute A = A + xy' where y' is the transpose of the vector y.
 *
 * @param io_A : matrix of size m x n;
 * @param i_x : vector of size m;
 * @param i_y : vector of size n;
 * @param p_iMax : if > 0, then the addition will be only done for (i,j) < (p_iMax, p_iMax).
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int addXYt(
    Matrix2 &io_A,
    std::vector<float> const& i_x,
    std::vector<float> const& i_y,
    const int p_iMax = -1);


/**
 * @brief Compute C = A * B, classic matrix product.
 *
 * @param i_A : matrix of size m x n;
 * @param i_B : matrix of size n x o;
 * @param o_C : will contain i_A * i_B of size m x o.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int productAB(
    Matrix2 const& i_A,
    Matrix2 const& i_B,
    Matrix2 &o_C);

/**
 * @brief Compute C = A * Bt, classic matrix product.
 *
 * @param i_A : matrix of size m x n;
 * @param i_B : matrix of size o x n;
 * @param o_C : will contain i_A * i_B of size m x o.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int productABt(
    Matrix2 const& i_A,
    Matrix2 const& i_Bt,
    Matrix2 &o_C);


/**
 * @brief Compute y = A * x.
 *
 * @param i_A : matrix of size m x n;
 * @param i_x : vector of size n;
 * @param o_y : will contain A * x of size m;
 * @param p_iMax : if > -1, the product will be done only for the p_iMax first line.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int productAx(
    Matrix2 const& i_A,
    std::vector<float> const& i_x,
    std::vector<float> &o_y,
    const int p_iMax = -1);


/**
 * @brief Compute y = At * x.
 *
 * @param i_A : matrix of size n x m;
 * @param i_x : vector of size n;
 * @param o_y : will contain At * x of size m;
 * @param p_iMax : if > -1, the product will be done only for the p_iMax first line.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int productAtx(
    Matrix2 const& i_At,
    std::vector<float> const& i_x,
    std::vector<float> &o_y,
    const int p_iMax = -1);


/**
 * @brief Compute C = A +/- B.
 *
 * @param i_A : matrix of size m x n;
 * @param i_B : matrix of size m x n;
 * @param o_C : will contain i_A +/- i_B of size m x n;
 * @param p_minus : if true, compute A - B, otherwise compute A + B.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int add(
    Matrix2 const& i_A,
    Matrix2 const& i_B,
    Matrix2 &o_C,
    const bool p_minus = false);


/**
 * @brief Compute z = x +/- y.
 *
 * @param i_x : vector of size m;
 * @param i_y : matrix of size m;
 * @param o_z : will contain i_x +/- i_y of size m;
 * @param p_minus : if true, compute x - y, otherwise compute x + y.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int add(
    std::vector<float> const& i_x,
    std::vector<float> const& i_y,
    std::vector<float> &o_z,
    const bool p_minus = false);


/**
 * @brief Compute the dot product C = A .* B
 *
 * @param i_A : matrix of size m x n;
 * @param i_B : matrix of size m x n;
 * @param o_C : will contain A .* B of size m x n.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int dotProduct(
    Matrix2 const& i_A,
    Matrix2 const& i_B,
    Matrix2 &o_C);


/**
 * @brief Compute the dot product z = x .* y
 *
 * @param i_x : vector of size m;
 * @param i_y : matrix of size m;
 * @param o_z : will contain x .* y of size m or p_iMax;
 * @param p_max : if > -1, the dot product is only done for the p_iMax first values.
 *
 * @return EXIT_FAILURE in case of size problem, EXIT_SUCCESS otherwise.
 **/
int dotProduct(
    std::vector<float> const& i_x,
    std::vector<float> const& i_y,
    std::vector<float> &o_z,
    const int p_max = -1);




#endif // LIBMATRIX_H_INCLUDED
