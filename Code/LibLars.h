/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LIBLARS_H_INCLUDED
#define LIBLARS_H_INCLUDED

#include <vector>


/**
 * @brief Auxiliary function for lasso
 *      solve min_{ alpha } | | alpha | | _1 s . t . | | x-Dalpha | | _2^2 <= lambda
 *
 * @param io_DtR : D'x
 * @param i_G : Gram Matrix = D'D
 * @param Gs, Ga, invGs, u, coeffs, ind, work, normX seems to be references that are filled in the programm below
 * @param mode ==> 2 = L2ERROR ==> solve min_{ alpha } | | alpha | | _1 s . t . | | x-Dalpha | | _2^2 <= lambda
 * @param pos : positivity constraint
 *
 * @return none.
 **/
void coreLARS2new(
    std::vector<float>& io_DtR
,   std::vector<float> const& i_G //abstract_matrix
//,   std::vector<float>& Gs // Matrix
//,   std::vector<float>& Ga // Matrix
//,   std::vector<float>& invGs // Matrix
//,   std::vector<float>& u
//,   std::vector<float>& coeffs
//,   std::vector<int>& ind
//,   std::vector<float>& work // Matrix
,   float& o_normX
,   const float p_constraint //! const constraint_type mode, ==> not needed cause only L2ERROR mode
,   const bool p_pos
);

/**
 * @brief Find the index of the (absolute) maximum value of a vector.
 *
 * @param i_vec : input vector;
 * @param p_abs : if true, look for the absolute maximum, otherwise the maximum.
 *
 * @return the index of the maximum value of i_vec.
 **/
int findMax(
    std::vector<float> const& i_vec
,   const bool p_abs
);

/**
 * @brief Find the index of the minimum absolute value between the N first values of a vector.
 *
 * @param i_vec : input vector;
 * @param i_N : number of the first values on which we are looking for.
 *
 * @return the index of the minimum absolute value of the i_N first values of i_vec.
 **/
int findMin(
    std::vector<float> const& i_vec
,   const unsigned i_N
);

/**
 * @brief Find the minimum of two values
 *
 * @param i_a : first value;
 * @param i_b : second value.
 *
 * @return min(i_a, i_b).
 **/
int getMin(
    const int i_a
,   const int i_b
);

/**
 * @brief Find the minimum of two values
 *
 * @param i_a : first value;
 * @param i_b : second value.
 *
 * @return min(i_a, i_b).
 **/
unsigned getMin(
    const unsigned i_a
,   const unsigned i_b
);

/**
 * @brief Find the minimum of two values
 *
 * @param i_a : first value;
 * @param i_b : second value.
 *
 * @return min(i_a, i_b).
 **/
float getMin(
    const float i_a
,   const float i_b
);

/**
 * @brief Print a matrix.
 *
 * @param i_mat : matrix to print on screen;
 * @param p_row : number of rows;
 * @param p_col : number of cols;
 * @param p_name : name of the matrix.
 *
 * @return none.
 **/
void printMat(
    std::vector<float> const& i_mat
,   const unsigned p_row
,   const unsigned p_col
,   const char* p_name
);

#endif // LIBLARS_H_INCLUDED
