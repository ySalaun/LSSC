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

/**
 * @file LibMatrix.cpp
 * @brief Small matrix library.
 *
 * @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibMatrix.h"

#include <stdlib.h>
#include <iostream>

using namespace std;


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
    const int p_iMax) {

    //! Initialization
    unsigned int m = i_x.size();
    unsigned int n = i_y.size();
    unsigned int r, c;
    io_A.getSize(r, c);

    //! Check sizes
    if (p_iMax > -1) {
        if (m < (unsigned int) p_iMax || n < (unsigned int) p_iMax
            || r < (unsigned int) p_iMax || c < (unsigned int) p_iMax) {
            cerr << "addXYt - error : p_iMax too large for vector/matrix size" << endl;
            return EXIT_FAILURE;
        }
        m = (unsigned int) p_iMax;
        n = (unsigned int) p_iMax;
    }
    else {
        if (r != m || c != n) {
            cerr << "addXYt - error : sizes of vectors does'nt match size of the matrix" << endl;
            return EXIT_FAILURE;
        }
    }

    //! Compute the addition
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            io_A(i, j) += i_x[i] * i_y[j];
        }
    }

    return EXIT_SUCCESS;
}


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
    Matrix2 &o_C) {

    //! Initialization
    unsigned int m, n, p, q;
    i_A.getSize(m, n);
    i_B.getSize(p, q);
    o_C.setSize(m, q);

    if (n != p) {
        cerr << "productAB - error : size of matrix aren't consistant" << endl;
        return EXIT_FAILURE;
    }

    //! Transpose i_B to simplify the calcul
    Matrix2 Bt;
    Bt.setTranspose(i_B);

    //! Compute the matrix product
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < q; j++) {
            float sum = 0.f;
            for (unsigned int k = 0; k < n; k++) {
                sum += ((Matrix2) i_A)(i, k) * Bt(j, k);
            }
            o_C(i, j) = sum;
        }
    }

    return EXIT_SUCCESS;
}


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
    Matrix2 &o_C) {

    //! Initialization
    unsigned int m, n, p, q;
    i_A .getSize(m, n);
    i_Bt.getSize(q, p);
    o_C .setSize(m, q);

    if (n != p) {
        cerr << "productAB - error : size of matrix aren't consistent" << endl;
        return EXIT_FAILURE;
    }

    //! Compute the matrix product
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < q; j++) {
            float sum = 0.f;
            for (unsigned int k = 0; k < n; k++) {
                sum += ((Matrix2) i_A)(i, k) * ((Matrix2) i_Bt)(j, k);
            }
            o_C(i, j) = sum;
        }
    }

    return EXIT_SUCCESS;
}


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
    const int p_iMax) {

    //! Initialization
    unsigned int m, n;
    i_A.getSize(m, n);

    //! Check size
    if ((p_iMax > -1 && p_iMax > (int) m) || n != i_x.size()) {
        cerr << "productAx - error : problem of size." << endl;
        return EXIT_FAILURE;
    }

    //! Initialization
    if (p_iMax > -1) {
        m = p_iMax;
    }
    if (o_y.size() < m) {
        o_y.resize(m);
    }

    //! Compute the product
    for (unsigned int i = 0; i < m; i++) {
        float sum = 0.f;
        for (unsigned int j = 0; j < n; j++) {
            sum += ((Matrix2) i_A)(i, j) * i_x[j];
        }
        o_y[i] = sum;
    }

    return EXIT_SUCCESS;
}


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
    const int p_iMax) {

    //! Initialization
    unsigned int m, n;
    i_At.getSize(n, m);

    //! Check size
    if ((p_iMax > -1 && p_iMax > (int) m) || n != i_x.size()) {
        cerr << "productAx - error : problem of size." << endl;
        return EXIT_FAILURE;
    }

    //! Initialization
    if (p_iMax > -1) {
        m = p_iMax;
    }
    o_y.assign(m, 0.f);

    //! Compute the product
    for (unsigned int j = 0; j < n; j++) {
        const float value = i_x[j];
        float* oY = &o_y[0];

        for (unsigned int i = 0; i < m; i++) {
            oY[i] += ((Matrix2) i_At)(j, i) * value;
        }
    }

    return EXIT_SUCCESS;
}


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
    const bool p_minus) {

    //! Initialization
    unsigned int m, n, p, q;
    i_A.getSize(m, n);
    i_B.getSize(p, q);
    const float sign = (p_minus ? -1.f : 1.f);

    //! Check size
    if (m != p || n != q) {
        cerr << "addAB - error : size of matrices aren't consistent" << endl;
    }
    o_C.setSize(m, n);

    //! Compute the addition
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            o_C(i, j) = ((Matrix2) i_A)(i, j) + sign * ((Matrix2) i_B)(i, j);
        }
    }

    return EXIT_SUCCESS;
}


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
    const bool p_minus) {

    //! Initialization
    const unsigned int m = i_x.size();
    const float sign     = (p_minus ? -1.f : 1.f);
    if (o_z.size() != m) {
        o_z.resize(m);
    }

    //! Check size
    if (i_y.size() != m) {
        cerr << "addxy - error : vector sizes aren't consistent" << endl;
        return EXIT_FAILURE;
    }

    //! Compute the addition
    for (unsigned int i = 0; i < m; i++) {
        o_z[i] = i_x[i] + sign * i_y[i];
    }

    return EXIT_SUCCESS;
}


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
    Matrix2 &o_C) {

    //! Initialization
    unsigned int m, n, p, q;
    i_A.getSize(m, n);
    i_B.getSize(p, q);
    o_C.setSize(m, n);

    //! Check size
    if (m != p || n != q) {
        cerr << "doPoduct - error : matrix sizes aren't consistent" << endl;
        return EXIT_FAILURE;
    }

    //! Compute the dot product
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            o_C(i, j) = ((Matrix2) i_A)(i, j) * ((Matrix2) i_B)(i, j);
        }
    }

    return EXIT_FAILURE;
}


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
    const int p_max) {

    //! Initialization
    const unsigned int m = (p_max > -1 ? p_max : i_x.size());
    if (o_z.size() != m) {
        o_z.resize(m);
    }

    //! Check size
    if (i_x.size() != i_y.size() || p_max > (int) i_x.size()) {
        cerr << "dotProduct - error : vector sizes aren't consistent" << endl;
        return EXIT_FAILURE;
    }

    //! Compute the dot product
    for (unsigned int i = 0; i < m; i++) {
        o_z[i] = i_x[i] * i_y[i];
    }

    return EXIT_SUCCESS;
}

