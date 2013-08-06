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
 * @file ClassMatrix.cpp
 * @brief Small matrix class.
 *
 * @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "ClassMatrix.h"

#include <iostream>

using namespace std;

//! Constructor
Matrix2::Matrix2() :
    m_row(0),
    m_col(0) {
    m_mat.assign(m_row * m_col, 0.f);
}

Matrix2::Matrix2(
    const unsigned int i_row,
    const unsigned int i_col) :
    m_row(i_row),
    m_col(i_col) {
    m_mat.assign(m_row * m_col, 0.f);

}

Matrix2::Matrix2(
    const unsigned int i_row,
    const unsigned int i_col,
    const float i_value) :
    m_row(i_row),
    m_col(i_col) {
    m_mat.assign(m_row * m_col, i_value);

}

//! Destructor
Matrix2::~Matrix2() {

}

//! Get the size of a matrix
void Matrix2::getSize(
    unsigned int &o_row,
    unsigned int &o_col) const {
    o_row = m_row;
    o_col = m_col;
}

//! Set the size of a matrix
void Matrix2::setSize(
    const unsigned int i_row,
    const unsigned int i_col) {

    m_row = i_row;
    m_col = i_col;
    m_mat.assign(m_row * m_col, 0.f);
}

//! Print a matrix to screen
void Matrix2::print() const {

    for (unsigned int i = 0; i < m_row; i++) {
        const float* iM = &m_mat[i * m_col];
        for (unsigned int j = 0; j < m_col; j++) {
            cout << iM[j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

//! Operator overloading
Matrix2& Matrix2::operator=(
    Matrix2 const& i_mat) {
    m_row = i_mat.m_row;
    m_col = i_mat.m_col;
    m_mat = i_mat.m_mat;

    return *this;
}

float& Matrix2::operator()(
    const unsigned int p_i,
    const unsigned int p_j) {
    return m_mat[p_i * m_col + p_j];
}


//! Get the transpose of a matrix.
void Matrix2::setTranspose(
    Matrix2 const& i_mat) {

    //! Initialization
    m_row = i_mat.m_col;
    m_col = i_mat.m_row;
    m_mat.assign(m_row * m_col, 0.f);

    for (unsigned int i = 0; i < m_row; i++) {
        float*       mM = &m_mat[i * m_col];
        const float* iM = &i_mat.m_mat[i];

        for (unsigned int j = 0; j < m_col; j++) {
            mM[j] = iM[j * m_row];
        }
    }
}


//! Set the Gram matrix of the input one.
void Matrix2::setGramMatrix(
    Matrix2 const& i_mat) {

    //! Initialization
    m_row = i_mat.m_col;
    m_col = i_mat.m_col;
    m_mat.assign(m_row * m_col, 0.f);
    Matrix2 trMat;
    trMat.setTranspose(i_mat);
    const unsigned int N = i_mat.m_row;

    for (unsigned int i = 0; i < m_col; i++) {
        const float* iT = &trMat.m_mat[i * N];
        float* iM       = &m_mat[i * m_col];

        for (unsigned int j = 0; j < m_col; j++) {
            const float* jT = &trMat.m_mat[j * N];
            float sum = 0.f;

            for (unsigned int k = 0; k < m_col; k++) {
                sum += iT[k] * jT[k];
            }
            iM[j] = sum;
        }
    }
}


//! Copy a row of another matrix into this one
void Matrix2::copyRow(
    Matrix2 const& i_mat,
    const unsigned int p_rowFrom,
    const unsigned int p_rowTo) {

    if(m_col != i_mat.m_col) {
        cerr << "copyRow : Error in col size" << endl;
    }
    else {
        const float* iM = &i_mat.m_mat[p_rowFrom * m_col];
        float* mM       = &m_mat[p_rowTo * m_col];

        for (unsigned int j = 0; j < m_col; j++) {
            mM[j] = iM[j];
        }
    }
}


//! Return the row of a matrix
void Matrix2::getRow(
    std::vector<float> &o_row,
    const unsigned int p_row) const{
    if (o_row.size() < m_col) {
        o_row.resize(m_col);
    }

    const float* mM = &m_mat[p_row * m_col];
    float* iR       = &o_row[0];

    for (unsigned int j = 0; j < m_col; j++) {
        iR[j] = mM[j];
    }
}

//! Symmetrize in-place the upper part of the matrix
void Matrix2::symmetrizeUpperPart(
    const unsigned int p_rowMax) {

    if (m_col != m_row) {
        cerr << "symmetrizeUpperPart : Error - the matrix is not square" << endl;
    }

    for (unsigned int i = 0; i < p_rowMax; i++) {
        const float* uM = &m_mat[i * m_col];
        float* dM       = &m_mat[i];

        for (unsigned int j = i + 1; j < p_rowMax; j++) {
            dM[j * m_row] = uM[j];
        }
    }
}


