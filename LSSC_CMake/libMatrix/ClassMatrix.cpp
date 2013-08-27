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
Matrix::Matrix() :
m_row(0),
  m_col(0) {
  m_mat.assign(m_row * m_col, 0.f);
}


Matrix::Matrix(
  const unsigned int i_row,
  const unsigned int i_col) :
m_row(i_row),
  m_col(i_col) {
  m_mat.assign(m_row * m_col, 0.f);

}


Matrix::Matrix(
  const unsigned int i_row,
  const unsigned int i_col,
  const float i_value) :
m_row(i_row),
  m_col(i_col) {
  m_mat.assign(m_row * m_col, i_value);

}


//! Destructor
Matrix::~Matrix() {

}


//! Get the size of a matrix
void Matrix::getSize(
  unsigned int &o_row,
  unsigned int &o_col) const {
  o_row = m_row;
  o_col = m_col;
}


//! Set the size of a matrix
void Matrix::setSize(
  const unsigned int i_row,
  const unsigned int i_col) {

  m_row = i_row;
  m_col = i_col;
  m_mat.assign(m_row * m_col, 0.f);
}


//! Print a matrix to screen
void Matrix::print() const {

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
Matrix& Matrix::operator=(
  Matrix const& i_mat) {
  m_row = i_mat.m_row;
  m_col = i_mat.m_col;
  m_mat = i_mat.m_mat;

  return *this;
}


float& Matrix::operator()(
  const unsigned int p_i,
  const unsigned int p_j) {
    return m_mat[p_i * m_col + p_j];
}


float Matrix::operator()(
  const unsigned int p_i,
  const unsigned int p_j)
  const {
    return m_mat[p_i * m_col + p_j];
}


//! Get the transpose of a matrix.
void Matrix::setTranspose(
  Matrix const& i_mat) {

  //! Initialization
  Matrix A(i_mat);
  (*this).setSize(A.m_col, A.m_row);

  for (unsigned int i = 0; i < m_row; i++) {
    float*       iM = &m_mat[i * m_col];
    const float* iA = &A.m_mat[i];

    for (unsigned int j = 0; j < m_col; j++) {
      iM[j] = iA[j * m_row];
    }
  }
}


//! Copy a row of another matrix into this one
void Matrix::copyRow(
  Matrix const& i_mat,
  const unsigned int p_rowFrom,
  const unsigned int p_rowTo) {

  if(m_col != i_mat.m_col) {
    cerr << "copyRow : Error in col size" << endl;
    return;
  }
  else {
    const float* iM = &i_mat.m_mat[p_rowFrom * m_col];
    float* mM       = &m_mat[p_rowTo * m_col];

    for (unsigned int j = 0; j < m_col; j++) {
      mM[j] = iM[j];
    }
  }
}


//! Fill a vector with the p_iMax first elements of the p_row-th row of the current matrix
void Matrix::getRow(
  vector<float> &o_row,
  const unsigned int p_row,
  const int p_iMax) const {

  //! Check size
  if (p_row >= m_row) {
    cerr << "getRow : Error - asked row is too big for the size of the matrix" << endl;
    return;
  }
  if (p_iMax > (int) m_col) {
    cerr << "getRow : Error - too much values wanted" << endl;
    return;
  }

  //! Initializations
  const unsigned int n = (p_iMax > -1 ? (unsigned int) p_iMax : m_col);
  if (o_row.size() != n) {
    o_row.resize(n);
  }
  const float* iM = &m_mat[p_row * m_col];
  float* oR       = &o_row[0];

  for (unsigned int j = 0; j < n; j++) {
    oR[j] = iM[j];
  }
}


//! Fill the i-th row of the current matrix by a vector
void Matrix::setRow(
  std::vector<float> const& i_row,
  const unsigned int p_row,
  const bool p_minus,
  const int p_iMax){

  //! Check size
  const unsigned int n = (p_iMax > -1 ? p_iMax : m_col);
  if (n > m_col || n > i_row.size() || p_row >= m_row) {
    cerr << "setRow : Error - sizes aren't consistent" << endl;
  }

  //! Initializations
  const float sign = (p_minus ? -1.f : 1.f);
  const float* iR  = &i_row[0];
  float* iM        = &m_mat[p_row * m_col];

  //! Fill the column
  for (unsigned int j = 0; j < n; j++) {
    iM[j] = sign * iR[j];
  }
}


//! Fill a vector with the p_iMax first elements of the j-th column of the current matrix
void Matrix::getCol(
  std::vector<float> &o_col,
  const unsigned int p_col,
  const int p_iMax) const {

  //! Check size
  if (p_col >= m_col) {
    cerr << "getCol : Error - asked column is too big for the size of the matrix" << endl;
    return;
  }
  if (p_iMax > (int) m_row) {
    cerr << "getCol : Error - too much values wanted" << endl;
    return;
  }

  //! Initializations
  const unsigned int m = (p_iMax > -1 ? (unsigned int) p_iMax : m_row);
  if (o_col.size() != m) {
    o_col.resize(m);
  }
  const float* iM = &m_mat[p_col];
  float* oC       = &o_col[0];

  for (unsigned int i = 0; i < m; i++) {
    oC[i] = iM[i * m_col];
  }
}


//! Fill the j-th column of the current matrix by a vector
void Matrix::setCol(
  std::vector<float> const& i_col,
  const unsigned int p_col){

  //! Check size
  if (p_col >= m_col) {
    cerr << "setCol : Error - asked column is too big for the size of the matrix" << endl;
    return;
  }
  if (i_col.size() > (int) m_row) {
    cerr << "setCol : Error - too much values" << endl;
    return;
  }

  //! Initializations
  float* iM = &m_mat[p_col];
  const float* iC = &i_col[0];

  for (unsigned int i = 0; i < i_col.size(); i++) {
    iM[i * m_col] = iC[i];
  }
}


//! Copy a column of another matrix into this one
void Matrix::copyCol(
  Matrix const& i_mat,
  const unsigned int p_colFrom,
  const unsigned int p_colTo) {

  if(m_row != i_mat.m_row) {
    cerr << "copyCol : Error in row size" << endl;
    return;
  }
  else {
    for (unsigned int i = 0; i < m_row; i++) {
      m_mat[i * m_col + p_colTo] = i_mat.m_mat[i * m_col + p_colFrom];
    }
  }
}


//! Symmetrize in-place the upper part of the matrix
void Matrix::symmetrizeUpperPart(
  const unsigned int p_rowMax) {

  if (m_col != m_row) {
    cerr << "symmetrizeUpperPart : Error - the matrix is not square" << endl;
    return;
  }

  for (unsigned int i = 0; i < p_rowMax; i++) {
    const float* uM = &m_mat[i * m_col];
    float* dM       = &m_mat[i];

    for (unsigned int j = i + 1; j < p_rowMax; j++) {
      dM[j * m_row] = uM[j];
    }
  }
}


//! Add A and B in the current matrix (sizes must be equals)
void Matrix::add(
  Matrix const& i_A,
  Matrix const& i_B,
  const bool p_minus){
  if (i_A.m_col != i_B.m_col || i_A.m_row != i_B.m_row) {
    cerr << "add : Error - matrix haven't the same size" << endl;
    return;
  }
  if (m_row != i_A.m_row || m_col != i_A.m_col) {
    (*this).setSize(i_A.m_row, i_A.m_col);
  }

  const float sign = (p_minus ? -1.f : 1.f);

  for (unsigned int i = 0; i < m_row; i++) {
    const float* iA = &i_A.m_mat[i * m_col];
    const float* iB = &i_B.m_mat[i * m_col];
    float* iM       = &m_mat[i * m_col];

    for (unsigned int j = 0; j < m_col; j++) {
      iM[j] = iA[j] + sign * iB[j];
    }
  }
}


//! Get the product A * B into the current matrix
void Matrix::productAB(
  Matrix const& i_A,
  Matrix const& i_B,
  const int p_iMax){

  //! Check sizes
  if (i_A.m_col != i_B.m_row) {
    cerr << "productAB : Error - sizes aren't consistent" << endl;
    return;
  }
  if (p_iMax > -1) {
    if (p_iMax > (int) i_A.m_col || p_iMax > (int) i_A.m_row
        || p_iMax > (int) i_B.m_row || p_iMax > (int) i_B.m_col) {
      cerr << "productAB : Error - p_iMax too big" << endl;
    }
  }

  //! Transpose i_B to simplify the calcul
  Matrix Bt;
  Bt.setTranspose(i_B);

  //! Initialization
  const unsigned int m = (p_iMax > -1 ? p_iMax : i_A.m_row);
  const unsigned int n = (p_iMax > -1 ? p_iMax : i_B.m_col);
  (*this).setSize(m, n);

  //! Do the product
  for (unsigned int i = 0; i < m_row; i++) {
    const float* iA = &i_A.m_mat[i * i_A.m_col];
    float* iM       = &m_mat[i * m_col];

    for (unsigned int j = 0; j < m_col; j++) {
      const float* iBt = &Bt.m_mat[j * i_B.m_row];
      float sum = 0.f;

      for (unsigned int k = 0; k < i_A.m_col; k++) {
        sum += iA[k] * iBt[k];
      }

      iM[j] = sum;
    }
  }
}


//! Get the product At * B into the current matrix
void Matrix::productAtB(
  Matrix const& i_A,
  Matrix const& i_B){

  //! Check sizes
  if (i_A.m_row != i_B.m_row) {
    cerr << "productAB : Error - sizes aren't consistent" << endl;
    return;
  }

  //! Transpose i_A and i_B to simplify the calcul
  Matrix Bt, At;
  At.setTranspose(i_A);
  Bt.setTranspose(i_B);

  //! Initialization
  (*this).setSize(i_A.m_col, i_B.m_col);

  for (unsigned int i = 0; i < m_row; i++) {
    float* iM = &m_mat[i * m_col];
    const float* iAt = &At.m_mat[i * At.m_col];

    for (unsigned int j = 0; j < m_col; j++) {
      float sum = 0.f;
      const float* iBt = &Bt.m_mat[j * Bt.m_col];

      for (unsigned int k = 0; k < i_A.m_row; k++) {
        sum += iAt[k] * iBt[k];
      }

      iM[j] = sum;
    }
  }
}


//! Get the dot product A .* B into the current matrix
void Matrix::dotProduct(
  Matrix const& i_A,
  Matrix const& i_B){

  //! Check sizes
  if (i_A.m_row != i_B.m_row || i_A.m_col != i_B.m_col) {
    cerr << "dotProduct : Error - matrix sizes aren't consistent" << endl;
    return;
  }

  //! Initialization
  if (m_row != i_A.m_row || m_col != i_A.m_col) {
    (*this).setSize(i_A.m_row, i_A.m_col);
  }

  //! Get the dot product
  for (unsigned int i = 0; i < m_row; i++) {
    const float* iA = &i_A.m_mat[i * m_col];
    const float* iB = &i_B.m_mat[i * m_col];
    float* iM       = &m_mat[i * m_col];

    for (unsigned int j = 0; j < m_col; j++) {
      iM[j] = iA[j] * iB[j];
    }
  }
}


//! Apply the product y = M * x, where M is the current matrix
void Matrix::productAx(
  std::vector<float> const& i_x,
  std::vector<float> &o_y,
  const int p_iMax) const {

  //! Initialization
  unsigned int m = m_row;
  unsigned int n = m_col;

  //! Check size
  if (p_iMax > -1){
    if(p_iMax > (int) m || p_iMax > (int) n || p_iMax > (int) i_x.size()){
      cerr << "productAx - error : problem of size, p_iMax is too large for the matrix/vector sizes." << endl;
      return;
    }
    m = p_iMax;
    n = p_iMax;
  }
  else if (n != i_x.size()) {
    cerr << "productAx - error : size of matrix do not correspond." << endl;
    return;
  }

  //! Initialization
  if (o_y.size() < m) {
    o_y.resize(m);
  }

  //! Compute the product
  for (unsigned int i = 0; i < m; i++) {
    const float* iM = &m_mat[i * m_col];
    const float* ix = &i_x[0];
    float* oy       = &o_y[0];
    float sum = 0.f;

    for (unsigned int j = 0; j < n; j++) {
      sum += iM[j] * ix[j];
    }
    oy[i] = sum;
  }
}


//! Apply the product y = Mt * x, where M is the current matrix.
void Matrix::productAtx(
  std::vector<float> const& i_x,
  std::vector<float> &o_y,
  const int p_iMax) const {

  //! Initialization
  unsigned int m = m_row;
  unsigned int n = m_col;

  //! Check size
  if (p_iMax > -1){
    if(p_iMax > (int) m || p_iMax > (int) n || p_iMax > (int) i_x.size()){
      cerr << "productAx - error : problem of size, p_iMax is too large for the matrix/vector sizes." << endl;
      return;
    }
    m = p_iMax;
    n = p_iMax;
  }
  else if (m != i_x.size()) {
    cerr << "productAx - error : size of matrix do not correspond." << endl;
    return;
  }

  //! Initialization
  if (o_y.size() < n) {
    o_y.resize(n);
  }

  //! Compute the product
  for (unsigned int j = 0; j < m; j++) {
    const float* iM   = &m_mat[j * m_col];
    float* oy         = &o_y[0];
    const float value = i_x[j];

    for (unsigned int i = 0; i < n; i++) {
      oy[i] += iM[i] * value;
    }
  }
}


//! Compute A = A + xy' where y' is the transpose of the vector y.
void Matrix::addXYt(
  std::vector<float> const& i_x,
  std::vector<float> const& i_y,
  const int p_iMax){

  //! Initialization
  unsigned int m = i_x.size();
  unsigned int n = i_y.size();

  //! Check sizes
  if (p_iMax > -1) {
    if (m < (unsigned int) p_iMax || n < (unsigned int) p_iMax
      || m_row < (unsigned int) p_iMax || m_col < (unsigned int) p_iMax) {
        cerr << "addXYt - error : p_iMax too large for vector/matrix size" << endl;
        return;
    }
    m = (unsigned int) p_iMax;
    n = (unsigned int) p_iMax;
  }
  else {
    if (m_row != m || m_col != n) {
      cerr << "addXYt - error : sizes of vectors does'nt match size of the matrix" << endl;
      return;
    }
  }

  //! Compute the addition
  for (unsigned int i = 0; i < m; i++) {
    float*       iM  = &m_mat[i * m_col];
    const float xVal = i_x[i];
    const float* iy  = &i_y[0];

    for (unsigned int j = 0; j < n; j++) {
      iM[j] += xVal * iy[j];
    }
  }
}


//! Remove the iMax first elements of the i-th row and the jMax first elements of the j-th column
// TODO Marc : validé
void Matrix::removeRowCol(
  const unsigned int p_row,
  const unsigned int p_col,
  const int p_iMax,
  const int p_jMax){

  //! Initializations
  const unsigned int m = (p_iMax > -1 ? std::min(p_iMax + 1, (int) m_row - 1) : m_row - 1);
  const unsigned int n = (p_jMax > -1 ? std::min(p_jMax + 1, (int) m_col - 1) : m_col - 1);

  //! Check sizes
  if (p_iMax > (int) m_row || p_jMax > (int) m_col) {
    cerr << "removeRowCol : Error - indexes are too big for the matrix" << endl;
    return;
  }
  if (p_row > m || p_col > n) {
    cerr << "removeRowCol : Error - line or column number is too big for the matrix" << endl;
    return;
  }

  //! Delete row
  for (unsigned int i = p_row; i < m; i++) {
    float* iM1       = &m_mat[ i      * m_col];
    const float* iM2 = &m_mat[(i + 1) * m_col];

    for (unsigned int j = 0; j < n; j++) {
      iM1[j] = iM2[j];
    }
  }

  //! Delete col
  for (unsigned int i = 0; i < m; i++) {
    float* iM1       = &m_mat[i * m_col + p_col    ];
    const float* iM2 = &m_mat[i * m_col + p_col + 1];

    for (unsigned int j = 0; j < n - p_col; j++) {
      iM1[j] = iM2[j];
    }
  }
}















