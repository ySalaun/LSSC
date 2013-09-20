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

#ifndef CLASSMATRIX_H_INCLUDED
#define CLASSMATRIX_H_INCLUDED

#include <vector>

class Matrix
{
public :
  //! Constructor
  Matrix();
  Matrix(
    const unsigned int i_row,
    const unsigned int i_col);
  Matrix(
    const unsigned int i_row,
    const unsigned int i_col,
    const float i_value);

  //! Destructor
  ~Matrix();

  //! Get the size of a matrix
  void getSize(
    unsigned int &o_row,
    unsigned int &o_col) const;

  //! Set the size of a matrix
  void setSize(
    const unsigned int i_row,
    const unsigned int i_col);

  //! Print a matrix to screen
  void print() const;

  //! Operator overloading
  Matrix& operator=(
    Matrix const& i_mat);

  float& operator()(
    const unsigned int p_i,
    const unsigned int p_j);

  float operator() (
    const unsigned int p_i,
    const unsigned int p_j) const;

  //! Set the transpose of a matrix.
  void setTranspose(
    Matrix const& i_mat);

  //! Copy a row of another matrix into this one
  void copyRow(
    Matrix const& i_mat,
    const unsigned int p_rowFrom,
    const unsigned int p_rowTo);

  //! Fill a vector with the p_iMax first elements of the p_row-th row of the current matrix
  void getRow(
    std::vector<float> &o_row,
    const unsigned int p_row,
    const int p_iMax = -1) const;

  //! Fill the i-th row of the current matrix by a vector
  void setRow(
    std::vector<float> const& i_row,
    const unsigned int p_row,
    const bool p_minus = false,
    const int p_iMax = -1);

  //! Fill a vector with the p_iMax first elements of the j-th column of the current matrix
  void getCol(
    std::vector<float> &o_col,
    const unsigned int p_col,
    const int p_iMax = -1) const;

  //! Fill the j-th column of the current matrix by a vector
  void setCol(
    std::vector<float> const& i_col,
    const unsigned int p_col,
    const bool p_minus = false,
    const int p_iMax = -1);

  //! Copy a column of another matrix into this one
  void copyCol(
    Matrix const& i_mat,
    const unsigned int p_colFrom,
    const unsigned int p_colTo);

  //! Symmetrize in-place the upper part of the matrix
  void symmetrizeUpperPart(
    const unsigned int p_rowMax);

  //! Add A and B in the current matrix (sizes must be equals)
  void add(
    Matrix const& i_A,
    Matrix const& i_B,
    const bool p_minus = false);

  //! Get the product A * B into the current matrix
  void productAB(
    Matrix const& i_A,
    Matrix const& i_B,
    const int p_iMax = -1);

  //! Get the product At * B into the current matrix
  void productAtB(
    Matrix const& i_A,
    Matrix const& i_B);

  //! Get the dot product A .* B into the current matrix
  void dotProduct(
    Matrix const& i_A,
    Matrix const& i_B);

  //! Apply the product y = M * x, where M is the current matrix
  //! if p_max > -1, the product will be done only for the p_jMax x p_iMax submatrice
  //! and the p_iMax subvector.
  void productAx(
    std::vector<float> const& i_x,
    std::vector<float> &o_y,
    const int p_iMax = -1,
    const int p_jMax = -1) const;

  //! Apply the product y = Mt * x, where M is the current matrix
  //! if p_max > -1, the product will be done only for the p_iMax x p_jMax submatrice
  //! and the p_iMax subvector.
  void productAtx(
    std::vector<float> const& i_x,
    std::vector<float> &o_y,
    const int p_iMax = -1,
    const int p_jMax = -1) const;

  //! Compute A = A + xy' where y' is the transpose of the vector y.
  //! if p_iMax > -1, then the addition will be only done for (i,j) < (p_iMax, p_iMax).
  void addXYt(
    std::vector<float> const& i_x,
    std::vector<float> const& i_y,
    const int p_iMax = -1);

  //! Remove the iMax first elements of the i-th row and the jMax first elements
  //! of the j-th column
  void removeRowCol(
    const unsigned int p_row,
    const unsigned int p_col,
    const int p_iMax = -1,
    const int p_jMax = -1);

private :
  unsigned int m_row;
  unsigned int m_col;
  std::vector<float> m_mat;
};


#endif // CLASSMATRIX_H_INCLUDED
