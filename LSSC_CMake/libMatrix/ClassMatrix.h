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

using namespace std;

class Matrix2
{
public :
  //! Constructor
  Matrix2();
  Matrix2(
    const unsigned int i_row,
    const unsigned int i_col);
  Matrix2(
    const unsigned int i_row,
    const unsigned int i_col,
    const float i_value);

  //! Destructor
  ~Matrix2();

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
  Matrix2& operator=(
    Matrix2 const& i_mat);

  float& operator()(
    const unsigned int p_i,
    const unsigned int p_j);

  float operator() (
    const unsigned int p_i,
    const unsigned int p_j)
    const;

  //! Set the Gram matrix of the input one.
  void setGramMatrix(
    Matrix2 const& i_mat);

  //! Set the transpose of a matrix.
  void setTranspose(
    Matrix2 const& i_mat);

  //! Copy a row of another matrix into this one
  void copyRow(
    Matrix2 const& i_mat,
    const unsigned int p_rowFrom,
    const unsigned int p_rowTo);

  //! Fill a vector with the first elements of the row of a matrix
  void getRow(
    vector<float> &o_row,
    const unsigned int p_row)
    const;

  //! Copy a column of another matrix into this one
  void copyCol(
    Matrix2 const& i_mat,
    const unsigned int p_colFrom,
    const unsigned int p_colTo);

  //! Symmetrize in-place the upper part of the matrix
  void symmetrizeUpperPart(
    const unsigned int p_rowMax);

private :
  unsigned int m_row;
  unsigned int m_col;
  vector<float> m_mat;
};


#endif // CLASSMATRIX_H_INCLUDED