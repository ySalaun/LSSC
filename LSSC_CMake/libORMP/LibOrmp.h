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
#ifndef LIBORMP_H_INCLUDED
#define LIBORMP_H_INCLUDED

#include <vector>


/**
 * @brief Compute the Orthogonal Recursive Matching Pursuit over i_X.
 *
 * @param i_X : vector of size (L x N) to represent by i_D;
 * @param i_D : dictionary of size (M x N);
 * @param o_indV : will contain p_L lists of indexes of atoms of i_D (one for each rows of i_X)
 * @param o_valV : will contain p_L lists of values corresponding to the index of atoms;
 * @param p_eps : break condition;
 * @param p_L : number of rows of i_X;
 * @param p_M : number of atoms in the dictionary;
 * @param p_N : number of columns of i_X and i_D (i.e. size of patches).
 *
 * @return none.
 **/
void computeORMP(
  std::vector<double> const& i_X,
  std::vector<double> const& i_D,
  std::vector<std::vector<unsigned int> > &o_indV,
  std::vector<std::vector<double> > &o_valV,
  const double p_eps,
  const unsigned int p_L,
  const unsigned int p_M,
  const unsigned int p_N);


/**
 * @brief Sub routine of computeORMP.
 *
 * @param i_D : dictionary of size (p_M x p_N);
 * @param i_DDt : precomputed matrix D * D' of size (p_M x p_M);
 * @param io_Xt : Xt[i] = <X, t_i^j> = <X, D_i> - \sum_{p=0}^{j-1} <D_i, e_l_p> <X, e_l_p>;
 * @param o_ind : will contain index of atoms of D at the end of the j-loop;
 * @param o_val : \alpha_l_q = \sum_{p=q}^j <x, e_l_p> a_{pq} : coordinate
 *                 of x on d_l_q at the end of the j-loop;
 * @param p_eps : break condition;
 * @param p_normr : ||x||^2 - \sum_{p=0}^j <x, e_l_p>^2;
 * @param p_M : number of elements of the dictionary;
 * @param p_N : size of each element of the dictionary;
 * @param p_K : maximum number of iterations.
 *
 * @return none.
 **/
void coreORMP(
  std::vector<double> const& i_D,
  std::vector<double> const& i_DDt,
  std::vector<double> &io_Xt,
  std::vector<unsigned int> &o_ind,
  std::vector<double> &o_val,
  const double p_eps,
  const double p_normr,
  const unsigned int p_M,
  const unsigned int p_N,
  const unsigned int p_K);


/**
 * @brief Find the index of the maximum absolute value of a table
 *
 * @param i_vec : table of values.
 *
 * @return : return the index of the maximum absolute
 *          value of V.
 **/
unsigned int findMax(
  std::vector<double> const& i_vec);


#endif // LIBORMP_H_INCLUDED
