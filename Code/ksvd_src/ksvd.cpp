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
 * @file ksvd.cpp
 * @brief K-SVD denoising functions
 *
 * @author Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ksvd.h"
#include "lib_ormp.h"
#include "utilities.h"
#include "lib_svd.h"

#include "addnoise_function.h"
#include "io_png.h"

//! For randperm function
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/**
 * @brief Run K-SVD
 *
 * @param sigma: noise value;
 * @param img_noisy : pointer to an allocated array of pixel,
 *					  containing the noisy image;
 * @param img_denoised : pointer to an allocated array which
 *                       will contain the final denoised image;
 * @param width : width of the image;
 * @param height : height of the image;
 * @param chnls : number of channels of the image
 * @param useT : if true, use the trick acceleration in order to
 *               speed up the algorithm by learning the dictionary
 *               on a sub sample of the full patches.
 *
 * @return none.
 **/
void ksvd_ipol(const double   sigma       ,
               double *       img_noisy   ,
               float *        img_denoised,
               const unsigned width       ,
               const unsigned height      ,
               const unsigned chnls       ,
               const bool     useT        )
{
    //! Initializations
    const unsigned N1          = (sigma * 255.0l <= 20 ? 5 :
                                 (sigma * 255.0l <= 60 ? 7 : 9)); //! Size of patches
    const unsigned N2          = 256;    //! size of the dictionary
    const unsigned N1_2        = N1 * N1;
    const unsigned N_iter      = 15; //! Number of iterations
    const unsigned T           = (sigma * 255.0l > 40.0l ? 8 :
                                 (sigma * 255.0l > 10.0l ? 16 : 32));
    const double   gamma       = 5.25l;
    const unsigned num_patches = (width - N1 + 1) * (height - N1 + 1);
    //! C = sqrt(1/(chnls * N1 * N1)*chi2inv(0.93,chnls * N1 * N1));
	const double   C           = (chnls == 1 ? (N1 == 5 ? 1.2017729876383829 : (N1 == 7 ? 1.1456550151825420 : 1.1139195378939404))
											 : (N1 == 5 ? 1.1182997771678573 : (N1 == 7 ? 1.0849724948297015 : 1.0662877194412401)));

	//! Assuming that img_noisy \isin [0, 255]
    for (unsigned k = 0; k < width * height * chnls; k++)
        img_noisy[k] /= 255.0l;

    //! Declarations
    matD_t patches(num_patches, vecD_t(N1_2 * chnls));
    matD_t dictionary(N2, vecD_t(N1_2 * chnls));

    //! Decompose the image into patches
    im2patches(img_noisy, patches, width, height, chnls, N1);

    //! Keep (1 / T) patches to learn the dictionary
    if (useT || floor(num_patches / T) > N2)
    {
        //! Obtain less patches (divided by T)
        const unsigned num_sub_patches = (unsigned) floor(num_patches / T);
        matD_t sub_patches(num_sub_patches, vecD_t(N1_2 * chnls));
        for (unsigned k = 0; k < num_sub_patches; k++)
            sub_patches[k] = patches[k * T];

        //! Obtain the initial dictionary
        obtain_dict(dictionary, sub_patches);

        //! Learn dictionary
		ksvd_process(img_noisy, img_denoised, sub_patches, dictionary, \
					 sigma, N1, N2, N_iter, gamma, C, width, \
					 height, chnls, false);

        //! Denoising
        ksvd_process(img_noisy, img_denoised, patches, dictionary, \
					 sigma, N1, N2, 1, gamma, C, width, \
					 height, chnls, true);
    }
    else
    {
        //! Obtain the initial dictionary
        obtain_dict(dictionary, patches);

        //! Denoising
        ksvd_process(img_noisy, img_denoised, patches, dictionary, sigma, \
                     N1, N2, N_iter, gamma, C, width, height, chnls, true);
    }

	//! Back to the [0, 255] for the image value
    for (unsigned k = 0; k < width * height * chnls; k++)
        img_denoised[k] *= 255.0f;
}

/**
 * @brief Decompose the image into patches
 *
 * @param img : pointer to an allocated array containing
 *              the image to decompose;
 * @param patches : will contain all patches (whose size
 *                  is N x N) of img;
 * @param width : width of img;
 * @param height : height of img;
 * @param chnls : number of channels of img.
 *
 * @return none.
 **/
void im2patches(const double * img,
                matD_t        &patches,
                const unsigned width,
                const unsigned height,
                const unsigned chnls,
                const unsigned N)
{
    //! Declarations
    const unsigned h_p = height - N + 1;
    const unsigned wh_p = (width - N + 1) * h_p;

    matD_t::iterator it_p = patches.begin();
    for (unsigned k = 0; k < wh_p; k++, it_p++)
    {
        unsigned dk = k / h_p + (k % h_p) * width;
        for (unsigned c = 0; c < chnls; c++)
        {
            unsigned dc = c * width * height + dk;
            iterD_t it = (*it_p).begin() + c * N * N;
            for (unsigned p = 0; p < N; p++, dc++)
            for (unsigned q = 0; q < N; q++, it++)
                (*it) = (double) img[dc + q * width];
        }
    }
}

/**
 * @brief Reconstruct images by aggregating patches of size N x N
 *
 * @param patches : matrix containing patches;
 * @param img : pointer to an allocated array
 *              which will contain the denoised image;
 * @param img_ref : pointer to an allocated array which
 *                  contains the noisy image;
 * @param with : width of both images;
 * @param height : height of both images;
 * @param chnls : number of channels of both images;
 * @param lambda : weighting coefficient used for the addition
 *                 of img_ref to img;
 * @param N : size of patches (N x N).
 *
 * @return none.
 **/
void patches2im(matD_t         &patches,
                float *         img,
                const double *  img_ref,
                const unsigned  width,
                const unsigned  height,
                const unsigned  chnls,
                const double    lambda,
                const unsigned  N)
{
    //! Declarations
    const unsigned h_p = height - N + 1;
    const unsigned wh_p = (width - N + 1) * h_p;

    for (unsigned c = 0; c < chnls; c++)
    {
        vecD_t denominator(height * width, 0.0f);
        vecD_t numerator  (height * width, 0.0f);
        matD_t::iterator it_p = patches.begin();

        //! Aggregation
        for (unsigned k = 0; k < wh_p; k++, it_p++)
        {
            unsigned ind = (k % h_p) * width + k / h_p;
            iterD_t it = (*it_p).begin() + c * N * N;
            for (unsigned q = 0; q < N; q++, ind++)
                for (unsigned p = 0; p < N; p++, it++)
                {
                    numerator  [p * width + ind] += (*it);
                    denominator[p * width + ind] ++;
                }
        }

        //! Weighting
        unsigned dc_i = c * width * height;
        iterD_t it_d = denominator.begin();
        iterD_t it_n = numerator.begin();
        for (unsigned k = 0; k < height * width; k++, dc_i++, it_d++, it_n++)
            img[dc_i] = ((*it_n) + lambda * img_ref[dc_i]) / ((*it_d) + lambda);
    }
}

/**
 * @brief Obtain random permutation for a tabular
 *
 * @param perm : will contain a random sequence of [1, ..., N]
 *               where N is the size of perm.
 *
 * @return none.
 **/
void randperm(vecU_t &perm)
{
    //! Initializations
    const unsigned N = perm.size();
    vecU_t tmp(N + 1, 0);
    tmp[1] = 1;
    srand(unsigned(time(NULL)));

    for (unsigned i = 2; i < N + 1; i++)
    {
        unsigned j = rand() % i + 1;
        tmp[i] = tmp[j];
        tmp[j] = i;
    }

    iterU_t it_t = tmp.begin() + 1;
    for (iterU_t it_p = perm.begin(); it_p < perm.end(); it_p++, it_t++)
        (*it_p) = (*it_t) - 1;
}

/**
 * @brief Obtain the initial dictionary, which
 *        its columns are normalized
 *
 * @param dictionary : will contain random patches from patches,
 *                     with its columns normalized;
 * @param patches : contains all patches in the noisy image.
 *
 * @return none.
 **/
void obtain_dict(matD_t         &dictionary,
                 matD_t const&   patches)
{
    //! Declarations
    vecU_t perm(patches.size());

    //! Obtain random indices
    randperm(perm);

    //! Getting the initial random dictionary from patches
    iterU_t it_p = perm.begin();
    for (matD_t::iterator it_d = dictionary.begin(); it_d < dictionary.end();
                                                                    it_d++, it_p++)
        (*it_d) = patches[*it_p];

    //! Normalize column
    double norm;
    for (matD_t::iterator it = dictionary.begin(); it < dictionary.end(); it++)
    {
        norm = 0.0l;
        for (iterD_t it_d = (*it).begin(); it_d < (*it).end(); it_d++)
            norm += (*it_d) * (*it_d);

        norm = 1 / sqrtl(norm);
        for (iterD_t it_d = (*it).begin(); it_d < (*it).end(); it_d++)
            (*it_d) *= norm;
    }
}

/**
 * @brief Apply the whole algorithm of K-SVD
 *
 * @param img_noisy : pointer to an allocated array containing
 *                    the original noisy image;
 * @param img_denoised : pointer to an allocated array which
 *                       will contain the final denoised image;
 * @param patches : matrix containing all patches including in
 *                  img_noisy;
 * @param dictionary : initial random dictionary, which will be
 *                     updated in each iteration of the algo;
 * @param sigma : noise value;
 * @param N1 : size of patches (N1 x N1);
 * @param N2 : number of atoms in the dictionary;
 * @param N_iter : number of iteration;
 * @param gamma : value used in the correction matrix in the
 *                case of color image;
 * @param C : coefficient used for the stopping criteria of
 *            the ORMP;
 * @param width : width of both images;
 * @param height : height of both images;
 * @param chnls : number of channels of both images;
 * @param doReconstruction : if true, do the reconstruction of
 *                           the final denoised image from patches
 *                           (only in the case of the acceleration
 *                            trick).
 *
 * @return none.
 **/
void ksvd_process(const double  *img_noisy,
                  float         *img_denoised,
                  matD_t        &patches,
                  matD_t        &dictionary,
                  const double   sigma,
                  const unsigned N1,
                  const unsigned N2,
                  const unsigned N_iter,
                  const double   gamma,
                  const double   C,
                  const unsigned width,
                  const unsigned height,
                  const unsigned chnls,
				  const bool     doReconstruction)
{
    //! Declarations
    const unsigned N1_2 = N1 * N1;
    const double   corr = (sqrtl(1.0l + gamma) - 1.0l) / ((double) N1_2);
    const double   eps  = ((double) (chnls * N1_2)) * C * C * sigma * sigma;
    const unsigned h_p  = patches[0].size();
    const unsigned w_p  = patches.size();

    //! Mat & Vec initializations
    matD_t dict_ormp   (N2 , vecD_t(h_p, 0.0l));
    matD_t patches_ormp(w_p, vecD_t(h_p, 0.0l));
    matD_t tmp         (h_p, vecD_t(N2, 0.0l));
    vecD_t normCol     (N2);
    matD_t Corr        (h_p, vecD_t(h_p, 0.0l));
    vecD_t U           (h_p);
    vecD_t V;
    matD_t E           (w_p, vecD_t(h_p));

    //! Vector for ORMP
    matD_t ormp_val        (w_p, vecD_t ());
    matU_t ormp_ind        (w_p, vecU_t ());
    matD_t res_ormp        (N2, vecD_t (w_p));
    matU_t omega_table     (N2, vecU_t ());
    vecU_t omega_size_table(N2, 0);
    matD_t alpha           (N2, vecD_t ());

    //! To avoid reallocation of memory
    for (unsigned k = 0; k < w_p; k++)
    {
        ormp_val[k].reserve(N2);
        ormp_ind[k].reserve(N2);
    }

    for (matU_t::iterator it = omega_table.begin(); it < omega_table.end(); it++)
        it->reserve(w_p);

    V.reserve(w_p);

    //! Correcting matrix
    for (unsigned i = 0; i < h_p; i++)
        Corr[i][i] = 1.0l;

    for (unsigned c = 0; c < chnls; c++)
    {
        matD_t::iterator it_Corr = Corr.begin() + N1_2 * c;
        for (unsigned i = 0; i < N1_2; i++, it_Corr++)
        {
            iterD_t it = it_Corr->begin() + N1_2 * c;
            for (unsigned j = 0; j < N1_2; j++, it++)
                (*it) += corr;
        }
    }

    #pragma omp parallel for
        for (unsigned j = 0; j < w_p; j++)
        {
            for (unsigned c = 0; c < chnls; c++)
            {
                iterD_t it_ormp = patches_ormp[j].begin() + c * N1_2;
                iterD_t it = patches[j].begin() + c * N1_2;
                for (unsigned i = 0; i < N1_2; i++, it++, it_ormp++)
                {
                    double val = 0.0l;
                    iterD_t it_tmp = patches[j].begin() + c * N1_2;
                    for (unsigned k = 0; k < N1_2; k++, it_tmp++)
                        val += corr * (*it_tmp);
                    (*it_ormp) = val + (*it);
                }
            }
        }

    //! Big loop
    for (unsigned iter = 0; iter < N_iter; iter++)
    {
        //! Sparse coding
        if (doReconstruction)
            cout << "Final Step :" << endl;
        else
            cout << "Step " << iter + 1 << ":" << endl;
        cout << " - Sparse coding" << endl;

        for (unsigned i = 0; i < h_p; i++)
        {
            iterD_t it_tmp = tmp[i].begin();
            for (unsigned j = 0; j < N2; j++, it_tmp++)
            {
                double val = 0.0l;
                iterD_t it_corr_i = Corr[i].begin();
                iterD_t it_dict_j = dictionary[j].begin();
                for (unsigned k = 0; k < h_p; k++, it_corr_i++, it_dict_j++)
                    val += (*it_corr_i) * (*it_dict_j);
                (*it_tmp) = val * val;
            }
        }

        iterD_t it_normCol = normCol.begin();
        for (unsigned j = 0; j < N2; j++, it_normCol++)
        {
            double val = 0.0l;
            for (unsigned i = 0; i < h_p; i++)
                val += tmp[i][j];
            (*it_normCol) = 1.0l / sqrtl(val);
        }

        for (unsigned i = 0; i < h_p; i++)
        {
            iterD_t it_normCol_j = normCol.begin();
            for (unsigned j = 0; j < N2; j++, it_normCol_j++)
            {
                double val = 0.0l;
                iterD_t it_corr_i  = Corr[i].begin();
                iterD_t it_dict_j = dictionary[j].begin();
                for (unsigned k = 0; k < h_p; k++, it_corr_i++, it_dict_j++)
                    val += (*it_corr_i) * (*it_dict_j);
                dict_ormp[j][i] = val * (*it_normCol_j);
            }
        }

        //! ORMP process
        cout << " - ORMP process" << endl;
        ormp_process(patches_ormp, dict_ormp, ormp_ind, ormp_val, N2, eps);

        for (unsigned i = 0; i < w_p; i++)
        {
            iterU_t it_ind = ormp_ind[i].begin();
            iterD_t it_val = ormp_val[i].begin();
            const unsigned size = ormp_val[i].size();
            for (unsigned j = 0; j < size; j++, it_ind++, it_val++)
                (*it_val) *= normCol[*it_ind];
        }

        //! Residus
        for (unsigned i = 0; i < N2; i++)
        {
            omega_size_table[i] = 0;
            omega_table[i].clear();
            alpha[i].clear();
            for (iterD_t it = res_ormp[i].begin(); it < res_ormp[i].end(); it++)
                *it = 0.0l;
        }

        for (unsigned i = 0; i < w_p; i++)
        {
            iterU_t it_ind = ormp_ind[i].begin();
            iterD_t it_val = ormp_val[i].begin();
            for (unsigned j = 0; j < ormp_val[i].size(); j++, it_ind++, it_val++)
            {
                omega_table[*it_ind].push_back(i);
                omega_size_table[*it_ind]++;
                alpha[*it_ind].push_back(*it_val);
                res_ormp[*it_ind][i] = *it_val;
            }
        }

        //! Dictionary update
        cout << " - Dictionary update" << endl;
        for (unsigned l = 0; l < N2; l++)
        {
            //! Initializations
            const unsigned omega_size = omega_size_table[l];
            iterD_t it_dict_l = dictionary[l].begin();
            iterD_t it_alpha_l = alpha[l].begin();
            iterU_t it_omega_l = omega_table[l].begin();
            U.assign(U.size(), 0.0l);

            if (omega_size > 0)
            {
                iterD_t it_a = it_alpha_l;
                iterU_t it_o = it_omega_l;
                for (unsigned j = 0; j < omega_size; j++, it_a++, it_o++)
                {
                    iterD_t it_d = it_dict_l;
                    iterD_t it_e = E[j].begin();
                    iterD_t it_p = patches[*it_o].begin();
                    for (unsigned i = 0; i < h_p; i++, it_d++, it_e++, it_p++)
                        (*it_e) = (*it_p) + (*it_d) * (*it_a);
                }

                matD_t::iterator it_res = res_ormp.begin();
                for (unsigned k = 0; k < N2; k++, it_res++)
                {
                    iterU_t it_o = it_omega_l;
                    iterD_t it_dict_k = dictionary[k].begin();
                    for (unsigned j = 0; j < omega_size; j++, it_o++)
                    {
                        const double val = (*it_res)[*it_o];
                        if (fabs(val) > 0.0l)
                        {
                            iterD_t it_d = it_dict_k;
                            iterD_t it_e = E[j].begin();
                            for (unsigned i = 0; i < h_p; i++, it_d++, it_e++)
                                (*it_e) -= (*it_d) * val;
                        }
                    }
                }

                //! SVD truncated
                V.resize(omega_size);
                double S = svd_trunc(E, U, V);

                dictionary[l] = U;

                it_a = it_alpha_l;
                iterD_t it_v = V.begin();
                it_o = it_omega_l;
                for (unsigned j = 0; j < omega_size; j++, it_a++, it_v++, it_o++)
                    res_ormp[l][*it_o] = (*it_a) = (*it_v) * S;
            }
        }
        cout << " - done." << endl;
    }

	if (doReconstruction)
	{
		//! Final patches estimations
		cout << " - Aggregation" << endl;
		for (matD_t::iterator it = patches.begin(); it < patches.end(); it++)
			for (iterD_t it_p = it->begin(); it_p < it->end(); it_p++)
				*it_p = 0.0l;

		#pragma omp parallel for
			for (unsigned l = 0; l < N2; l++)
			{
				iterD_t it_a = alpha[l].begin();
				iterU_t it_o = omega_table[l].begin();
				const unsigned omega_size = omega_size_table[l];
				for (unsigned j = 0; j < omega_size; j++, it_a++, it_o++)
				{
					const double val = (*it_a);
					iterD_t it_d = dictionary[l].begin();
					iterD_t it_p = patches[*it_o].begin();
					for (unsigned k = 0; k < h_p; k++, it_d++, it_p++)
						(*it_p) += (*it_d) * val;
				}
			}

        //! First, obtention of the denoised image without weighting with lambda
        patches2im(patches, img_denoised, img_noisy, width, height, chnls, 0.0l, N1);

        //! Second, obtention of lambda from the norm between img_denoised and img_noisy
        double d = 0.0l;
        for (unsigned k = 0; k < width * height * chnls; k++)
            d += (img_denoised[k] - img_noisy[k]) * (img_denoised[k] - img_noisy[k]);
        d /= (height * width * chnls * sigma * sigma);
        const double lambda = abs(sqrtl(d) - 1.0l);

        //! Finally, obtention of the denoised image with lambda
		patches2im(patches, img_denoised, img_noisy, width, height, chnls, lambda, N1);
	}
}


