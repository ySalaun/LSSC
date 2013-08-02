/*
 * Copyright (c) 2011, Marc Lebrun <marc.lebrun.ik@gmail.com>
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
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"
#include "mt19937ar.h"
#ifdef __linux__
    #include "../libIO/io_png.h"
#else
    #include "libIO/io_png.h"
#endif

#define _USE_MATH_DEFINES
#include <cmath>			// for M_PI
#ifndef __linux__
    #include <process.h>		// for getpid() function
#endif
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

 using namespace std;

/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
int loadImage(
    char* p_name
,   std::vector<float> &o_im
,   ImageSize &o_imSize
,   const bool p_verbose
){
    //! read input image
    if (p_verbose) {
        cout << endl << "Read input image...";
    }
	float *imTmp = NULL;
	size_t w, h, c;
	imTmp = read_png_f32(p_name, &w, &h, &c);
	if (!imTmp) {
		cout << "error :: " << p_name << " not found or not a correct png image" << endl;
		return EXIT_FAILURE;
	}
	if (p_verbose) {
        cout << "done." << endl;
	}

	//! test if image is really a color image and exclude the alpha channel
	if (c > 2) {
	    unsigned k = 0;
	    while (k < w * h && imTmp[k] == imTmp[w * h + k] && imTmp[k] == imTmp[2 * w * h + k]) {
            k++;
	    }
        c = (k == w * h ? 1 : 3);
	}

	//! Some image informations
	if (p_verbose) {
        cout << "image size :" << endl;
        cout << " - width          = " << w << endl;
        cout << " - height         = " << h << endl;
        cout << " - nb of channels = " << c << endl;
	}

	//! Initializations
	o_imSize.width      = w;
	o_imSize.height     = h;
	o_imSize.nChannels  = c;
	o_imSize.wh         = w * h;
	o_imSize.whc        = w * h * c;
	o_im.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
        o_im[k] = imTmp[k];

    return EXIT_SUCCESS;
}

/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 * @param p_min, p_max : range of data (usually [0, 255]).
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int saveImage(
    char* p_name
,   std::vector<float> const& i_im
,   const ImageSize &p_imSize
,   const float p_min
,   const float p_max
){
    //! Allocate Memory
    float* imTmp = new float[p_imSize.whc];

    //! Check for boundary problems
    for (unsigned k = 0; k < p_imSize.whc; k++) {
        imTmp[k] = (i_im[k] < p_min ? p_min : (i_im[k] > p_max ? p_max : i_im[k]));
    }

    if (write_png_f32(p_name, imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels) != 0) {
        cout << "... failed to save png image " << p_name << endl;
        return EXIT_FAILURE;
    }

    //! Free Memory
    delete[] imTmp;

    return EXIT_SUCCESS;
}

/**
 * @brief add noise to img.
 *
 * @param i_im : original noise-free image;
 * @param o_imNoisy = im + noise;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
void addNoise(
    std::vector<float> const& i_im
,   std::vector<float> &o_imNoisy
,   const float p_sigma
,   const bool p_verbose
){
    if (p_verbose) {
        cout << "Add noise [sigma = " << p_sigma << "] ...";
    }

	//! Initialization
    o_imNoisy = i_im;
#ifdef __linux__
    mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid());
#else
    mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) _getpid());
#endif

    //! Add noise
    for (unsigned k = 0; k < i_im.size(); k++) {
        const double a = mt_genrand_res53();
        const double b = mt_genrand_res53();

        o_imNoisy[k] += p_sigma * (float) (sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
    }

    if (p_verbose) {
        cout << "done." << endl;
    }
}

/**
 * @brief Compute PSNR and RMSE between i_im1 and i_im2
 *
 * @param i_im1 : pointer to an allocated array of pixels;
 * @param i_im2 : pointer to an allocated array of pixels;
 * @param o_psnr  : will contain the PSNR;
 * @param o_rmse  : will contain the RMSE;
 * @param p_imageName: name of the image;
 * @param p_verbose: if true, print values of PSNR and RMSE.
 *
 * @return EXIT_FAILURE if both images haven't the same size.
 **/
int computePsnr(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   float &o_psnr
,   float &o_rmse
,   const char* p_imageName
,   const bool p_verbose
){
    if (i_im1.size() != i_im2.size()) {
        cout << "Can't compute PSNR & RMSE: images have different sizes: " << endl;
        cout << "i_im1 : " << i_im1.size() << endl;
        cout << "i_im2 : " << i_im2.size() << endl;
        return EXIT_FAILURE;
    }

    float sum = 0.f;
    for (unsigned k = 0; k < i_im1.size(); k++)
        sum += (i_im1[k] - i_im2[k]) * (i_im1[k] - i_im2[k]);

    o_rmse = sqrtf(sum / (float) i_im1.size());
    o_psnr = 20.f * log10f(255.f / o_rmse);

    if (p_verbose) {
        cout << p_imageName << endl;
        cout << "PSNR = " << o_psnr << endl;
        cout << "RMSE = " << o_rmse << endl;
    }

    return EXIT_SUCCESS;
}

/**
 * @brief Compute a difference image between i_im1 and i_im2.
 *
 * @param i_im1: reference image;
 * @param i_im2: image to compare;
 * @param o_imDiff: will contain the difference;
 * @param p_min, p_max : range of data (usually [0, 255]);
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    std::vector<float> const& i_im1
,   std::vector<float> const& i_im2
,   std::vector<float> &o_imDiff
,   const float p_min
,   const float p_max
,   const bool p_verbose
){
    if (i_im1.size() != i_im2.size()) {
        cout << "Can't compute difference, i_im1 and i_im2 don't have the same size" << endl;
        cout << "i_im1 : " << i_im1.size() << endl;
        cout << "i_im2 : " << i_im2.size() << endl;
        return EXIT_FAILURE;
    }

    if (p_verbose) {
        cout << "Compute difference..." << endl;
    }

    const unsigned size = i_im1.size();
    if (o_imDiff.size() != size) {
        o_imDiff.resize(size);
    }

    //! Auto-adjust by computing the RMSE
    float psnr, rmse;
    computePsnr(i_im1, i_im2, psnr, rmse, "", false);

    const float s = (rmse > 5.f ? 4.f * rmse : 5.f);

    for (unsigned k = 0; k < size; k++) {
        float value =  (i_im1[k] - i_im2[k] + s) * p_max / (2.f * s);
        o_imDiff[k] = (value < p_min ? p_min : (value > p_max ? p_max : value));
    }

    if (p_verbose) {
        cout << "done." << endl;
    }

    return EXIT_SUCCESS;
}
