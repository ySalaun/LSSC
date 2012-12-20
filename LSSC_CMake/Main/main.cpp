/*
 * Copyright (c) 2012, Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
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
 * @file main.cpp
 * @brief main function.
 *
 * @author Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
 **/

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "libImages/LibImages.h"

using namespace std;

int main(
    int argc
,   char **argv
){
    //! Check if there is the right call for the algorithm
	if (argc < 9) {
		cout << "usage: LSSC image sigma noisy denoised difference \
		bias diff_bias computeBias" << endl;
		return EXIT_FAILURE;
	}

	 //! Variables initialization
	const float sigma = float(atof(argv[2]));
	const bool doBias = bool(atof(argv[8]) != 0);

    //! Declarations
	vector<float> im, imNoisy, imFinal, imDiff;
	vector<float> imBias, imDiffBias;
	ImageSize imSize;
    const bool verbose = false;

    //! Read Image
    if (loadImage(argv[1], im, imSize, verbose) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

	//! Add noise
    addNoise(im, imNoisy, sigma, verbose);

    //! LSSC
    for (unsigned c = 0; c < imSize.nChannels; c++) {
        for (unsigned i = 0; i < imSize.height; i++) {
            for (unsigned j = 0; j < imSize.width; j++) {
                imFinal[c * imSize.wh + i * imSize.width + j] =
                    imNoisy[c * imSize.wh + i * imSize.width + j];
            }
        }
    }

    //! Compute PSNR and RMSE
    float psnr, rmse;
    computePsnr(im, imFinal, psnr, rmse, "imFinal", verbose);

    float psnrBias, rmseBias;
    if (doBias) {
        computePsnr(im, imBias, psnrBias, rmseBias, "imBiasFinal", verbose);
    }

    //! save noisy, denoised and differences images
	if (verbose) {
	    cout << "Save images...";
	}
	if (saveImage(argv[3], imNoisy, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	if (saveImage(argv[4], imFinal, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

    if (saveImage(argv[5], imDiff, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
    }

    if (doBias) {
        if (saveImage(argv[6], imBias, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }

        if (saveImage(argv[7], imDiffBias, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    if (verbose) {
        cout << "done." << endl;
    }

    //! Exit the Main Function
    return EXIT_SUCCESS;
}
