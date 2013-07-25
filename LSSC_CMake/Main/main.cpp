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
#include <ctime>

#include "libImages/LibImages.h"
#include "libLSSC/LibLSSC.h"

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
    const bool verbose = true;

    //! Read Image
    if (loadImage(argv[1], im, imSize, verbose) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
    }

	//! Add noise
    addNoise(im, imNoisy, sigma, verbose);
	
	//! Parameters
	// TODO remark, for now the pictures are considered to be B&W
	// TODO the parameters have to be smarter
	Parameters params;
	params.h = imSize.height;
	params.w = imSize.width;
	params.sPatch = 16;
	params.m = params.sPatch * params.sPatch;
	params.k = 512;
	params.nPatch = imSize.wh/params.m;
	params.nRowPatches = params.w/params.sPatch;
	params.nColPatches = params.h/params.sPatch;
	params.reg = 1e7; // TODO: compute the real value
	params.update_iteration = 1; // TODO see into Mairal's code

  //! LSSC
  Matrix dict(params.m, params.k);
	// PART 1
  // LEARNING PART WITH LARS
  unsigned nRandomPatches = 10;//unsigned(floor(.2 * params.nPatch));
	trainL1(dict, imNoisy, nRandomPatches, params);
  
  // PART 2
  // ORMP from KSVD

  // TODO

  // PART 3
  // CLUSTERING

  // TODO

  // PART 4
  // LEARNING WITH SIMULTANEOUS LARS

  // TODO

  // PART 5
  // SIMULTANEOUS ORMP

  // TODO

	for (unsigned c = 0; c < imSize.nChannels; c++) {
        for (unsigned i = 0; i < imSize.height; i++) {
            for (unsigned j = 0; j < imSize.width; j++) {
				imFinal.push_back(imNoisy[0 * imSize.wh + i * imSize.width + j]);
				imDiff.push_back(imNoisy[0 * imSize.wh + i * imSize.width + j]);
				imBias.push_back(imNoisy[0 * imSize.wh + i * imSize.width + j]);
				imDiffBias.push_back(imNoisy[0 * imSize.wh + i * imSize.width + j]);
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
