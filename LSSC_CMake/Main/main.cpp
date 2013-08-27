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
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctime>

#ifdef __linux__
#include "../libImages/LibImages.h"
#include "../libLSSC/LibLSSC.h"
#include "../libORMP/LibOrmp.h"
#else
#include "libImages/LibImages.h"
#include "libLSSC/LibLSSC.h"
#include "libORMP/LibORMP.h"
#endif


using namespace std;

int main(
int argc
,   char **argv
){
  //! Check if there is the right call for the algorithm
  if (argc < 10) {
    cout << "usage: LSSC image sigma noisy denoised difference \
            bias diff_bias computeBias dict.txt" << endl;
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
  Parameters params(imSize.height, imSize.width);

  //! read initial dictionary
  Matrix dict(params.m, params.k);
  ifstream txtDict(argv[9], ios::in);

  for(unsigned int i = 0; i < params.m; i++){
    for(unsigned int j = 0; j < params.k; j++){
      txtDict >> dict(i, j);
    }
  }

  //! LSSC
  display("----------------------------------------------", params);
  display("PART 1 - LEARNING PART WITH LARS", params);
  unsigned nRandomPatches = 50;//unsigned(floor(.2 * params.nPatch));
  trainL1(dict, imNoisy, nRandomPatches, params);

  display("PART 2 - ORMP from KSVD", params);
  // TODO

  display("PART 3 - CLUSTERING", params);
  // TODO

  display("PART 4 - LEARNING WITH SIMULTANEOUS LARS", params);
  // TODO

  display("PART 5 - SIMULTANEOUS ORMP", params);
  // TODO

  // TODO Marc: ne jamais faire de push_back lorsqu'on peut éviter, c'est ultra lent.
  // TODO Marc: c'est normal de ne copier que le premier canal ?
  // TODO Yohann: concrètement ça sert à rien ces lignes et faudra les changer pour en faire ce qui nous intéresse ;)
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
