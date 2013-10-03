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

void testFunction(Matrix& dict, Parameters& params)
{
    bool testFailed = false;
    float epsilon = 0.005;

    //! UPDATE GRAM
    Matrix G    (params.k, params.k);
    G.productAtB(dict, dict);

    //! Test with random values if G = D^t * D
    srand((unsigned int)(time(time_t(NULL))));

    for(unsigned int test = 0; test < 25; test++)
    {
        const unsigned int i = rand() % params.k;
        const unsigned int j = rand() % params.k;

        float sum = 0.f;
        for(unsigned k=0; k<params.m; k++)
        {
            sum += dict(k,i)*dict(k,j);
        }

        if(fabs(sum - G(i,j)) > epsilon)
        {
            cout << "TEST ERROR: the productAtB function does not work for Gram matrix" << endl;
            cout << "G(" << i << "," << j << ") = " << G(i,j) << " != " << sum << endl;
            testFailed = true;
        }
    }

    if(testFailed)
    {
        cout << "---------------FAILED---------------" << endl;
    }

    // Test of update Gram matrix
    epsilon = 0.000005;

    Matrix invG    (params.k, params.k);
    Matrix Id1     (params.k, params.k);
    Matrix Id2     (params.k, params.k);
    Matrix Gs      (params.k, params.k);
    Matrix Ga      (params.k, params.k);

    const unsigned int iterMax = 70;
    for(unsigned int iter = 0; iter<iterMax; iter++)
    {
        updateGram(invG, G, iter);

        Id1.productAB(G, invG, iter+1);
        Id2.productAB(invG, G, iter+1);

        for(unsigned int i=0; i<iter+1; i++)
        {
            for(unsigned int j=0; j<iter+1; j++)
            {
                const float ref = (i==j)? 1.f : 0.f;

                if(fabs(ref - Id1(i,j)) > epsilon)
                {
                    cout << "TEST ERROR: the updateGram function does not work" << endl;
                    cout << "at " << iter << "-th iteration, G*invG(" << i << "," << j << ") = " << Id1(i,j) << " != " << ref << endl;
                    testFailed = true;
                }

                if(fabs(ref - Id2(i,j)) > epsilon)
                {
                    cout << "TEST ERROR: the updateGram function does not work" << endl;
                    cout << "at " << iter << "-th iteration, invG*G(" << i << "," << j << ") = " << Id2(i,j) << " != " << ref << endl;
                    testFailed = true;
                }

            }
        }

        for(unsigned int i=0; i<iter+1; i++)
        {
            for(unsigned int j=0; j<i; j++)
            {

                if(fabs(G(i,j) - G(j,i)) > epsilon)
                {
                    cout << "TEST ERROR: Gram matrix is not symmetric" << endl;
                    cout << "at " << iter << "-th iteration, G(" << i << "," << j << ") = " << G(i,j) << " != " << " G(" << j << "," << i << ") = " << G(j,i) << endl;
                    testFailed = true;
                }

                if(fabs(invG(i,j) - invG(j,i)) > epsilon)
                {
                    cout << "TEST ERROR: the updateGram function does not work, inverse of Gram matrix is not symmetric" << endl;
                    cout << "at " << iter << "-th iteration, invG(" << i << "," << j << ") = " << invG(i,j) << " != " << " invG(" << j << "," << i << ") = " << invG(j,i) << endl;
                    testFailed = true;
                }

            }
        }

        if(testFailed)
        {
            /*cout << "-------------------------------" << endl;
            display(G, iter+1);
            cout << "-------------------------------" << endl;
            display(invG, iter+1);
            cout << "-------------------------------" << endl;
            display(Id1, iter+1);
            cout << "-------------------------------" << endl;
            display(Id2, iter+1);
            cout << "-------------------------------" << endl;*/
            cout << "---------------FAILED---------------" << endl;
        }
    }
    cout << "------------------" << endl;
    vector<int> A(iterMax +1, 0);
    for(unsigned int iter = iterMax; iter>0; iter--)
    {
        downdateGram(invG, G, Ga, A, iter-1, 0);
        Id1.productAB(G, invG, iter-1);
        Id2.productAB(invG, G, iter-1);

        for(unsigned int i=0; i<iter-1; i++)
        {
            for(unsigned int j=0; j<iter-1; j++)
            {
                const float ref = (i==j)? 1.f : 0.f;

                if(fabs(ref - Id1(i,j)) > epsilon)
                {
                    cout << "TEST ERROR: the downdateGram function does not work" << endl;
                    cout << "at " << iter << "-th iteration, G*invG(" << i << "," << j << ") = " << Id1(i,j) << " != " << ref << endl;
                    testFailed = true;
                }

                if(fabs(ref - Id2(i,j)) > epsilon)
                {
                    cout << "TEST ERROR: the downdateGram function does not work" << endl;
                    cout << "at " << iter << "-th iteration, invG*G(" << i << "," << j << ") = " << Id2(i,j) << " != " << ref << endl;
                    testFailed = true;
                }

            }
        }

        for(unsigned int i=0; i<iter-1; i++)
        {
            for(unsigned int j=0; j<i; j++)
            {

                if(fabs(G(i,j) - G(j,i)) > epsilon)
                {
                    cout << "TEST ERROR: Gram matrix is not symmetric" << endl;
                    cout << "at " << iter << "-th iteration, G(" << i << "," << j << ") = " << G(i,j) << " != " << " G(" << j << "," << i << ") = " << G(j,i) << endl;
                    testFailed = true;
                }

                if(fabs(invG(i,j) - invG(j,i)) > epsilon)
                {
                    cout << "TEST ERROR: the downdateGram function does not work, inverse of Gram matrix is not symmetric" << endl;
                    cout << "at " << iter << "-th iteration, invG(" << i << "," << j << ") = " << invG(i,j) << " != " << " invG(" << j << "," << i << ") = " << invG(j,i) << endl;
                    testFailed = true;
                }

            }
        }

        if(testFailed)
        {
            /*cout << "-------------------------------" << endl;
            display(G, iter+1);
            cout << "-------------------------------" << endl;
            display(invG, iter+1);
            cout << "-------------------------------" << endl;
            display(Id1, iter+1);
            cout << "-------------------------------" << endl;
            display(Id2, iter+1);
            cout << "-------------------------------" << endl;*/
            cout << "---------------FAILED---------------" << endl;
        }
    }
}

void denoiseL0(
  const vector<float>& i_imNoisy,
  vector<double>& o_imDenoised,
  const Matrix& i_dict,
  const Parameters& params)
{
    vector<unsigned int> median(i_imNoisy.size(), 0);

    //! size of patches
    const unsigned int sP = params.sPatch;

    //! write the picture into patch form
    vector<double> i_X(params.nPatch * params.m);

    for(unsigned int i = 0; i < params.nRowPatches; i++)
    {
        for(unsigned int j = 0; j < params.nColPatches; j++)
        {
            //! patch with top left corner im_noisy(i,j)
            const float* iI = &i_imNoisy[i * params.w + j];
            double* iP      = &i_X[(i * params.nColPatches + j) * params.m];

            for (unsigned int p = 0; p < sP; p++)
            {
                for (unsigned int q = 0; q < sP; q++)
                {
                    iP[p * sP + q] = (double) iI[p * params.w + q];
                }
            }
        }
    }

    //! write the dictionary into vector form
    vector<double> i_D(params.m * params.k);

    for(unsigned int i = 0; i < params.m; i++)
    {
        for(unsigned int j = 0; j < params.k; j++)
        {
            i_D[j * params.m + i] = (double) i_dict(i, j);
        }
    }

    //! denoised patch output
    vector<vector<unsigned int> > o_indV;
    vector<vector<double> > o_valV;

    //! compute ORMP
    computeORMP(i_X, i_D, o_indV, o_valV, params.epsORMP, params.nPatch, params.k, params.m);

    //! denoise each patch separately
    for(unsigned int i = 0; i < params.nRowPatches; i++)
    {
        for(unsigned int j = 0; j < params.nColPatches; j++)
        {
            const unsigned int numPatch = i * params.nColPatches + j;
            vector<double> coeff        = o_valV[numPatch];
            vector<unsigned int> code   = o_indV[numPatch];

            //! patch with top left corner im_noisy(i,j)
            double* iP      = &i_X[numPatch * params.m];

            for (unsigned int p = 0; p < sP; p++)
            {
                for (unsigned int q = 0; q < sP; q++)
                {
                    double denoisedPixel          = 0;
                    const unsigned int indexPatch = p * sP + q;
                    for(unsigned int r = 0; r < code.size(); r++)
                    {
                        denoisedPixel += coeff[r] * i_D[code[r] * params.m + indexPatch];
                    }
                    iP[indexPatch] = denoisedPixel;
                    o_imDenoised[(i + p) * params.w + (j + q)] += denoisedPixel;
                    median[(i + p) * params.w + (j + q)] += 1;
                }
            }
        }
    }

    //! merge dneoised patch information
    for(unsigned int i = 0; i < params.nRowPatches; i++)
    {
        for(unsigned int j = 0; j < params.nColPatches; j++)
        {
            o_imDenoised[i * params.w + j] /= min(float(median[i * params.w + j]), 1.f);
        }
    }
}

int main(
    int argc,
    char **argv
)
{
    //! Check if there is the right call for the algorithm
    if (argc < 10)
    {
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
    if (loadImage(argv[1], im, imSize, verbose) != EXIT_SUCCESS)
    {
        return EXIT_FAILURE;
    }

    //! Add noise
    addNoise(im, imNoisy, sigma, verbose);

    //! Parameters
    // TODO remark, for now the pictures are considered to be B&W
    // TODO the parameters have to be smarter
    Parameters params(imSize.height, imSize.width, 1, sigma);

    //! read initial dictionary
    Matrix dict(params.m, params.k);
    ifstream txtDict(argv[9], ios::in);

    for(unsigned int i = 0; i < params.m; i++)
    {
        for(unsigned int j = 0; j < params.k; j++)
        {
            txtDict >> dict(i, j);
        }
    }

    //! TEST PART
    bool test = false;
    if(test)
    {
        display("----------------------------------------------", params);
        display("-------------------TEST-----------------------", params);
        testFunction(dict, params);
    }

    //! LSSC
    display("----------------------------------------------", params);

    display("PART 1 - LEARNING PART WITH LARS", params);
    //! number of random patches for dictionary update
    unsigned nRandomPatches = 50;//unsigned(floor(.2 * params.nPatch));

    //! dictionary training with LARS/update
    //trainL1(dict, imNoisy, nRandomPatches, params);

    display("PART 2 - ORMP from KSVD", params);

    //! roughly denoised picture
    vector<double> denoi1(imNoisy.size(), (double)0);

    //! pre-denoise picture with ORMP and previously learned dictionary
    denoiseL0(imNoisy, denoi1, dict, params);

    display("PART 3 - CLUSTERING", params);
    // TODO

    display("PART 4 - LEARNING WITH SIMULTANEOUS LARS", params);
    // TODO

    display("PART 5 - SIMULTANEOUS ORMP", params);
    // TODO

    // TODO Marc: ne jamais faire de push_back lorsqu'on peut éviter, c'est ultra lent.
    // TODO Marc: c'est normal de ne copier que le premier canal ?
    // TODO Yohann: concrètement ça sert à rien ces lignes et faudra les changer pour en faire ce qui nous intéresse ;)
    for (unsigned c = 0; c < imSize.nChannels; c++)
    {
        for (unsigned i = 0; i < imSize.height; i++)
        {
            for (unsigned j = 0; j < imSize.width; j++)
            {
                imFinal.push_back((float) denoi1[0 * imSize.wh + i * imSize.width + j]);
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
    if (doBias)
    {
        computePsnr(im, imBias, psnrBias, rmseBias, "imBiasFinal", verbose);
    }

    //! save noisy, denoised and differences images
    if (verbose)
    {
        cout << "Save images...";
    }
    if (saveImage(argv[3], imNoisy, imSize, 0.f, 255.f) != EXIT_SUCCESS)
    {
        return EXIT_FAILURE;
    }

    if (saveImage(argv[4], imFinal, imSize, 0.f, 255.f) != EXIT_SUCCESS)
    {
        return EXIT_FAILURE;
    }

    if (saveImage(argv[5], imDiff, imSize, 0.f, 255.f) != EXIT_SUCCESS)
    {
        return EXIT_FAILURE;
    }

    if (doBias)
    {
        if (saveImage(argv[6], imBias, imSize, 0.f, 255.f) != EXIT_SUCCESS)
        {
            return EXIT_FAILURE;
        }

        if (saveImage(argv[7], imDiffBias, imSize, 0.f, 255.f) != EXIT_SUCCESS)
        {
            return EXIT_FAILURE;
        }
    }
    if (verbose)
    {
        cout << "done." << endl;
    }

    //! Exit the Main Function
    return EXIT_SUCCESS;
}
