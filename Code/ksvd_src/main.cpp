#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "addnoise_function.h"
#include "io_png.h"
#include "utilities.h"
#include "ksvd.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author MARC LEBRUN  <marc.lebrun@cmla.ens-cachan.fr>
 */


int main(int argc, char **argv)
{
    //! Check if there is the right call for the algorithm
	if (argc < 10)
	{
		cout << "usage: K-SVD image sigma noisy denoised difference bias diff_bias \
                                                            doBias doSpeedUp" << endl;
		return EXIT_FAILURE;
	}

	//! read input image
	cout << "Read input image...";
	size_t height, width, chnls;
	float *img = NULL;
	img = io_png_read_f32(argv[1], &width, &height, &chnls);
	if (!img)
	{
		cout << "error :: " << argv[1] << " not found  or not a correct png image" << endl;
		return EXIT_FAILURE;
	}
	cout << "done." << endl << endl;

	//! test if image is really a color image and exclude the alpha channel
	if (chnls > 2)
	{
	    unsigned k = 0;
	    while (k < width * height \
            && img[k] == img[width * height + k] \
            && img[k] == img[2 * width * height + k])
            k++;
        chnls = (k == width * height ? 1 : 3);
	}

	//! Printing some characterics of the input image
    cout << "image size : " << endl;
    cout << "-width    = " << width << endl;
    cout << "-height   = " << height << endl;
    cout << "-channels = " << chnls << endl << endl;

    //! Variables initialization
	double fSigma   = atof(argv[2]);
	unsigned wh     = (unsigned) width * height;
	unsigned whc    = (unsigned) wh * chnls;

	//! Add noise
	cout << "Add noise [sigma = " << fSigma << "] ...";
	double *img_noisy    = new double[whc];
	float  *img_denoised = new float [whc];
	float  *img_bias     = new float [whc];

	for (unsigned c = 0; c < chnls; c++)
		fiAddNoise(&img[c * wh], &img_noisy[c * wh], fSigma, c, wh);
    cout << "done." << endl;

	//! Denoising
	bool useAcceleration = atof(argv[9]);
	cout << "Applying K-SVD to the noisy image :" << endl;
    ksvd_ipol((double) fSigma / 255.0l, img_noisy, img_denoised, width, height, chnls,
              useAcceleration);
    cout << endl;

    if (atof(argv[8]))
    {
        cout << "Applying K-SVD to the original image :" << endl;
        double *img_noisy_bias = new double[whc];
        for (unsigned k = 0; k < whc; k++)
            img_noisy_bias[k] = (double) img[k];
        ksvd_ipol((double) fSigma / 255.0l, img_noisy_bias, img_bias, width, height,
                  chnls, useAcceleration);
        delete[] img_noisy_bias;
        cout << endl;
    }

    //! Compute RMSE and PSNR
    float rmse, rmse_bias, psnr, psnr_bias;
    psnr_rmse(img, img_denoised, &psnr, &rmse, whc);
    cout << endl;
    cout << "For noisy image :" << endl;
    cout << "PSNR: " << psnr << endl;
    cout << "RMSE: " << rmse << endl << endl;
    if (atof(argv[8]))
    {
        psnr_rmse(img, img_bias, &psnr_bias, &rmse_bias, whc);
        cout << "For original image :" << endl;
        cout << "PSNR: " << psnr_bias << endl;
        cout << "RMSE: " << rmse_bias << endl << endl;
    }

	//! writing measures
    char path[13] = "measures.txt";
    ofstream file(path, ios::out);
    if(file)
    {
        file << "************" << endl;
        file << "-sigma = " << fSigma << endl;
        file << "-PSNR  = " << psnr << endl;
        file << "-RMSE  = " << rmse << endl;
        if (atof(argv[8]))
        {
            file << "-PSNR_bias  = " << psnr_bias << endl;
            file << "-RMSE_bias  = " << rmse_bias << endl << endl;
        }
        file.close();
    }
    else
        return EXIT_FAILURE;

	//! Compute Difference
	cout << "Compute difference...";
	fSigma *= 4.0f;
	float *img_diff      = new float[whc];
	float *img_diff_bias = new float[whc];
	float fValue, fValue_bias;

    #pragma omp parallel for
        for (unsigned k = 0; k < whc; k++)
        {
            fValue =  (img[k] - img_denoised[k] + fSigma) * 255.0f / (2.0f * fSigma);
            fValue_bias =  (img[k] - img_bias[k] + fSigma) * 255.0f / (2.0f * fSigma);
            img_diff[k] = (fValue < 0.0f ? 0.0f : (fValue > 255.0f ? 255.0f : fValue));
            img_diff_bias[k] = (fValue_bias < 0.0f ? 0.0f : \
                                (fValue_bias > 255.0f ? 255.0f : fValue_bias));
        }
	cout << "done." << endl << endl;

	//! save noisy, denoised and differences images
	cout << "Save images...";
	float * img_noisy_f = new float[whc];
	for (unsigned k = 0; k < whc; k++)
        img_noisy_f[k] = (float) (img_noisy[k] * 255.0l);

	if (io_png_write_f32(argv[3], img_noisy_f, width, height, chnls) != 0)
		cout << "... failed to save png image " << argv[3] << endl;

	if (io_png_write_f32(argv[4], img_denoised, width, height, chnls) != 0)
		cout << "... failed to save png image " << argv[4] << endl;

    if (io_png_write_f32(argv[5], img_diff, width, height, chnls) != 0)
		cout << "... failed to save png image " << argv[5] << endl;
    if (atof(argv[8]))
    {
        if (io_png_write_f32(argv[6], img_bias, width, height, chnls) != 0)
            cout << "... failed to save png image " << argv[6] << endl;

        if (io_png_write_f32(argv[7], img_diff_bias, width, height, chnls) != 0)
            cout << "... failed to save png image " << argv[7] << endl;
    }
    cout << "done." << endl;

	//! Free Memory
	delete[] img_denoised;
	delete[] img_noisy;
	delete[] img_noisy_f;
	delete[] img_diff;
	delete[] img_bias;
	delete[] img_diff_bias;

	return EXIT_SUCCESS;
}



