#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <vector>

//! Compute PSNR and RMSE
void psnr_rmse(const float *  img_1,
               const float *  img_2,
               float *        psnr,
               float *        rmse,
               const unsigned size);

#endif // UTILITIES_H_INCLUDED
