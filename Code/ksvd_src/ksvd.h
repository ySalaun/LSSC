/**
    ksvd.h
**/

#ifndef KSVD_H_INCLUDED
#define KSVD_H_INCLUDED

#include<vector>

typedef std::vector<std::vector<double> > matD_t;
typedef std::vector<std::vector<unsigned> > matU_t;
typedef std::vector<double> vecD_t;
typedef std::vector<unsigned> vecU_t;
typedef std::vector<double>::iterator iterD_t;
typedef std::vector<unsigned>::iterator iterU_t;

//! Run K-SVD
void ksvd_ipol(const double   sigma       ,
               double *       img_noisy   ,
               float *        img_denoised,
               const unsigned width       ,
               const unsigned height      ,
               const unsigned chnls       ,
               const bool     useT        );

//! Decompose the image into patches of size NxN
void im2patches(const double * img,
                matD_t        &patches,
                const unsigned width,
                const unsigned height,
                const unsigned chnls,
                const unsigned N);

//! Reconstruct images by aggregating patches of size N x N
void patches2im(matD_t        &patches,
                 float *        img,
                 const double * img_ref,
                 const unsigned width,
                 const unsigned height,
                 const unsigned chnls,
                 const double   lambda,
                 const unsigned N);

//! Obtain random permutation for a tabular of size N
void randperm(vecU_t &perm);

 //! Obtain the initial dictionary, which
 //! columns are normalized
void obtain_dict(matD_t        &dictionary,
                 matD_t const&  patches);

//! Process K-SVD denoising
void ksvd_process(const double *img_noisy,
                  float       * img_denoised,
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
				  const bool     doReconstruction);

#endif // KSVD_H_INCLUDED
