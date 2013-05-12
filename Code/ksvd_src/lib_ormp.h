#ifndef LIBORMP_H_INCLUDED
#define LIBORMP_H_INCLUDED

#include <vector>

typedef std::vector<std::vector<double> > matD_t;
typedef std::vector<std::vector<unsigned> > matU_t;
typedef std::vector<double> vecD_t;
typedef std::vector<unsigned> vecU_t;
typedef std::vector<double>::iterator iterD_t;
typedef std::vector<unsigned>::iterator iterU_t;

//! ORMP process
void ormp_process(matD_t      &X,
                  matD_t      &D,
                  matU_t      &ind_v,
                  matD_t      &val_v,
                  unsigned     L,
                  const double eps);

//! sub function for ORMP process
void coreORMP(matD_t      &D,
              matD_t      &D_D,
              vecD_t      &scores,
              vecD_t      &norm,
              matD_t      &A,
              matD_t      &D_ELj,
              matD_t      &D_DLj,
              vecD_t      &x_T,
              vecU_t      &ind,
              vecD_t      &coord,
              const double eps2,
              double       normr);

//! Find best ind
unsigned ind_fmax(vecD_t &V);


#endif // LIBORMP_H_INCLUDED
