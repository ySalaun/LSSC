/*
* Copyright (c) 2013, Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
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
* @file params.h
* @brief Library for parameters used in LSSC algorithm
*
* @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
**/

#ifndef PARAMETERS_INCLUDED
#define PARAMETERS_INCLUDED

#include <limits>

struct Parameters{
  // Picture informations
  unsigned int w;           // width
  unsigned int h;           // height
  unsigned int wh;          // number of pixels
  unsigned int chnls;       // number of channels

  double sigma;             // noise s.d.

  // Dictionnary informations
  unsigned int k;           // number of elements
  unsigned int m;           // squared of patch size (m = sPatch^2)

  // Patches informations
  unsigned int nPatch;      // number
  unsigned int sPatch;      // size
  unsigned int nRowPatches; // number of patches in a row
  unsigned int nColPatches; // number of patches in a column

  // LARS algorithm
  float reg;                // regularization coefficient TODO: but isn t it chosen wrt the current patch/picture ?

  // Update algorithm
  unsigned updateIteration; // TODO: Mairal seems to have set it to 1 (and 5 in batch case <== ?)

  // ORMP algorithm
  double C; // from KSV IPOL algo
  double epsORMP;

  // Verbose option
  bool verbose;

  // Constants
  float infinity;
  float epsilon;

  // Debug
  bool debug;

  Parameters(
    unsigned int imH,
    unsigned int imW,
    unsigned int imC,
    double s){
    h               = imH;
    w               = imW;
    wh              = w*h;
    chnls           = imC;

    sigma           = s;

    // from KSVD IPOL
    sPatch          =  (sigma  <= 20 ? 5 :
                       (sigma  <= 60 ? 7 : 9)); 
    m               = sPatch * sPatch;
    k               = 512;
    nRowPatches     = h-sPatch+1;
    nColPatches     = w-sPatch+1;
    nPatch          = nRowPatches*nColPatches;

    reg             = 0; // TODO: compute the real value

    updateIteration = 1; // TODO Mairal used this parameter as default


    C               = (chnls == 1 ? (sPatch == 5 ? 1.2017729876383829 : (sPatch == 7 ? 1.1456550151825420 : 1.1139195378939404))					
                                                                      : (sPatch == 5 ? 1.1182997771678573 : (sPatch == 7 ? 1.0849724948297015 : 1.0662877194412401)));
    epsORMP         = ((double) (chnls * m)) * C * C * sigma / 255.l * sigma / 255.l;

    verbose         = true;

    infinity        = std::numeric_limits<float>::max();
    epsilon         = 1e-15f;  // TODO: change ? float / double ? 1e-15 or another number ?

    debug           = false;
  }
};

#endif // PARAMETERS_INCLUDED
