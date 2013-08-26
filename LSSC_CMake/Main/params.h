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

  // Verbose option
  bool verbose;

  // Constants
  float INFINITY;
  float EPSILON;

  Parameters(int imH, int imW){
    h               = imH;
    w               = imW;
    
    sPatch          = 9;
    m               = sPatch * sPatch;
    k               = 512;
    nPatch          = w*h/m;
    nRowPatches     = w/sPatch;
    nColPatches     = h/sPatch;
    
    reg             = 1e7; // TODO: compute the real value
    
    updateIteration = 1; // TODO Mairal used this parameter as default
    
    verbose         = true;

    INFINITY        = std::numeric_limits<float>::max();
    EPSILON         = 1e-15;  // TODO: change ? float / double ? 1e-15 or another number ?
  }
};

#endif // PARAMETERS_INCLUDED
