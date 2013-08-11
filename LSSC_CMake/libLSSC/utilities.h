/*
* Copyright (c) 2013, Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun.ik@gmail.com>
* All rights reserved.
*
* This program is free software: you can use, modify and/or
* redistribute it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later
* version. You should have received a copy of this license along
* this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <stdlib.h>
#include <iostream>
#include <vector>

#ifdef __linux__
#include "../Main/params.h"
#include "../libMatrix/ClassMatrix.h"
#else
#include "Main/params.h"
#include "libMatrix/ClassMatrix.h"
#endif

using namespace std;

void display(const char* msg, const Parameters &params, bool endline = true);
void display(const vector<float>& vec);
void display(const Matrix & mat, const int iMax = -1);

#endif // UTILITIES_H_INCLUDED
