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

struct Parameters{
	// Picture informations
	unsigned w;				// width
	unsigned h;				// height

	// Dictionnary informations
	unsigned k;				// number of elements
	unsigned m;				// squared of patch size (m = sPatch^2)

	// Patches informations
	unsigned nPatch;		// number
	unsigned sPatch;		// size
	unsigned nRowPatches;	// number of patches in a row
	unsigned nColPatches;	// number of patches in a column
	
	// LARS algorithm
	float reg;				// regularization coefficient TODO: but isn t it chosen wrt the current patch/picture ?

	// Update algorithm
	unsigned update_iteration;	// TODO: Mairal seems to have set it to 1 (and 5 in batch case <== ?)
};

#endif // PARAMETERS_INCLUDED
    