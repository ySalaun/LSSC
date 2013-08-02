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
* @file LibLSSC.cpp
* @brief Library for LSSC algorithm
*
* @author Yohann Salaun <yohann.salaun@polytechnique.org> & Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
**/
#include "LibLSSC.h"

using namespace std;

#ifdef __linux__
    #include <limits>
    #include <math.h>
#else
    #define INFINITY	std::numeric_limits<float>::max()
#endif
#define EPSILON		1e-15

//! train the algorithm with the L1 norm
void trainL1(
    Matrix &io_dict,
    const std::vector<float> &i_imNoisy,
    const unsigned int p_nPatch,
    const Parameters &params){

    //! For convenience
    const unsigned int sP = params.sPatch;
    const unsigned int w  = params.w;

    //! matrices used for the dictionnary update
    Matrix A(params.k, params.k);
    Matrix B(params.m, params.k);

    //! sparse coefficients of the dictionnary
    vector<float> alpha(params.k, 0.f);

    //! generation of p_nPatch iid patches
    display("- generation of iid patches", params, false);
    vector<unsigned int> iidPatches;
    getRandList(p_nPatch, params.nPatch, iidPatches);
    display("...done", params);

    for (unsigned int n = 0; n < p_nPatch; n++) {
        vector<float> patch(params.m, 0.f);
        const unsigned int i = iidPatches[n] / params.nRowPatches * params.sPatch;
        const unsigned int j = iidPatches[n] % params.nRowPatches * params.sPatch;
        const float* iI      = &i_imNoisy[i * params.w + j];
        float* iP            = &patch[0];

        //! initialize patch vector
        for (unsigned int p = 0; p < sP; p++) {
            for (unsigned int q = 0; q < sP; q++) {
                iP[p * sP + q] = iI[p * w + q];
            }
        }

        //! LARS
        display("- LARS", params);
        vector<float> alpha;
        computeLars(io_dict, patch, params, alpha);

        add_xyT(A, alpha, alpha);
        add_xyT(B, patch, alpha);

        //! Update
        display("- dictionary update", params, false);
        updateDictionary(io_dict, A, B, params);
        display("...done", params);
    }
}

//! generation of a list of random number
void getRandList(
    const unsigned int p_nb,
    const unsigned int p_maxRange,
    std::vector<unsigned int> &o_randList){

    //! list with all possibilities
    vector<unsigned int> list(p_nb);
    for(unsigned int n = 0; n < p_maxRange; n++) {
        list[n] = n;
    }

    //! list of random numbers
    o_randList.resize(p_nb);

    //! seed the random function
    srand(time(NULL));

    //! find the p_nb random numbers
    for(unsigned int n = 0; n < p_nb; n++){
        const unsigned int index = rand() % list.size();
        o_randList[n] = list[index];
        list.erase(list.begin() + index);
    }
}

//! lars algorithm that minimizes ||alpha||_1 s.t. ||i_noisy - dict*alpha||_2 < lambda
void computeLars(
    const Matrix &p_dict,
    const std::vector<float> &p_patch,
    const Parameters &p_params,
    std::vector<float> &o_alpha){

    //! For convenience
    const unsigned int nb = p_params.k;
    const float reg       = p_params.reg;

    //! Initialization
    display("-- initialization", p_params);
    o_alpha.assign(nb, 0);

    //! active set indexes
    vector<int> activeIndexes(0);

    //! noisy picture norm
    float norm = dotProduct(p_patch, p_patch);
    //! if the norm is lower than the regularization parameter, stop the algorithm
    if (norm > reg) {
        return;
    }

    //! most correlated element
    vector<float> correlation = product_Ax(p_dict, p_patch, true);
    float cMax = abs(correlation[0]);
    unsigned int currentIndex = 0;
    for (unsigned int n = 1; n < correlation.size(); n++) {
        if (fabs(correlation[n]) > cMax) {
            cMax = abs(correlation[n]);
            currentIndex = n;
        }
    }

    //! begin by adding a new atom in the code
    bool newAtom = true;

    //! matrices parameters
    Matrix G    (nb, nb);
    Matrix Ga   (nb, nb);
    Matrix Gs   (nb, nb);
    Matrix invGs(nb, nb);
    product_AB(p_dict, p_dict, G, true);

    //! add the new most correlated element at each iteration
    for (unsigned int i = 0; i < nb; i++) {

        //! UPDATE
        if (newAtom) {
            display("-- update", p_params);
            activeIndexes.push_back(currentIndex);
            Ga.copyRow(G , currentIndex, i);
            Gs.copyCol(Ga, currentIndex, i);
            Gs.symmetrizeUpperPart(i);
            updateGram(invGs, Gs, i);
        }

        //! VARIABLES UPDATE
        //! compute sign vector
        //! sgn = sgn(c)
        vector<float> sgn(i + 1, 1.f);
        for (unsigned int j = 0; j <= i; j++) {
            if (correlation[j] < 0) {
                sgn[j] = -1.f;
            }
        }
        display("sgn: ", p_params, false);
        display(sgn);

        //! compute direction vector
        //! Ua = invGs * sgn
        vector<float> Ua = product_Ax(invGs, sgn, false, i + 1);
        display("Ua: ", p_params, false);
        display(Ua);
        display("Gs: ", p_params, false);
        display(Gs, i + 1);
        display("invGs: ", p_params, false);
        display(invGs, i + 1);

        //! STEP
        display("-- compute step", p_params);
        float step = INFINITY;
        //float value; // c'est une variable locale, ça sert à rien de la déclarer ici
        //! current maximum of correlation
        cMax = correlation[0];
        vector<float> tGaUa = product_Ax(Ga, Ua, true, i + 1);

        for (unsigned int j = 0; j <= i; j++) {
            tGaUa[activeIndexes[j]] = INFINITY;
        }
        for (unsigned int j = 0; j < nb; j++) {
            if (tGaUa[j] == INFINITY) {
                continue;
            }

            if (tGaUa[j] > -1) {
                const float value = (cMax + correlation[j]) / (1 + tGaUa[j]);
                if (value < step) {
                    step = value;
                    currentIndex = j;
                }
            }

            if (tGaUa[j] < 1) {
                const float value = (cMax - correlation[j]) / (1 - tGaUa[j]);
                if (value < step) {
                    step = value;
                    currentIndex = j;
                }
            }
        }

        //! DOWNDATE STEP
        display("-- compute downdate step", p_params);
        //! if this step is reached, downdate
        float stepDowndate = INFINITY;
        //! Marc: WARNING : indexDowndate not used after that !!!
        //int indexDowndate;
        for (unsigned int j = 0; j <= i; j++) {
            const float ratio = -o_alpha[i] / Ua[i];
            if (ratio < stepDowndate && ratio >= 0) {
                stepDowndate = ratio;
                //indexDowndate = j;
            }
        }
        if (stepDowndate < 0) {
            cout << "step downdate sign issue" << endl;
            //! Marc: il manque pas un truc là ?
        }

        //! STOPPING STEP
        display("-- compute stopping step", p_params);
        cout << sgn.size() << "/" << Ua.size() << "/" << correlation.size() << "/" << i + 1 << endl;
        const float a = dotProduct(sgn, Ua);
        const float b = dotProduct(correlation, Ua, i + 1);
        const float c = norm - reg;
        const float delta = b * b - a * c;
        const float stepStop = min((delta > 0)? (b - sqrtf(delta)) / a : INFINITY, cMax);

        //! TAKE THE STEP
        display("-- take the step", p_params);
        float finalStep = min(step, min(stepDowndate, stepStop));
        for (unsigned int j = 0; j <= i; j++) {
            o_alpha[activeIndexes[j]] += finalStep * Ua[j];
        }
        for (unsigned int j = 0; j < nb; j++) {
            correlation[j] -= finalStep * tGaUa[j];
        }
        norm += (a * finalStep - 2 * b) * finalStep;

        //! STOPPING CONDITIONS
        if (fabs(finalStep) < EPSILON ||
            finalStep == stepStop		||
            norm < EPSILON				||
            norm - reg < EPSILON	){
            //! Marc: pourquoi ce break est commenté ? à quoi sert ce test du coup ?
            //break;
        }
        cout << "final step " << finalStep << endl;
        //! DOWNDATE
        if (false && finalStep == stepDowndate) {
            display("-- downdate", p_params);
            // TODO
            i -= 2;
            newAtom = false;
        }
        else {
            newAtom = true;
        }
    }
}

//! Update the inverse of the Gram matrix
void updateGram(
    Matrix &io_invGs,
    Matrix const& i_Gs,
    const unsigned int p_iter){

    //! case when the matrix is 1x1
    if (p_iter == 0) {
        cout << "beware:" << i_Gs.matrix[0] << endl;
        io_invGs.matrix[0] = 1.f / i_Gs.matrix[0];
    }
    //! General case
    else {
        //! iRow = Gs[i-th row]
        vector<float> iRow = ((Matrix&) i_Gs).row(p_iter, p_iter);
        cout << "iRow: ";
        display(iRow);

        //! u = Gs^-1 Gs[i-th row]
        vector<float> u = product_Ax(io_invGs, iRow, false, p_iter);
        cout << "u: ";
        display(u);

        //! sigma = 1/(Gs_ii - u. Gs[i-th row])
        const float sigma = 1.f / (((Matrix&)i_Gs)(p_iter, p_iter) - dotProduct(u, iRow));
        cout << "sigma: " << sigma << endl;

        //!Gs^-1_iter,iter = sigma;
        io_invGs(p_iter, p_iter) = sigma;

        //! Gs^-1[i-th row] = - sigma u
        for (unsigned int j = 0; j < p_iter; j++) {
            io_invGs(p_iter, j) = -sigma * u[j];
        }

        //! Gs^-1 = Gs^-1 + sigma u u^t
        add_xyT(io_invGs, u, u, p_iter);
    }
}

//! dictionary update algorithm
// Marc : travailler sur u^t tout du long, afin d'optimiser les boucles !!
void updateDictionary(
    Matrix &io_D,
    const Matrix &i_A,
    const Matrix &i_B,
    const Parameters &params){

    //! D = Matrix(D.nRow, D.nCol);
    for(unsigned int iter = 0; iter < params.updateIteration; iter++) {

        //! u = DA
        Matrix u(params.m, params.k);
        product_AB(io_D, i_A, u);

        //! u = B - u
        add(i_B, u, u, true);

        //! uj = uj/ajj for each column j of u
        for (unsigned int j = 0; j < params.k; j++) {
            for (unsigned int i = 0; i < params.m; i++) {
                // TODO Beware, this condition should really be there ?
                if (i_A.matrix[j * params.k + j] != 0){
                    u.matrix[i * params.k + j] /= i_A.matrix[j * params.k + j];
                }
            }
        }
        /*for (unsigned int j = 0; j < params.k; j++) {
            if (i_A(j, j) != 0) {
                const float ajjInv = 1.f / i_A(j ,j);
                for (unsigned int i = 0; i < params.m; i++) {
                    u(i, j) *= ajjInv;
                }
            }
        }*/

        //! u = u + D
        add(u, io_D, u);

        //! uj = uj / max(||u||_2, 1) for each column j
        for (unsigned int j = 0; j < params.k; j++) {
            float norm2 = 0.f;

            for (unsigned int i = 0; i < params.m; i++) {
                norm2 += u.matrix[i * params.k + j] * u.matrix[i * params.k + j];
            }

            if (norm2 > 1.f) {
                const float norm2Inv = 1.f / sqrtf(norm2);
                for (unsigned int i = 0; i < params.m; i++) {
                    u.matrix[i * params.k + j] *= norm2Inv;
                }
            }
        }

        //! D = u
        io_D = u;
    }
}













