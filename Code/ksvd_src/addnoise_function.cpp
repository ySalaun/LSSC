#include "mt19937ar.h"
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

void fiAddNoise(float *u, double *v, double std, long int randinit, unsigned size)
{
	mt_init_genrand((unsigned long int) time (NULL) + (unsigned long int) getpid()  + (unsigned long int) randinit);

    for (unsigned i = 0; i < size; i++)
    {
        double a = mt_genrand_res53();
        double b = mt_genrand_res53();
        double z = (double)(std) * sqrtl(-2.0 * log(a)) * cos(2.0 * M_PI * b);

        v[i] =  (double) u[i] + z;
    }
}

