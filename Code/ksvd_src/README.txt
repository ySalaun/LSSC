% K-SVD image denoising.

# ABOUT

* Author    : Marc Lebrun <marc.lebrun@cmla.ens-cachan.fr>
* Copyright : (C) 2011 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the K-SVD image denoising.

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the make program.

- Library. 
This code requires the libpng library.

- Image format. 
Only the PNG format is supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Compile the source code (on Unix/Linux/Mac OS). 
There are two ways to compile the code. 
(1) RECOMMENDED, with Open Multi-Processing multithread parallelization 
(http://openmp.org/). Roughly speaking, it accelerates the program using the 
multiple processors in the computer. Run
make OMP=1

OR
(2) If the complier does not support OpenMp, run 
make

3. Run K-SVD image denoising.
./ksvd
(1) you can moreover want to compute the bias (algorithm applied to the original
image). To do this, use doBias = true.
There are multiple ways to run the code:
(2) if you want to use the speed-up trick, use doSpeedUp = true (recommanded)
(3) if you want to obtain the best PSNR result, use doBestPSNR = true
 
Example, run
./ksvd cinput.png 10 ImNoisy.png ImDenoised.png ImDiff.png ImBias.png ImDiffBias.png 0 1 0



# ABOUT THIS FILE

Copyright 2011 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
