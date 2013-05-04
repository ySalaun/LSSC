template <typename T>
inline void Image<T>::denoise6b(const T sigma, Matrix<T>& D, Image<T>& noisy, const Image<T>& orig, Trainer<T>& trainer,
								const int n, const T J,
								const T* param_C, const T thrsb,
								const int window,
								const int num_threads) {
	const T thrs = abs<T>(thrsb);
	Timer time;
	time.start();
	PRINT_F(param_C[_m])
	PRINT_F(J)
	const int batch_size=static_cast<int>(floor(256*((num_threads-1)/2+1)));
	int JJ;
	printf("Matching\n");
	list_groups listsl;
	this->grouping4(listsl,n,window,thrs,MIN(4,1));
	vector_groups lists(listsl.begin(),listsl.end());
	time.printElapsed();
	if (J != 0) {
		if (J > 0)  {
			JJ = static_cast<int>(ceil(J*static_cast<T>(_n)/batch_size))+1;
			printf("Number of patches processed : %d\n",JJ*batch_size);
		} else {
			JJ = static_cast<int>(floor(J));
		}
		printf("Group learning\n");
		trainer.setLambda(sigma*sigma);
		noisy.setDataVariables(n);
		list_groups listsl2;
		this->grouping5(listsl2,n,window,thrs,MIN(4,1));
		vector_groups listsv(listsl2.begin(),listsl2.end());
		trainer.train(noisy,listsv,JJ,L2ERROR,true,param_C,1,false);
		trainer.getD(D);
	}
	this->copy(noisy);
	int NUM_THREADS=1;
#ifdef _OPENMP
	NUM_THREADS=num_threads;
	omp_set_num_threads(NUM_THREADS);
	omp_set_nested(0);
	omp_set_dynamic(0);
#endif
	Image<T> I2(_w,_h,_V);
	I2.setZeros();
	Image<T> count(_w,_h,_V);
	count.setZeros();
	const int numPatchesX =_h-n+1;
	const int numPatchesY =_w-n+1;
	const int numPatches = numPatchesX*numPatchesY;
	const int K=D.n();
	const int L = MIN(K,_m);

	Matrix<T> G;
	D.XtX(G);
	int Mcount=0;
	int i;
	int num_groups=lists.size();
	PRINT_I(num_groups)
#pragma omp parallel for private(i)
	for (i = 0; i<num_groups; ++i) {
#ifdef _OPENMP
		int numT=omp_get_thread_num();
#else
		int numT=0;
#endif
		group& list = lists[i];
		const int M = list.size();
		Matrix<T> X(_m,M);
		Vector<T> Xj;
		/// extract data
		int j=0;
		for (group::iterator it = list.begin(); it != list.end(); ++it) {
			X.refCol(j++,Xj);
			const int py = (*it)/numPatchesX;
			const int px = (*it)-py*numPatchesX;
			this->getPatch(Xj,px,py,n);
		}
		Vector<T> mean(_V);
		X.whiten(mean,false);
		const int ind2 = MIN(9999,M*_m);
		const T eps = param_C[ind2]*param_C[ind2]*sigma*sigma*_m*M;

		if (thrsb > 0) {
			/// perform sparse decomposition
			Vector<int> ind;
			Matrix<T> vM;

			coreSOMP(X,D,G,vM,ind,L,eps);

			/// reconstruct
			Vector<T> d;
			X.setZeros();
			for (int j = 0; j<vM.m(); ++j) {
				D.refCol(ind[j],d);
				cblas_ger<T>(CblasColMajor,_m,M,T(1.0),
							d.rawX(),1,vM.rawX()+j,vM.m(),X.rawX(),_m);
			}
		} else {
			Matrix<T> RtD;
			X.mult(D,RtD,true,false);
			Matrix<T> alphat(M,K);
			coreGroupISTConstrained2(X,G,RtD,alphat,
										eps,100,T(0.001));
			D.mult(alphat,X,false,true);
		}

		X.unwhiten(mean,false);
		/// update image and count
		j = 0;
		for (group::iterator it = list.begin(); it != list.end(); ++it) {
			X.refCol(j++,Xj);
			const int py = (*it)/numPatchesX;
			const int px = (*it)-py*numPatchesX;
			I2.addPatch(Xj,px,py,n);
			count.incrPatch(px,py,n);
		}
		list.clear();
#pragma omp critical
		Mcount += M;
	}
	PRINT_I(numPatches)
	PRINT_I(Mcount)
	I2.scalByInv(count);
	I2.clip();
	printf("\t\t\tPSNR : %g\n",I2.psnr(orig));
	this->copy(I2);
	time.stop();
	time.printElapsed();
};