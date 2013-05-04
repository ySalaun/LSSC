template <typename T>
inline void Image<T>::denoise4(const T sigma, Matrix<T>& D,
							const Image<T>& orig, const int n, const T eps,
							const int num_threads) {
	int NUM_THREADS=1;
#ifdef _OPENMP
	NUM_THREADS=num_threads;
	omp_set_num_threads(NUM_THREADS);
	omp_set_nested(0);
	omp_set_dynamic(0);
#endif
	Timer time_global;
	time_global.start();
	Image<T> I2(_w,_h,_V);
	I2.setZeros();
	const int numPatchesX =_h-n+1;
	const int numPatchesY =_w-n+1;
	const int numPatches = numPatchesX*numPatchesY;
	const int sizePatch=n*n*_V;
	const int K=D.n();
	if (sizePatch != D.m()) {
		cerr << "Bad dictionary" << endl;
		return;
	}
	const int L = MIN(K,sizePatch);
	Vector<T>* XT=new Vector<T>[NUM_THREADS];
	Vector<T>* scoresT=new Vector<T>[NUM_THREADS];
	Vector<T>* normT=new Vector<T>[NUM_THREADS];
	Vector<T>* tmpT=new Vector<T>[NUM_THREADS];
	Vector<T>* RdnT=new Vector<T>[NUM_THREADS];
	Vector<T>* RUnT=new Vector<T>[NUM_THREADS];
	Vector<int>* indT=new Vector<int>[NUM_THREADS];
	Matrix<T>* UnT=new Matrix<T>[NUM_THREADS];
	Matrix<T>* UndnT=new Matrix<T>[NUM_THREADS];
	Matrix<T>* UndsT=new Matrix<T>[NUM_THREADS];
	Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
	for (int i = 0; i<NUM_THREADS; ++i) {
		scoresT[i].resize(K);
		normT[i].resize(K);
		tmpT[i].resize(K);
		RdnT[i].resize(K);
		RUnT[i].resize(L);
		indT[i].resize(L);
		UnT[i].resize(L,L);
		UnT[i].setZeros();
		UndnT[i].resize(K,L);
		UndsT[i].resize(L,L);
		GsT[i].resize(K,L);
		XT[i].resize(sizePatch);
	}
	int i;
	Matrix<T> G;
	D.XtX(G);
#pragma omp parallel for private(i)
	for (i = 0; i<numPatches; ++i) {
#ifdef _OPENMP
		int numT=omp_get_thread_num();
#else
		int numT=0;
#endif
		const int py = i / numPatchesX;
		const int px = i-py*numPatchesX;
		Vector<T>& X = XT[numT];
		Vector<T> mean(_V);
		this->getPatch(X,px,py,n);
		X.whiten(mean,false);
		T normX = X.nrm2sq();
		Vector<T>& Rdn=RdnT[numT];
		Vector<int>& ind=indT[numT];
		Vector<T>& RUn=RUnT[numT];
		ind.set(-1);
		D.multTrans(X,Rdn);
		coreORMP(scoresT[numT],normT[numT],tmpT[numT],UnT[numT],UndnT[numT],UndsT[numT],
					GsT[numT],Rdn,G,ind,RUn,eps, normX);
		X.setZeros();
		for (int j = 0; j<L; ++j) {
			Vector<T> d;
			if (ind[j] == -1) break;
			D.refCol(ind[j],d);
			X.add(d,RUn[j]);
		}
		X.unwhiten(mean,false);
		I2.addPatch(X,px,py,n);
	}
	const int maw= MIN(_w-n+1,n);
	const int mah= MIN(_h-n+1,n);
	for (int i = 0; i<_w; ++i) {
		int num1 = MIN(MIN(maw,i+1),_w-i);
		for (int j = 0; j<_h; ++j) {
			int num2 = MIN(MIN(mah,j+1),_h-j);
			int num=num1*num2;
			for (int k = 0; k<_V; ++k) {
				I2(j,i,k) /= static_cast<T>(num);
			}
		}
	}

	delete[](XT);
	delete[](scoresT);
	delete[](normT);
	delete[](tmpT);
	delete[](RdnT);
	delete[](RUnT);
	delete[](indT);
	delete[](UnT);
	delete[](UndnT);
	delete[](UndsT);
	delete[](GsT);
	I2.clip();
	printf("\t\t\tPSNR : %g\n",I2.psnr(orig));
	this->copy(I2);
	time_global.stop();
	time_global.printElapsed();
}
