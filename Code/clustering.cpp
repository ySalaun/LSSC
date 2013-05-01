/// perform grouping
template <typename T>
inline void Image<T>::grouping(list_groups& lists,
		const int n, const int window, const T thrs,
		const int num_threads) const {
	PRINT_F(thrs)
	lists.clear();
#ifdef _OPENMP
	omp_set_num_threads(num_threads);
	omp_set_nested(0);
	omp_set_dynamic(0);
#endif
	const int numPatchesX =_h-n+1;
	const int numPatchesY =_w-n+1;
	const int numPatches = numPatchesX*numPatchesY;
	std::vector<bool> used(numPatchesX*numPatchesY,false);
	for (int j =0; j<numPatches; ++j)
		used[j]=false;
	const int sizePatch=n*n*_V;
	Vector<T>* currentT = new Vector<T>[num_threads];
	Vector<T>* tmpT= new Vector<T>[num_threads];
	for (int i = 0; i<num_threads; ++i) {
		currentT[i].resize(sizePatch);
		tmpT[i].resize(sizePatch);
	}
	const T eps = thrs/(sizePatch);
	int i;
	#pragma omp parallel for private(i)
	for (i = 0; i<numPatches; ++i) {
	#ifdef _OPENMP
		int numT=omp_get_thread_num();
	#else
		int numT=0;
	#endif
	Vector<T>& current = currentT[numT];
	Vector<T>&  tmp = tmpT[numT];
	const int jj = i/ numPatchesX;
	const int ii = i-jj*numPatchesX;
	if (used[ii+jj*numPatchesX]) continue;
	group list;
	this->getPatch(current,ii,jj,n);
	// const int maxg = 200;
		const int xm = MAX(0,ii-window);
		const int xM = MIN(_h-n,ii+window);
		const int ym = MAX(0,jj-window);
		const int yM = MIN(_w-n,jj+window);
		int sizegroup=0;
		for (int kk = xm; kk<=xM; ++kk) {
			for (int ll = ym; ll<=yM; ++ll) {
				// if (kk != ii && ll != jj && sizegroup > maxg) continue;
				this->getPatch(tmp,kk,ll,n);
				tmp.sub(current);
				if (tmp.nrm2sq() < eps) {
					list.push_back(kk+ll*numPatchesX);
					used[kk+ll*numPatchesX]=true;
					++sizegroup;
				}
			}
		}
		if (list.size() > 0)
		#pragma omp critical
			 lists.push_back(list);
	}
	delete[](tmpT);
	delete[](currentT);
};