/* ******************
 * Simultaneous OMP 
 * *****************/

template <typename T>
void somp(const Matrix<T>* X, const Matrix<T>& D, SpMatrix<T>* spalpha, 
      const int Ngroups, const int L, const T eps,const int numThreads) {
   somp(X,D,spalpha,Ngroups,L,&eps,false,numThreads);
}

template <typename T>
void somp(const Matrix<T>* XT, const Matrix<T>& D, SpMatrix<T>* spalphaT, 
      const int Ngroups, const int LL, const T* eps, const bool adapt,
      const int numThreads) {
   if (LL <= 0) return;
   const int K = D.n();
   const int L = MIN(D.m(),MIN(LL,K));

   if (!D.isNormalized()) {
      cerr << "Current implementation of OMP does not support non-normalized dictionaries" << endl;
      return;
   }

   /// compute the Gram Matrix G=D'D
   Matrix<T> G;
   D.XtX(G);

   int NUM_THREADS=init_omp(numThreads);

   int i;
#pragma omp parallel for private(i) 
   for (i = 0; i< Ngroups; ++i) {
      const Matrix<T>& X = XT[i];
      const int M = X.n();
      SpMatrix<T>& spalpha = spalphaT[i];
      spalpha.clear();
      Vector<int> rv;
      Matrix<T> vM;
      T thrs = adapt ? eps[i] : M*(*eps);
      coreSOMP(X,D,G,vM,rv,L,thrs);
      spalpha.convert2(vM,rv,K);   
   }
}
