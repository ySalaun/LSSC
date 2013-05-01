/* *********************
 * Implementation of OMP
 * *********************/

/// Forward Selection (or Orthogonal matching pursuit) 
/// Address the problem of:
/// \forall i, \min_{\alpha_i} ||X_i-D\alpha_i||_2^2 
///                        s.t. ||\alphai||_0 <= L or
/// \forall i, \min_{\alpha_i} ||\alpha_i||_0 
///                        s.t. ||\X_i-D\alpha_i||_2^2 <= epsilon
/// This function is 
///   * efficient (Cholesky-based)
///   * parallel
///   * optimized for a big number of signals (precompute the Gramm matrix

template <typename T>
void omp(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha, 
      const int* pL, const T* peps, const T* pLambda, 
      const bool vecL, const bool vecEps,
      const bool vecLambda, const int numThreads, Matrix<T>* path) {
   int L;
   if (!vecL) {
      L=*pL;
   } else {
      Vector<int> vL(const_cast<int*>(pL),X.n());
      L=vL.maxval();
   }
   spalpha.clear();
   if (L <= 0) return;
   const int M = X.n();
   const int K = D.n();
   L = MIN(X.m(),MIN(L,K));
   Matrix<T> vM(L,M);
   Matrix<int> rM(L,M);

   ProdMatrix<T> G(D, K < 25000 && M > 10);

   int NUM_THREADS=init_omp(numThreads);

   Vector<T>* scoresT=new Vector<T>[NUM_THREADS];
   Vector<T>* normT=new Vector<T>[NUM_THREADS];
   Vector<T>* tmpT=new Vector<T>[NUM_THREADS];
   Vector<T>* RdnT=new Vector<T>[NUM_THREADS];
   Matrix<T>* UnT=new Matrix<T>[NUM_THREADS];
   Matrix<T>* UndnT=new Matrix<T>[NUM_THREADS];
   Matrix<T>* UndsT=new Matrix<T>[NUM_THREADS];
   Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
   for (int i = 0; i<NUM_THREADS; ++i) {
      scoresT[i].resize(K);
      normT[i].resize(K);
      tmpT[i].resize(K);
      RdnT[i].resize(K);
      UnT[i].resize(L,L);
      UnT[i].setZeros();
      UndnT[i].resize(K,L);
      UndsT[i].resize(L,L);
      GsT[i].resize(K,L);
   }

   int i;
#pragma omp parallel for private(i) 
   for (i = 0; i< M; ++i) {
#ifdef _OPENMP
      int numT=omp_get_thread_num();
#else
      int numT=0;
#endif
      Vector<T> Xi;
      X.refCol(i,Xi);
      T normX = Xi.nrm2sq();

      Vector<int> ind;
      rM.refCol(i,ind);
      ind.set(-1);

      Vector<T> RUn;
      vM.refCol(i,RUn);

      Vector<T>& Rdn=RdnT[numT];
      D.multTrans(Xi,Rdn);
      coreORMP(scoresT[numT],normT[numT],tmpT[numT],UnT[numT],UndnT[numT],UndsT[numT],
            GsT[numT],Rdn,G,ind,RUn, normX, vecEps ? peps+i : peps,
            vecL ? pL+i : pL, vecLambda ? pLambda+i : pLambda, 
            path && i==0 ? path->rawX() : NULL);
   }

   delete[](scoresT);
   delete[](normT);
   delete[](tmpT);
   delete[](RdnT);
   delete[](UnT);
   delete[](UndnT);
   delete[](UndsT);
   delete[](GsT);

   /// convert the sparse matrix into a proper format
   spalpha.convert(vM,rM,K);
};