void ist_groupLasso(const Matrix<T>* XT, const Matrix<T>& D,
      Matrix<T>* alphaT, const int Ngroups, 
      const T lambda, const constraint_type mode,
      const int itermax,
      const T tol, const int numThreads) {
   int K=D.n();
   int n = D.m();

   if (!D.isNormalized()) {
      cerr << "Current implementation of block coordinate descent does not support non-normalized dictionaries" << endl;
      return;
   }

   if (mode == L1COEFFS) {
      std::cerr << "Mode not implemented" << std::endl;
      return;
   }


   /// compute the Gram Matrix G=D'D
   Matrix<T> G;
   D.XtX(G);

   int NUM_THREADS=init_omp(numThreads);

   Matrix<T>* RtDT = new Matrix<T>[NUM_THREADS];
   Matrix<T>* alphatT = new Matrix<T>[NUM_THREADS];

   int i;
#pragma omp parallel for private(i) 
   for (i = 0; i< Ngroups; ++i) {
#ifdef _OPENMP
      int numT=omp_get_thread_num();
#else
      int numT=0;
#endif
      const Matrix<T>& X = XT[i];
      int M = X.n();
      Matrix<T>& alphat = alphatT[numT];
      alphaT[i].transpose(alphat);
      Matrix<T>& RtD = RtDT[numT];
      X.mult(D,RtD,true,false);


      Vector<T> col, col2;
      T norm1 = alphat.asum();
      T normX2;

      if (!norm1) {
         Vector<T> DtR_mean(K);
         Vector<T> coeffs_mean(K);
         coeffs_mean.setZeros();
         RtD.meanRow(DtR_mean);
         coeffs_mean.setZeros();
         if (mode == PENALTY) {
            coreIST(G,DtR_mean,coeffs_mean,lambda/T(2.0),itermax,tol);
         } else {
            Vector<T> meanVec(n);
            X.meanCol(meanVec);
            normX2=meanVec.nrm2sq(); 
            coreISTconstrained(G,DtR_mean,coeffs_mean,normX2,
                  lambda,itermax,tol);
            SpVector<T> spalpha(K);
            normX2-=computeError(normX2,G,DtR_mean,coeffs_mean,spalpha);
            normX2=X.normFsq()-M*normX2;
         }
         alphat.fillRow(coeffs_mean);         
      }

      if (M > 1) {
         for (int j = 0; j<K; ++j) {
            alphat.refCol(j,col);
            const T nrm=col.nrm2sq();
            if (nrm) {
               G.refCol(j,col2);
               RtD.rank1Update(col,col2,T(-1.0));
            }
         }

         if (mode == PENALTY) {
            coreGroupIST(G,RtD,alphat,sqr<T>(M)*lambda/T(2.0),itermax,sqr<T>(M)*tol);
         } else  {
            coreGroupISTConstrained(G,RtD,alphat,normX2,M*lambda,itermax,sqr<T>(M)*tol);
         }
      }
      alphat.transpose(alphaT[i]);
   }

   delete[](RtDT);
   delete[](alphatT);
};