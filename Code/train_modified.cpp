template <typename T>
void coreLARS(Vector<T>& Rdnv, Vector<T>& Xdnv, Vector<T>& Av,
      Vector<T>& uv, Vector<T>& sigv, Vector<T>& avv, Vector<T>& RUnv,
      Matrix<T>& Unm, Matrix<T>& Undsm, Matrix<T>& Gsm,
      Matrix<T>& Gsam, Matrix<T>& workm, Matrix<T>& Rm, 
      const AbstractMatrix<T>& Gm,T& normX, 
      Vector<int>& indv,Vector<T>& coeffsv,const T constraint,
      const bool ols,const bool pos, constraint_type mode,
      T* path, int length_path) {
   if (mode == L2ERROR && normX < constraint) return;

   const int LL = Gsm.n();
   const int K = Gsm.m();
   const int L = MIN(LL,K);
   if (length_path <= 1) length_path=4*L;
   // permit unsafe fast low level access
   T* const Rdn = Rdnv.rawX();
   T* const Xdn = Xdnv.rawX();
   T* const A = Av.rawX();
   T* const u = uv.rawX();
   T* const sig = sigv.rawX();
   T* const av = avv.rawX();
   T* const RUn = RUnv.rawX();
   T* const Un = Unm.rawX();
   T* const Unds = Undsm.rawX();
   T* const Gs = Gsm.rawX();
   T* const Gsa = Gsam.rawX();
   T* const work = workm.rawX();
   //T* const G = Gm.rawX();
   T* const R = Rm.rawX();
   int* ind = indv.rawX();
   T* coeffs = coeffsv.rawX();

   coeffsv.setZeros();
   indv.set(-1);

   if (ols) Xdnv.copy(Rdnv);
   int currentInd= pos ? Rdnv.max() : Rdnv.fmax();
   bool newAtom=true;
   T Cmax;
   int iter=1;
   T thrs = 0.0;

   int* const ind_orig = ind;
   T* const coeffs_orig = coeffs;

   int j;
   for (j = 0; j<L; ++j) {
      if (newAtom) {
         ind[j]=currentInd;

         Cmax = abs<T>(Rdn[currentInd]);
         sig[j] = SIGN(Rdn[currentInd]);

         for (int k = 0; k<=j; ++k) Un[j*L+k]=0.0;
         Un[j*L+j]=1.0;
         Gm.extract_rawCol(currentInd,Gs+K*j);
         for (int k = 0; k<j; ++k) Gs[K*j+ind[k]] *= sig[k];
         if (sig[j] < 0) {
            Rdn[currentInd]=-Rdn[currentInd];
            if (ols) Xdn[currentInd]=-Xdn[currentInd];
            cblas_scal<T>(K,sig[j],Gs+K*j,1);
            cblas_scal<T>(j+1,sig[j],Gs+currentInd,K);
         }
         cblas_copy<T>(j+1,Gs+currentInd,K,Gsa+j*L,1);
         for (int k = 0; k<j; ++k) Gsa[k*L+j]=Gsa[j*L+k];

         // <d_j,d_i>
         cblas_copy<T>(j,Gsa+j*L,1,Unds+j,L);
         // <U_j final,d_i>
         cblas_trmv<T>(CblasColMajor,CblasUpper,CblasTrans,CblasNonUnit,
               j+1,Un,L,Unds+j,L);
         // norm2
         T norm2=Gsa[j*L+j];
         for (int k = 0; k<j; ++k) norm2 -= Unds[k*L+j]*Unds[k*L+j];
         if (norm2 < 1e-15) {
            ind[j]=-1;
      //      cerr << "bad exit" << endl;
            break;
         }
      
      //   int iter2 = norm2 < 0.5 ? 2 : 1;
      //   for(int k = 0; k<iter2; ++k) {
      //      for (int l = 0; l<j; ++l) {
      //         T scal=-cblas_dot<T>(j+1-l,Un+j*L+l,1,Unds+l*L+l,1);
      //         cblas_axpy<T>(l+1,scal,Un+l*L,1,Un+j*L,1);
      //      }
      //   }
         Un[j*L+j]=-T(1.0);
         cblas_copy<T>(j,Unds+j,L,Un+j*L,1);
         cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,j,Un,L,Un+j*L,1);

         /// Un is the orthogonalized vectors in the D basis
         T invNorm=1.0/sqrt(norm2);
         cblas_scal<T>(j+1,-invNorm,Un+j*L,1);
         Unds[j*L+j]=cblas_dot<T>(j+1,Un+j*L,1,Gsa+j*L,1);
      }

      for (int k = 0; k<=j; ++k) u[k]=T(1.0);
      cblas_trmv<T>(CblasColMajor,CblasUpper,CblasTrans,CblasNonUnit,
            j+1,Un,L,u,1);

      T a = T(1.0)/cblas_nrm2<T>(j+1,u,1);

      cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
            j+1,Un,L,u,1);
      cblas_scal<T>(j+1,a,u,1);

      cblas_gemv<T>(CblasColMajor,CblasNoTrans,K,j+1,T(1.0),Gs,K,u,1,T(0.0),A,1);

      T potentNorm=0.0;
      if (!ols) {
         for (int k = 0; k<=j; ++k)  potentNorm += Rdn[ind[k]]*u[k];
      }

      if (pos) {
         for (int k = 0; k<K; ++k) {
            T diff = a-A[k];
            work[k]= diff <= 0 ? INFINITY : (Cmax-Rdn[k])/diff;
         }
         for (int k = 0; k<=j; ++k) {
            work[ind[k]]=INFINITY; 
         }
         for (int k = 0; k<K; ++k) 
            if (work[k] <=0) work[k]=INFINITY;
         currentInd =cblas_iamin<T>(K,work,1);
      } else {
         memset(work,0,2*K*sizeof(T));
         for (int k = 0; k<=j; ++k) {
            const int index=2*ind[k];
            work[index]=INFINITY; 
            work[index+1]=INFINITY; 
         }
         for (int k = 0; k<K; ++k) {
            const int index=2*k;
            if (!work[index]) {
               const T diff1=a-A[k];
               work[index]= diff1 <= 0 ? INFINITY : (Cmax-Rdn[k])/diff1;
               const T diff2=a+A[k];
               work[index+1]=diff2 <= 0 ? INFINITY : (Cmax+Rdn[k])/diff2;
            }
         }
         currentInd =cblas_iamin<T>(2*K,work,1);
      }
      T gamma=work[currentInd];
      T gammaMin=0;
      int minBasis=0;

      //if (j == L-1) gamma=potentNorm;

      if (mode == PENALTY) {
         gamma=MIN(gamma,(Cmax-constraint)/a);
      }

//      if (j > 0) {
         vDiv<T>(j+1,coeffs,u,work);
         cblas_scal<T>(j+1,-T(1.0),work,1);
         /// voir pour petites valeurs
         for (int k=0; k<=j; ++k) 
            if (coeffs[k]==0 || work[k] <=0) work[k]=INFINITY;
         minBasis=cblas_iamin<T>(j+1,work,1);
         gammaMin=work[minBasis];
         if (gammaMin < gamma) gamma=gammaMin;
 //     }

      if (mode == L1COEFFS) {
         T Tu = 0.0;
         for (int k = 0; k<=j; ++k) Tu += u[k];

         if (Tu > EPSILON) 
            gamma= MIN(gamma,(constraint-thrs)/Tu);
         thrs+=gamma*Tu;
      }

      // compute the norm of the residdual

      if (ols == 0) {
         const T t = gamma*gamma - 2*gamma*potentNorm;
         if (t > 0 || isnan(t) || isinf(t)) {
      //      cerr << "bad bad exit" << endl;
     //       cerr << t << endl;
            ind[j]=-1;
            break;
         }
         normX += t;
      } else {
         // plan the last orthogonal projection
         if (newAtom) {
            RUn[j]=0.0;
            for (int k = 0; k<=j; ++k) RUn[j] += Xdn[ind[k]]*
               Un[j*L+k];
            normX -= RUn[j]*RUn[j];
         }
      }

      // Update the coefficients
      cblas_axpy<T>(j+1,gamma,u,1,coeffs,1);

      if (pos) {
         for (int k = 0; k<j+1; ++k)
            if (coeffs[k] < 0) coeffs[k]=0;
      }

      cblas_axpy<T>(K,-gamma,A,1,Rdn,1);
      if (!pos) currentInd/= 2;
      if (path) {
         for (int k = 0; k<=j; ++k) 
            path[iter*K+ind[k]]=coeffs[k]*sig[k];
      }

      if (gamma == gammaMin) {
         downDateLasso<T>(j,minBasis,normX,ols,pos,Rdnv,ind,coeffs,sigv,
               avv,Xdnv, RUnv, Unm, Gsm, Gsam,Undsm,Rm);
         newAtom=false;
         Cmax=abs<T>(Rdn[ind[0]]);
         --j;
      } else {
         newAtom=true;
      }
      ++iter;

      if (mode == PENALTY) {
         thrs=abs<T>(Rdn[ind[0]]);
      }

      if ((j == L-1) || 
            (mode == PENALTY && (thrs - constraint < 1e-15)) ||
            (mode == L1COEFFS && (thrs - constraint > -1e-15)) || 
            (newAtom && mode == L2ERROR && (normX - constraint < 1e-15)) ||
            (normX < 1e-15) ||
            (iter >= length_path)) {
     //       cerr << "exit" << endl;
     //       PRINT_F(thrs)
     //       PRINT_F(constraint)
     //       PRINT_F(normX)
         break;
      }

   }
   if (ols) {
      cblas_copy<T>(j+1,RUn,1,coeffs,1);
      cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
            j+1,Un,L,coeffs,1);
   }
   vMul<T>(j+1,coeffs,sig,coeffs);
};


/// Auxiliary function for lasso 
/// solve min_{ alpha } | | alpha | | _1 s . t . | | x-Dalpha | | _2^2 <= lambda
/// @brief
/// @param DtR = Dtx
/// @param G = Gram Matrix = DtD
/// @param Gs, Ga, invGs, u, coeffs, ind, work, normX seems to be references that are filled in the programm below
/// @param mode ==> 2 = L2ERROR ==> solve min_{ alpha } | | alpha | | _1 s . t . | | x-Dalpha | | _2^2 <= lambda
/// @param pos ==> positivity constraint
template <typename T>
void coreLARS2(Vector<T>& DtR, const AbstractMatrix<T>& G,
				Matrix<T>& Gs,
				Matrix<T>& Ga,
				Matrix<T>& invGs,
				Vector<T>& u,
				Vector<T>& coeffs,
				Vector<int>& ind,
				Matrix<T>& work,
				T& normX,
				//const constraint_type mode, ==> not needed cause only L2ERROR mode
				const T constraint,
				const bool pos,
				T* path, int length_path) {
	const int LL = Gs.n();
	const int K = G.n();
	const int L = MIN(LL,K);
	if (length_path <= 1) length_path=4*L;

	coeffs.setZeros();
	ind.set(-1);

	// rawX ==> modifiable ref of the data
	T* const pr_Gs = Gs.rawX();
	T* const pr_invGs = invGs.rawX();
	T* const pr_Ga = Ga.rawX();
	T* const pr_work = work.rawX();
	T* const pr_u = u.rawX();
	T* const pr_DtR = DtR.rawX();
	T* const pr_coeffs = coeffs.rawX();
	int* const pr_ind = ind.rawX();

	// Find the most correlated element ==> c_j in LARS article
	// fmax ==> max in magnitude, seems to be max(abs())
	int currentInd = pos ? DtR.max() : DtR.fmax();
	
	// why issue if normX < constraint ?
	if (normX < constraint) return;
	
	// begin by adding a new atom
	bool newAtom = true;
	
	// iteration parameters
	int i;
	int iter=0;
	// loop
	// add one atom at each iteration except when gone too far and need to come back
	// stop when the criterion is no more checked
	for (i = 0; i<L; ++i) {
		++iter;
		// new atom, update Ga, Gs, invGs
		if (newAtom) {
			pr_ind[i] = currentInd;
			// j-th column of Gram Matrix where j is the current index
			// Ga = DatDa where Da is the dictionnary with only selected index
			G.extract_rawCol(pr_ind[i],pr_Ga+i*K);
			// what is Gs ?
			// seems to be triangular sup version of Ga
			for (int j = 0; j<=i; ++j){
				pr_Gs[i*LL+j]=pr_Ga[i*K+pr_ind[j]];
			}
			// Update inverse of Gs
			if (i == 0) {
				pr_invGs[0]=T(1.0)/pr_Gs[0];
			} else {
				// not seen in details...
				cblas_symv<T>(CblasColMajor,CblasUpper,i,T(1.0),
								pr_invGs,LL,pr_Gs+i*LL,1,T(0.0),pr_u,1);
				const T schur = T(1.0)/(pr_Gs[i*LL+i]-cblas_dot<T>(i,pr_u,1,pr_Gs+i*LL,1));
				pr_invGs[i*LL+i]=schur;
				cblas_copy<T>(i,pr_u,1,pr_invGs+i*LL,1);
				cblas_scal<T>(i,-schur,pr_invGs+i*LL,1);
				cblas_syr<T>(CblasColMajor,CblasUpper,i,schur,pr_u,1,
							pr_invGs,LL);
			}
		}
		// Compute the path direction 
		// pr work[j] = sgn(c_j)
		for (int j = 0; j<=i; ++j){
			pr_work[j]= pr_DtR[pr_ind[j]] > 0 ? T(1.0) : T(-1.0);
		}
		// what does it do ???
		cblas_symv<T>(CblasColMajor,CblasUpper,i+1,T(1.0),pr_invGs,LL,
					pr_work,1,T(0.0),pr_u,1);
		
		// Compute the step on the path
		T step_max = INFINITY;
		int first_zero = -1;
		for (int j = 0; j<=i; ++j) {
			T ratio = -pr_coeffs[j]/pr_u[j];
			if (ratio > 0 && ratio <= step_max) {
				step_max=ratio;
				first_zero=j;
			}
		}
		
		T current_correlation = abs<T>(pr_DtR[pr_ind[0]]);
		cblas_gemv<T>(CblasColMajor,CblasNoTrans,K,i+1,T(1.0),pr_Ga,
						K,pr_u,1,T(0.0),pr_work+2*K,1);
		cblas_copy<T>(K,pr_work+2*K,1,pr_work+K,1);
		cblas_copy<T>(K,pr_work+2*K,1,pr_work,1);

		for (int j = 0; j<=i; ++j) {
			pr_work[pr_ind[j]]=INFINITY;
			pr_work[pr_ind[j]+K]=INFINITY;
		}
		for (int j = 0; j<K; ++j) {
			pr_work[j] = ((pr_work[j] < INFINITY) && (pr_work[j] > T(-1.0))) ? (pr_DtR[j]+current_correlation)/(T(1.0)+pr_work[j]) : INFINITY;
		}
		
		for (int j = 0; j<K; ++j) {
			pr_work[j+K] = ((pr_work[j+K] < INFINITY) && (pr_work[j+K] < T(1.0))) ? (current_correlation-pr_DtR[j])/(T(1.0)-pr_work[j+K]) : INFINITY;
		}
		
		if (pos) {
			for (int j = 0; j<K; ++j) {
				pr_work[j]=INFINITY;
			}
		}
		
		int index = cblas_iamin<T>(2*K,pr_work,1);
		T step = pr_work[index];

		// Choose next element
		currentInd = index % K;

		// compute the coefficients of the polynome representing normX^2
		T coeff1 = 0;
		for (int j = 0; j<=i; ++j){
			coeff1 += pr_DtR[pr_ind[j]] > 0 ? pr_u[j] : -pr_u[j];
		}
		T coeff2 = 0;
		for (int j = 0; j<=i; ++j){
			coeff2 += pr_DtR[pr_ind[j]]*pr_u[j];
		}
		T coeff3 = normX-constraint;
		T step_max2;
		
		/// L2ERROR
		const T delta = coeff2*coeff2-coeff1*coeff3;
		step_max2 = delta < 0 ? INFINITY : (coeff2-sqrt(delta))/coeff1;
		step_max2 = MIN(current_correlation,step_max2);
		
		step = MIN(MIN(step,step_max2),step_max);
		if (step == INFINITY) break; // stop the path
		// Update coefficients
		cblas_axpy<T>(i+1,step,pr_u,1,pr_coeffs,1);

		if (pos) {
			for (int j = 0; j<i+1; ++j)
				if (pr_coeffs[j] < 0) pr_coeffs[j]=0;
		}
		
		// Update correlations
		cblas_axpy<T>(K,-step,pr_work+2*K,1,pr_DtR,1);
		
		// Update normX
		normX += coeff1*step*step-2*coeff2*step;
		
		if (path) {
			for (int k = 0; k<=i; ++k){
				path[iter*K+ind[k]]=pr_coeffs[k];
			}
		}
		
		// Choose next action
		// case 1: when need to go back on the path cause of sign issues <== check
		// case 2: next index, new atom
		if (step == step_max) {
		/// Downdate, remove first_zero
		/// Downdate Ga, Gs, invGs, ind, coeffs
			for (int j = first_zero; j<i; ++j) {
				cblas_copy<T>(K,pr_Ga+(j+1)*K,1,pr_Ga+j*K,1);
				pr_ind[j]=pr_ind[j+1];
				pr_coeffs[j]=pr_coeffs[j+1];
			}
			pr_ind[i]=-1;
			pr_coeffs[i]=0;
			for (int j = first_zero; j<i; ++j) {
				cblas_copy<T>(first_zero,pr_Gs+(j+1)*LL,1,pr_Gs+j*LL,1);
				cblas_copy<T>(i-first_zero,pr_Gs+(j+1)*LL+first_zero+1,1,
								pr_Gs+j*LL+first_zero,1);
			}
			const T schur = pr_invGs[first_zero*LL+first_zero];
			cblas_copy<T>(first_zero,pr_invGs+first_zero*LL,1,pr_u,1);
			cblas_copy<T>(i-first_zero,pr_invGs+(first_zero+1)*LL+first_zero,LL,
							pr_u+first_zero,1);
			for (int j = first_zero; j<i; ++j) {
				cblas_copy<T>(first_zero,pr_invGs+(j+1)*LL,1,pr_invGs+j*LL,1);
				cblas_copy<T>(i-first_zero,pr_invGs+(j+1)*LL+first_zero+1,1,
								pr_invGs+j*LL+first_zero,1);
			}
			cblas_syr<T>(CblasColMajor,CblasUpper,i,T(-1.0)/schur,
							pr_u,1,pr_invGs,LL);
			newAtom=false;
			// why i-2 ? i-1 ??
			i=i-2;
		} else {
			newAtom=true;
		}
		if ((iter >= length_path-1) || abs(step) < 1e-15 ||
				step == step_max2 || (normX < 1e-15) ||
				(i == (L-1)) ||
				(mode == L2ERROR && normX - constraint < 1e-15)) {
			break;
		}
	}
}

template <typename T>
void Trainer<T>::train(const Data<T>& X, const ParamDictLearn<T>& param) {
	T rho = param.rho;
	T t0 = param.t0;
	int sparseD;
	if(param.modeD == L1L2){
		sparseD = 2;
	} else if(param.modeD == L1L2MU){
		sparseD = 7;
	} else{
		sparseD = 6;
	}
	int NUM_THREADS = init_omp(_NUM_THREADS);
	
/* VERBOSE
	if (param.verbose) {
		cout << "num param iterD: " << param.iter_updateD << endl;
		if (param.batch) {
			cout << "Batch Mode" << endl;
		} else if (param.stochastic) {
			cout << "Stochastic Gradient. rho : " << rho << ", t0 : " << t0 << endl;
		} else {
			if (param.modeParam == AUTO) {
				cout << "Online Dictionary Learning with no parameter " << endl;
			} else if (param.modeParam == PARAM1) {
				cout << "Online Dictionary Learning with parameters: " << t0 << " rho: " << rho << endl;
			} else {
				cout << "Online Dictionary Learning with exponential decay t0: " << t0 << " rho: " << rho << endl;
			}
		}
		if (param.posD){
			cout << "Positivity constraints on D activated" << endl;
		}
		if (param.posAlpha){
			cout << "Positivity constraints on alpha activated" << endl;
		 }
		if (param.modeD != L2) cout << "Sparse dictionaries, mode: " << param.modeD << ", gamma1: " << param.gamma1 << ", gamma2: " << param.gamma2 << endl;
			cout << "mode Alpha " << param.mode << endl;
		if (param.clean) cout << "Cleaning activated " << endl;
		if (param.log && param.logName) {
			cout << "log activated " << endl;
			cerr << param.logName << endl;
		}
		if (param.mode == PENALTY && param.lambda==0 && param.lambda2 > 0 && !param.posAlpha)
			cout << "L2 solver is used" << endl;
		if (_itercount > 0)
			cout << "Retraining from iteration " << _itercount << endl;
		flush(cout);
	}
*/

	const int M = X.n();
	const int K = _k;
	const int n = X.m();
	const int L = (param.mode == SPARSITY)? static_cast<int>(param.lambda)
				: (param.mode == PENALTY && param.lambda == 0 && param.lambda2 > 0 && !param.posAlpha)? K 
				: MIN(n,K);
	const int batchsize = param.batch ? M : MIN(_batchsize,M);
	
/* VERBOSE
	if (param.verbose) {
		cout << "batch size: " << batchsize << endl;
		cout << "L: " << L << endl;
		cout << "lambda: " << param.lambda << endl;
		cout << "mode: " << param.mode << endl;
		flush(cout);
	}
*/

	srandom(0);
	Vector<T> col(n);
	
	// detect if initial dictionnary and initializes it
	if (_D.m() != n || _D.n() != K){
		_D.resize(n,K);
		for (int i = 0; i<K; ++i) {
			const int ind = random() % M;
			Vector<T> d;
			_D.refCol(i,d);
			X.getData(col,ind);
			d.copy(col);
		}
	}

/* VERBOSE
   if (param.verbose) {
      cout << "*****Online Dictionary Learning*****" << endl;
      flush(cout);
   }
*/

	Vector<T> tmp(n);
	if (param.modeD != L2) {
		for (int i = 0; i<K; ++i) {
			Vector<T> d;
			_D.refCol(i,d);
			tmp.copy(d);
			tmp.sparseProject(d, T(1.0), sparseD,
				param.gamma1, param.gamma2, T(2.0), param.posD);
		}
	} else {
		if (param.posD){
			_D.thrsPos();
		}
		_D.normalize();
	}

	int count = 0;
	int countPrev = 0;
	T scalt0 =  abs<T>(t0);
	if (_itercount == 0) {
		_A.resize(K,K);
		_A.setZeros();
		_B.resize(n,K);
		_B.setZeros();
		if (!param.batch) {
			_A.setDiag(scalt0);
			_B.copy(_D);
			_B.scal(scalt0);
		}
	}

/* MATRICES INIT ?
	//Matrix<T> G(K,K);
	Matrix<T> Borig(n,K);
	Matrix<T> Aorig(K,K);
	Matrix<T> Bodd(n,K);
	Matrix<T> Aodd(K,K);
	Matrix<T> Beven(n,K);
	Matrix<T> Aeven(K,K);
	SpVector<T>* spcoeffT=new SpVector<T>[_NUM_THREADS];
	Vector<T>* DtRT=new Vector<T>[_NUM_THREADS];
	Vector<T>* XT=new Vector<T>[_NUM_THREADS];
	Matrix<T>* BT=new Matrix<T>[_NUM_THREADS];
	Matrix<T>* AT=new Matrix<T>[_NUM_THREADS];
	Matrix<T>* GsT=new Matrix<T>[_NUM_THREADS];
	Matrix<T>* GaT=new Matrix<T>[_NUM_THREADS];
	Matrix<T>* invGsT=new Matrix<T>[_NUM_THREADS];
	Matrix<T>* workT=new Matrix<T>[_NUM_THREADS];
	Vector<T>* uT=new Vector<T>[_NUM_THREADS];
	for (int i = 0; i<_NUM_THREADS; ++i) {
		spcoeffT[i].resize(K);
		DtRT[i].resize(K);
		XT[i].resize(n);
		BT[i].resize(n,K);
		BT[i].setZeros();
		AT[i].resize(K,K);
		AT[i].setZeros();
		GsT[i].resize(L,L);
		GsT[i].setZeros();
		invGsT[i].resize(L,L);
		invGsT[i].setZeros();
		GaT[i].resize(K,L);
		GaT[i].setZeros();
		workT[i].resize(K,3);
		workT[i].setZeros();
		uT[i].resize(L);
		uT[i].setZeros();
	}
*/

	Timer time, time2;
	time.start();
	srandom(0);
	Vector<int> perm;
	perm.randperm(M);

	Aodd.setZeros();
	Bodd.setZeros();
	Aeven.setZeros();
	Beven.setZeros();
	Aorig.copy(_A);
	Borig.copy(_B);

	// number of iteration ==> why change with param.iter
	int nbIter;
	if (param.iter < 0){
		nbIter = 100000000;
	} else {
		nbIter = param.iter;
	}
	
	bool even = true;
	int last_written = -40;
	int i;
	for (i = 0; i<nbIter; ++i) {
/* VERBOSE
		if (param.verbose) {
			cout << "Iteration: " << i << endl;
			flush(cout);
		}
*/
		time.stop();
		
		// early stop ==> why ?
		if (param.iter < 0 && time.getElapsed() > T(-param.iter)){
			break;
		}
/* LOG
		if (param.log) {
			int seconds=static_cast<int>(floor(log(time.getElapsed())*5));
			if (seconds > last_written) {
				last_written++;
				sprintf(buffer_string, "%s_%d.log", param.logName, last_written+40);
				writeLog(_D,T(time.getElapsed()),i,buffer_string);
				fprintf(stderr,"\r%d",i);
			}
		}
*/
		time.start();
      
		// Gram matrix G = DtD
		Matrix<T> G;
		_D.XtX(G);
		
		if (param.clean){
			this->cleanDict(X, G, param.posD, param.modeD, param.gamma1, param.gamma2);
		}
		// avoid null coeff on the diagonal ?
		G.addDiag(MAX(param.lambda2,1e-10));
		int j;
		
		// A & B initialization
		for (j = 0; j<_NUM_THREADS; ++j) {
			AT[j].setZeros();
			BT[j].setZeros();
		}

#pragma omp parallel for private(j)
		for (j = 0; j<batchsize; ++j) {
#ifdef _OPENMP
			int numT = omp_get_thread_num();
#else
			int numT = 0;
#endif
			const int index = perm[(j+i*batchsize) % M];
			Vector<T>& Xj = XT[numT];
			SpVector<T>& spcoeffj = spcoeffT[numT];
			Vector<T>& DtRj = DtRT[numT];
			//X.refCol(index,Xj);
			X.getData(Xj,index);
			if (param.whiten) {
				Xj.whiten(X.V());
				}
			_D.multTrans(Xj,DtRj);
			Matrix<T>& Gs = GsT[numT];
			Matrix<T>& Ga = GaT[numT];
			Matrix<T>& invGs = invGsT[numT];
			Matrix<T>& work= workT[numT];
			Vector<T>& u = uT[numT];
			Vector<int> ind;
			Vector<T> coeffs_sparse;
			spcoeffj.setL(L);
			spcoeffj.refIndices(ind);
			spcoeffj.refVal(coeffs_sparse);
			T normX=Xj.nrm2sq();
			coeffs_sparse.setZeros();
			if (param.mode < SPARSITY) {
				if (param.mode == PENALTY && param.lambda==0 && param.lambda2 > 0 && !param.posAlpha) {
					Matrix<T>& GG = G;
					u.set(0);
					GG.conjugateGradient(DtRj,u,1e-4,2*K);
					for (int k = 0; k<K; ++k) {
						ind[k]=k;
						coeffs_sparse[k]=u[k];
					}
				} else {
					coreLARS2(DtRj,G,Gs,Ga,invGs,u,coeffs_sparse,ind,work,normX,param.mode,param.lambda,param.posAlpha);
				}
			} else {
				if (param.mode == SPARSITY) {
					coreORMPB(DtRj,G,ind,coeffs_sparse,normX,L,T(0.0),T(0.0));
				} else if (param.mode==L2ERROR2) {
					coreORMPB(DtRj,G,ind,coeffs_sparse,normX,L,param.lambda,T(0.0));
				} else {
					coreORMPB(DtRj,G,ind,coeffs_sparse,normX,L,T(0.0),param.lambda);
				}
			}
			int count2=0;
			for (int k = 0; k<L; ++k){
				if (ind[k] == -1) {
					break;
				} else {
					++count2;
				}
			}
			sort(ind.rawX(),coeffs_sparse.rawX(),0,count2-1);
			spcoeffj.setL(count2);
			AT[numT].rank1Update(spcoeffj);
			BT[numT].rank1Update(Xj,spcoeffj);
		}

		if (param.batch) {
			_A.setZeros();
			_B.setZeros();
			for (j = 0; j<_NUM_THREADS; ++j) {
				_A.add(AT[j]);
				_B.add(BT[j]);
			}
			Vector<T> di, ai,bi;
			Vector<T> newd(n);
			for (j = 0; j<param.iter_updateD; ++j) {
				for (int k = 0; k<K; ++k) {
					if (_A[k*K+k] > 1e-6) {
						_D.refCol(k,di);
						_A.refCol(k,ai);
						_B.refCol(k,bi);
						_D.mult(ai,newd,T(-1.0));
						newd.add(bi);
						newd.scal(T(1.0)/_A[k*K+k]);
						newd.add(di);
						if (param.modeD != L2) {
							newd.sparseProject(di,T(1.0), sparseD, param.gamma1,
												param.gamma2, T(2.0), param.posD);
						} else {
							if (param.posD){
								newd.thrsPos();
							}
							newd.normalize2();
							di.copy(newd);
						}
					} else if (param.clean) {
						_D.refCol(k,di);
						di.setZeros();
					}
				}
			}
		} else if (param.stochastic) {
			_A.setZeros();
			_B.setZeros();
			for (j = 0; j<_NUM_THREADS; ++j) {
				_A.add(AT[j]);
				_B.add(BT[j]);
			}
			_D.mult(_A,_B,false,false,T(-1.0),T(1.0));
			T step_grad=rho/T(t0+batchsize*(i+1));
			_D.add(_B,step_grad);
			Vector<T> dj;
			Vector<T> dnew(n);
			if (param.modeD != L2) {
				for (j = 0; j<K; ++j) {
					_D.refCol(j,dj);
					dnew.copy(dj);
					dnew.sparseProject(dj, T(1.0), sparseD, param.gamma1,
										param.gamma2, T(2.0), param.posD);
				}
			} else {
				for (j = 0; j<K; ++j) {
					_D.refCol(j,dj);
					if (param.posD) dj.thrsPos();
						dj.normalize2();
					}
			}
		} else {
			/// Dictionary Update
			/// Check the epoch parity
			int epoch = (((i+1) % M)*batchsize) / M;
			if ((even && ((epoch % 2) == 1)) || (!even && ((epoch % 2) == 0))) {
				Aodd.copy(Aeven);
				Bodd.copy(Beven);
				Aeven.setZeros();
				Beven.setZeros();
				count=countPrev;
				countPrev=0;
				even=!even;
			}

			int ii=_itercount+i;
			int num_elem=MIN(2*M, ii < batchsize ? ii*batchsize :
					batchsize*batchsize+ii-batchsize);
			T scal2=T(T(1.0)/batchsize);
			T scal;
			int totaliter=_itercount+count;
			if (param.modeParam == PARAM2) {
				scal=param.rho;
			} else if (param.modeParam == PARAM1) {
				scal=MAX(0.95,pow(T(totaliter)/T(totaliter+1),-rho));
			} else {
				scal = T(_itercount+num_elem+1-
                  batchsize)/T(_itercount+num_elem+1);
			}
			Aeven.scal(scal);
			Beven.scal(scal);
			Aodd.scal(scal);
			Bodd.scal(scal);
			if ((_itercount > 0 && i*batchsize < M) 
               || (_itercount == 0 && t0 != 0 && 
                  i*batchsize < 10000)) {
				Aorig.scal(scal);
				Borig.scal(scal);
				_A.copy(Aorig);
				_B.copy(Borig);
			} else {
				_A.setZeros();
				_B.setZeros();
			}
			for (j = 0; j<_NUM_THREADS; ++j) {
				Aeven.add(AT[j],scal2);
				Beven.add(BT[j],scal2);
			}
			_A.add(Aodd);
			_A.add(Aeven);
			_B.add(Bodd);
			_B.add(Beven);
			++count;
			++countPrev;

			Vector<T> di, ai,bi;
			Vector<T> newd(n);
			for (j = 0; j<param.iter_updateD; ++j) {
				for (int k = 0; k<K; ++k) {
				if (_A[k*K+k] > 1e-6) {
					_D.refCol(k,di);
					_A.refCol(k,ai);
					_B.refCol(k,bi);
					_D.mult(ai,newd,T(-1.0));
					newd.add(bi);
					newd.scal(T(1.0)/_A[k*K+k]);
					newd.add(di);
					if (param.modeD != L2) {
						newd.sparseProject(di,T(1.0),sparseD,
                           param.gamma1,param.gamma2,T(2.0),param.posD);
					} else {
						if (param.posD){
							newd.thrsPos();
						}
						newd.normalize2();
						di.copy(newd);
					}
				} else if (param.clean && 
                     ((_itercount+i)*batchsize) > 10000) {
					_D.refCol(k,di);
					di.setZeros();
					}
				}
			}
		}
	}

	_itercount += i;
	if (param.verbose)
		time.printElapsed();
	delete[](spcoeffT);
	delete[](DtRT);
	delete[](AT);
	delete[](BT);
	delete[](GsT);
	delete[](invGsT);
	delete[](GaT);
	delete[](uT);
	delete[](XT);
	delete[](workT);
};