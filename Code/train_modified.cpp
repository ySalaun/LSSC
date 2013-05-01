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

	int JJ;
	if (param.iter < 0){
		JJ = 100000000;
	} else {
		JJ = param.iter;
	}
	bool even = true;
	int last_written = -40;
	int i;
	for (i = 0; i<JJ; ++i) {
		if (param.verbose) {
			cout << "Iteration: " << i << endl;
			flush(cout);
		}
		time.stop();
		if (param.iter < 0 && time.getElapsed() > T(-param.iter)){
			break;
		}
		if (param.log) {
			int seconds=static_cast<int>(floor(log(time.getElapsed())*5));
			if (seconds > last_written) {
				last_written++;
				sprintf(buffer_string, "%s_%d.log", param.logName, last_written+40);
				writeLog(_D,T(time.getElapsed()),i,buffer_string);
				fprintf(stderr,"\r%d",i);
			}
		}
		time.start();
      
		Matrix<T> G;
		_D.XtX(G);
		if (param.clean){
			this->cleanDict(X, G, param.posD, param.modeD, param.gamma1, param.gamma2);
		}
		G.addDiag(MAX(param.lambda2,1e-10));
		int j;
		for (j = 0; j<_NUM_THREADS; ++j) {
			AT[j].setZeros();
			BT[j].setZeros();
		}

#pragma omp parallel for private(j)
		for (j = 0; j<batchsize; ++j) {
#ifdef _OPENMP
			int numT=omp_get_thread_num();
#else
			int numT=0;
#endif
			const int index=perm[(j+i*batchsize) % M];
			Vector<T>& Xj = XT[numT];
			SpVector<T>& spcoeffj = spcoeffT[numT];
			Vector<T>& DtRj = DtRT[numT];
			//X.refCol(index,Xj);
			X.getData(Xj,index);
			if (param.whiten) {
				if (param.pattern) {
					Vector<T> mean(4);
					Xj.whiten(mean,param.pattern);
				} else {
					Xj.whiten(X.V());
				}
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