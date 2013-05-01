	template <typename T>
	inline void Image<T>::denoise6(const T sigma, Matrix<T>& D, const Image<T>& orig,
		const int n,
		const T C1b, const T C2, const T J, const T J2, const T thrs,
		const int window,
		const int num_threads, const int iter) {
		Image<T> noisy;
	bool pattern=false;
	T C1=abs<T>(C1b);
	const T C = C2;
	std::ifstream f;
	static T param_C[10000];
	char str[13];
	sprintf(str,"misc/C%d.data",static_cast<int>(C));
	f.open(str);
	if (!f) {
		sprintf(str,"../misc/C%d.data",static_cast<int>(C));
		f.open(str);
	}
	if (!f) {
		printf("Problem while opening C file\n");
		return;
	}
	for (int i = 0; i<10000; ++i) {
		f >> param_C[i];
	}
	f.close();
	static T param_C1[10000];
	char str[13];
	sprintf(str,"misc/C%d.data",static_cast<int>(C1));
	f.open(str);
	if (!f) {
		sprintf(str,"../misc/C%d.data",static_cast<int>(C1));
		f.open(str);
	}
	if (!f) {
		printf("Problem while opening C file\n");
		return;
	}
	for (int i = 0; i<10000; ++i) {
		f >> param_C1[i];
	}
	f.close();
	noisy.copy(*this);
	this->setDataVariables(n);
	const T eps = param_C1[MIN(9999,_m)]*param_C1[MIN(9999,_m)]*sigma*sigma*_m;
	const T eps1 = eps;
	const int batch_size=static_cast<int>(floor(256*((num_threads-1)/2+1)));
	int JJ = 0;
	if (J != 0) {
		if (J > 0)  {
			JJ = static_cast<int>(ceil(J*static_cast<T>(_n)/batch_size))+1;
			printf("Number of patches processed : %d\n",JJ*batch_size);
		} else {
			JJ = static_cast<int>(floor(J));
		}
	}
	Trainer<T> trainer(D,eps,batch_size,num_threads);
	if (JJ != 0)
		trainer.train(*this,JJ,L2ERROR,true,1,pattern);
	Matrix<T> Dold;
	Dold.copy(D);
	trainer.getD(D);
	this->denoise4(sigma,D,orig,n,eps1,num_threads,pattern);
	this->setDataVariables(n);
	for (int i = 0; i<iter; ++i)
		this->denoise6b(sigma,D,noisy,orig,trainer,n,J2,param_C,thrs,window,num_threads,pattern);