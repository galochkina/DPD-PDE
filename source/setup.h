#ifndef SETUP_H
#define	SETUP_H

//*******************************************************************************************************************************************************
class SetupAtoms
{
public:
	VecL2 n;		//number of atoms (x*y)
	DBL r;			//physical radius of each atom
	DBL rShow;		//drawing radius
	DBL m;			//total m of all atoms
	DBL v_maxPF;	//maximum velocity in x direction in Poiseuille flow (used to imply initial velocity profile)

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SATM);
		f.write(n.x);
		f.write(n.y);
		f.write(r);
		f.write(m);
		f.write(v_maxPF);
		f.write(rShow);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SATM)) return;
		f.read(n.x);
		f.read(n.y);
		f.read(r);
		f.read(m);
		f.read(v_maxPF);
		f.read(rShow);
		f.nextLine();
	}
};

//*******************************************************************************************************************************************************
class Analysis
{
public:

	int	density;	//density in y direction		
	int	anlStep;	//number of steps taken into account for measurements(for general analysis)
	int clotStep;	//platelets in clot count memorizing step

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SANL);
		f.write(density);
		f.write(anlStep);
		f.write(clotStep);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SANL)) return;
		f.read(density);
		f.read(anlStep);
		f.read(clotStep);
		f.nextLine();
	}
};

//*******************************************************************************************************************************************************
class SetupMethod
{
public:		
	VecF2	G;					//total body F (body F f = G / (size.x * size.y))
	DBL		maxV;				//maximum allowed speed

	//DPD
	DBL		gamma;				//F function gamma			(fType == 1)
	DBL		sigma;				//F function sigma			(fType == 1)
	DBL		TkB;				//F function T * kB			(fType == 1)
	DBL		rc;					//F function rc				(fType == 1)
	DBL		p;					//F function p				(fType == 1)
	DBL		A;					//F function coefficients	(fType == 1)
	DBL		sigmaSqrtdt;		//sigma / sqrt(dt)

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SMET);
		f.write(G.x);
		f.write(G.y);
		f.write(maxV);
		f.write(gamma);
		f.write(sigma);
		f.write(TkB);
		f.write(rc);
		f.write(p);
		f.write(A);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SMET)) return;
		f.read(G.x);
		f.read(G.y);
		f.read(maxV);
		f.read(gamma);
		f.read(sigma);
		f.read(TkB);
		f.read(rc);
		f.read(p);
		f.read(A);
		f.nextLine();
	}
};

//*******************************************************************************************************************************************************
class SetupSpace
{
public:
	VecF2	size;

	DBL		maxR;					//maximum radius for mesh

	//for preparation area
	bool	prepArea;
	DBL		sizePrepX;				//how long preparation area should be

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SSPC);
		f.write(size.x);
		f.write(size.y);
		f.write(maxR);
		f.write(prepArea);
		f.write(sizePrepX);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SSPC)) return;
		f.read(size.x);
		f.read(size.y);
		f.read(maxR);
		f.read(prepArea);
		f.read(sizePrepX);
		f.nextLine();
	}
};

//*******************************************************************************************************************************************************
class SetupRun
{
public:
	DBL		t;
	long	step;
	long	imgCount;
	DBL		dt;
	DBL		dt2;			//dt for spring forces (dt2=dt/n, where n is int)
	DBL		stop;
	string	name;			//name of simulation
	DBL		excT;			//data exchange between DPD and PDE
	DBL		pdeT;

	bool	newInstance;

	int		showStep;
	int		saveStep;
	int		txtStep;
	int		bmpStep;

	bool	preSimulation;
	bool	preAnimOver;
	DBL		preAnimTime;

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SRUN);
		f.write(t);
		f.write(step);
		f.write(dt);
		f.write(imgCount);
		f.write(stop);
		f.write(name);
		f.write(showStep);
		f.write(saveStep);
		f.write(txtStep);
		f.write(bmpStep);
		f.write(preSimulation);
		f.write(preAnimOver);
		f.write(preAnimTime);
		f.write(excT);
		f.write(pdeT);
		f.write(dt2);
		f.write(newInstance);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SRUN)) return;
		f.read(t);
		f.read(step);
		f.read(dt);
		f.read(imgCount);
		f.read(stop);
		f.read(name);
		f.read(showStep);
		f.read(saveStep);
		f.read(txtStep);
		f.read(bmpStep);
		f.read(preSimulation);
		f.read(preAnimOver);
		f.read(preAnimTime);
		f.read(excT);
		f.read(pdeT);
		f.read(dt2);
		f.read(newInstance);
		f.nextLine();
	}
};

//*******************************************************************************************************************************************************
class SetupPlatelet
{
public:
	bool	enabled;
	DBL		percent;
	DBL		F2weak;
	DBL		F2strong;
	DBL		F2t;			//time to pass from F2weak to F2strong
	DBL		dc;
	DBL		dd;
	bool	constInflow;
	DBL		critFib;
	DBL		F2fib;

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SPLT);
		f.write(enabled);
		f.write(percent);
		f.write(dc);
		f.write(dd);
		f.write(F2weak);
		f.write(F2strong);
		f.write(F2t);
		f.write(constInflow);
		f.write(critFib);
		f.write(F2fib);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SPLT)) return;
		f.read(enabled);
		f.read(percent);
		f.read(dc);
		f.read(dd);
		f.read(F2weak);
		f.read(F2strong);
		f.read(F2t);
		f.read(constInflow);
		f.read(critFib);
		f.read(F2fib);
		f.nextLine();
	}
};



//*******************************************************************************************************************************************************
class SetupUnits
{
public:
	DBL rho_in, mu_in, D_in, u_in, m_in, G_in, r_in, L_in;
	DBL rho_out, mu_out, D_out, u_out, m_out, G_out, r_out, L_out;
	int Nx, Ny;
	int mFrom, mTo, mFromOth, mToOth;
	int gFrom, gTo, gFromOth, gToOth;
	int sFrom, sTo, sFromOth, sToOth;

	SetupUnits()
	{
		rho_in = 1000;
		mu_in = 0.001;
		D_in = 0.00001;
		L_in = 0.00003;
		u_in = 0.001;
		Ny = 30;
		m_in = 0;
		G_in = 0;
		r_in = 0;

		Nx = 0;
		rho_out = 0;
		mu_out = 0;
		D_out = 0;
		L_out = 0;
		u_out = 0;
		m_out = 0;
		G_out = 0;
		r_out = 0;

		mFrom = mFromOth = gFromOth = sFrom = sFromOth = Const::NONE;
		gFrom = Const::KILO;
		mTo = mToOth = Const::MICRO;
		gTo = gToOth = Const::PICO;
		sTo = sToOth = Const::MILI;
	}

	DBL getMultiplier(int m, int g, int s)
	{
		return pow(10, (DBL)(m*(mTo-mFrom) + g*(gTo-gFrom) + s*(sTo-sFrom)));
	}

	DBL getDensityM()	{ return getMultiplier(-3, 1, 0); }
	DBL getVolumeM()	{ return getMultiplier(3, 0, 0); }
	DBL getVelocityM()	{ return getMultiplier(1, 0, -1); }
	DBL getMassM()		{ return getMultiplier(0, 1, 0); }
	DBL getLengthM()	{ return getMultiplier(1, 0, 0); }
	DBL getAreaM()		{ return getMultiplier(2, 0, 0); }
	DBL getTimeM()		{ return getMultiplier(0, 0, 1); }
	DBL getAccM()		{ return getMultiplier(1, 0, -2); }
	DBL getViscM()		{ return getMultiplier(-1, 1, -1); }

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SUNI);
		f.write(rho_in);
		f.write(mu_in);
		f.write(D_in);
		f.write(L_in);
		f.write(u_in);
		f.write(m_in);
		f.write(r_in);
		f.write(G_in);
		f.write(Ny);
		f.write(Nx);
		f.write(rho_out);
		f.write(mu_out);
		f.write(D_out);
		f.write(L_out);
		f.write(u_out);
		f.write(m_out);
		f.write(r_out);
		f.write(G_out);
		f.write(mFrom);
		f.write(mTo);
		f.write(mFromOth);
		f.write(mToOth);
		f.write(gFrom);
		f.write(gTo);
		f.write(gFromOth);
		f.write(gToOth);
		f.write(sFrom);
		f.write(sTo);
		f.write(sFromOth);
		f.write(sToOth);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::SUNI)) return;
		f.read(rho_in);
		f.read(mu_in);
		f.read(D_in);
		f.read(L_in);
		f.read(u_in);
		f.read(m_in);
		f.read(r_in);
		f.read(G_in);
		f.read(Ny);
		f.read(Nx);
		f.read(rho_out);
		f.read(mu_out);
		f.read(D_out);
		f.read(L_out);
		f.read(u_out);
		f.read(m_out);
		f.read(r_out);
		f.read(G_out);
		f.read(mFrom);
		f.read(mTo);
		f.read(mFromOth);
		f.read(mToOth);
		f.read(gFrom);
		f.read(gTo);
		f.read(gFromOth);
		f.read(gToOth);
		f.read(sFrom);
		f.read(sTo);
		f.read(sFromOth);
		f.read(sToOth);
		f.nextLine();
	}

};

//*******************************************************************************************************************************************************
class SetupOptions
{
public:
	int textSize;	//text size
	int textSpace;	//space between rows
	int	proc;		//number of processors
	int	imgW;		//image width in pixels
	string outputPath;	//path to the output folder

	COLORREF cpla;
	COLORREF cplt;
	COLORREF cfib;
	COLORREF cf1;
	COLORREF cf2;

	bool showPla;
	bool showPlt;
	bool showFib;
	bool showF1;
	bool showF2;

	void save(MyWriteFile &f)
	{
		f.write(textSize);
		f.write(textSpace);
		f.write(proc);
		f.write(showPla);
		f.write(showPlt);
		f.write(showFib);
		f.write(showF1);
		f.write(showF2);
		f.write(cpla);
		f.write(cplt);
		f.write(cfib);
		f.write(cf1);
		f.write(cf2);
		f.write(imgW);
		f.write(outputPath);
	}

	void load(MyReadFile &f)
	{
		f.read(textSize);
		f.read(textSpace);
		f.read(proc);
		f.read(showPla);
		f.read(showPlt);
		f.read(showFib);
		f.read(showF1);
		f.read(showF2);
		f.read(cpla);
		f.read(cplt);
		f.read(cfib);
		f.read(cf1);
		f.read(cf2);
		f.read(imgW);
		f.read(outputPath);
	}
};

//*******************************************************************************************************************************************************
class SetupFib
{
public:
	bool enabled;
	DBL dx;
	DBL dt;
	DBL time;
	DBL alpha1;		//thrombin(T) diffusion term coefficient
	DBL alpha2;		//fibrinogen(Fg) diffusion term coefficient
	DBL gamma3;		//thrombin(T) destruction term coefficient
	DBL T0;			//a constant in thrombin reaction term
	DBL C0;			//initial quantity of factor II and thrombin (C0 = II = T)
	DBL k1;			//thrombin(II->T) reaction term coefficient
	DBL k3;			//fibrin polymer(Fg->Fp) reaction term coefficient
	DBL gamma1;
	DBL alpha3;
	DBL gamma2;
	DBL beta1;
	DBL beta2;
	DBL beta3;
	DBL alpha4;
	DBL woundBeg;
	DBL woundEnd;
	DBL woundVal;
	DBL initFg;

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::FIBR);
		f.write(enabled);
		f.write(dx);
		f.write(dt);
		f.write(time);
		f.write(woundBeg);
		f.write(woundEnd);
		f.write(woundVal);
		f.write(gamma1);
		f.write(alpha1);
		f.write(alpha2);
		f.write(T0);
		f.write(C0);
		f.write(k1);
		f.write(k3);
		f.write(initFg);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::FIBR)) return;
		f.read(enabled);
		f.read(dx);
		f.read(dt);
		f.read(time);
		f.read(woundBeg);
		f.read(woundEnd);
		f.read(woundVal);
		f.read(gamma1);
		f.read(alpha1);
		f.read(alpha2);
		f.read(T0);
		f.read(C0);
		f.read(k1);
		f.read(k3);
		f.read(initFg);
		f.nextLine();
	}
};

//*******************************************************************************************************************************************************
class Setup
{
public:
	SetupSpace		spc;
	SetupMethod		met;
	SetupAtoms		atom;
	SetupRun		run;
	Analysis		anl;
	SetupPlatelet	plt;
	SetupPlatelet	cplt;
	SetupOptions	opt;
	SetupUnits		uni;
	SetupFib		fib;

	void save(MyWriteFile &f)
	{
		spc.save(f);
		met.save(f);
		atom.save(f);
		plt.save(f);
		anl.save(f);
		run.save(f);
		uni.save(f);
		fib.save(f);
	}

	void load(MyReadFile &f)
	{
		spc.load(f);
		met.load(f);
		atom.load(f);
		plt.load(f);
		anl.load(f);
		run.load(f);
		uni.load(f);
		fib.load(f);
	}

};


#endif