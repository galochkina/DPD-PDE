#ifndef PDESYS_H
#define	PDESYS_H

class PdeSys
{
public:
	Pde1 pdeT;	//thrombin
	Pde2 pdeFg;	//fibrinogen
	Pde3 pdeFp;	//fibrin polymer
	Pde4 pdeTf;	//fibrinogen
	Pde5 pdeTm;	//fibrin polymer
	DBL critU;

	PdeSys()
	{
		pdeT	= Pde1();
		pdeFg	= Pde2();
		pdeFp	= Pde3();
		pdeTf	= Pde4();
		pdeTm	= Pde5();
	}

	PdeSys(Setup &s)
	{
		pdeT	= Pde1(s);
		pdeFg	= Pde2(s);
		pdeFp	= Pde3(s);
		pdeTf	= Pde4(s);
		pdeTm	= Pde5(s);
		critU = s.plt.critFib;
	}

	void update(DBL coeff)
	{
		pdeT.update(coeff);
	}

	void update(VecF2 v)
	{
		pdeT.update(v);
	}

	void update(Mat<Stat> &M, DBL ignoreVal, CProgressCtrl &p, CStatic &desc)
	{
		pdeT.update(M, ignoreVal, p, desc);
	}

	void nextStep(Setup &s, bool closedWound)
	{
		//pdeT.nextStep(s, s.run.t < 100 ? false : true);
		pdeTf.nextStep(closedWound);
		pdeTm.nextStep(closedWound);
		pdeT.nextStep(closedWound, *(pdeTm.X), *(pdeTf.X));
		pdeFg.nextStep(pdeT.V, *(pdeT.X), *(pdeFg.X), pdeFp.X);	
		pdeFp.nextStep(*(pdeT.X), *(pdeFg.X));
	}

	bool isOld(VecF2 &p)
	{	
		if(pdeFp.getU(p) > critU) return true;
		return false;
	}

	bool oscillations()
	{
		for(int i=0; i<pdeT.X->m; i++)
			for(int j=0; j<pdeT.X->n; j++)
				if((*(pdeT.X))(i,j) < 0) return true;
		return false;
	}

	void draw(Setup &s, string name)
	{
		int step = s.run.step;
		int m = pdeT.V.m;
		int n = pdeT.V.n;
		VecF2 rx = VecF2(0, pdeT.size.x * MU::len);
		VecF2 ry = VecF2(0, pdeT.size.y * MU::len);

		int N = m*n;
		DBL *vx = new DBL[N];
		DBL *vy = new DBL[N];

		DBL maxX = 0, minX = 0, maxY = 0;
		for(int i=0; i<N; i++)
		{
			DBL t = pdeT.V.v[i].x * MU::vel;
			vx[i] = t;
			maxX = max(maxX, t);
			minX = min(minX, t);

			t = pdeT.V.v[i].y * MU::vel;
			vy[i] = t;
			maxY = max(maxY, fabs(t));
		}


		DBL* T = pdeT.X->arr();
		DBL* Fg = pdeFg.X->arr();
		DBL* Fp = pdeFp.X.arr();

		DBL minT = -0.1;
		DBL maxT = 1.1;

		string s1 = s.run.name + "_pdeD";
		string s2 = s.run.name + "_pdeS";

		Graph2 g = Graph2();
		g.create(FS::result, s1, step, pdeT.size.x / pdeT.size.y, 4);
		g.setLabel("x [m]", "y [m]"); g.setData(vx, m, n); g.setRange(rx, ry, VecF2(-minX, maxX)); g.dens(VecL3(1,5,0));
		g.setLabel("x [m]", "y [m]"); g.setData(vy, m, n); g.setRange(rx, ry, VecF2(-maxY, maxY)); g.dens(VecL3(1,5,1));
		g.setLabel("x [m]", "y [m]"); g.setData(T, m, n); g.setRange(rx, ry, VecF2(minT,maxT)); g.dens(VecL3(1,5,2));
		g.setLabel("x [m]", "y [m]"); g.setData(Fg, m, n); g.setRange(rx, ry, VecF2(minT,maxT)); g.dens(VecL3(1,5,3));
		g.setLabel("x [m]", "y [m]"); g.setData(Fp, m, n); g.setRange(rx, ry, VecF2(minT,maxT)); g.dens(VecL3(1,5,4));
		g.draw();

		g.create(FS::result, s2, step, 2, 2);
		g.setLabel("x [m]", "y [m]", "\\v_x(x,y) [m/s]"); g.setData(vx, m, n); g.setRange(rx, ry, VecF2(-minX, maxX)); g.surf(VecL3(2,3,0), false);
		g.setLabel("x [m]", "y [m]", "\\v_y(x,y) [m/s]"); g.setData(vy, m, n); g.setRange(rx, ry, VecF2(-maxY, maxY)); g.surf(VecL3(2,3,1), false);
		g.setLabel("x [m]", "y [m]", "T(x,y)"); g.setData(T, m, n); g.setRange(rx, ry, VecF2(minT,maxT)); g.surf(VecL3(2,3,2), false);
		g.setLabel("x [m]", "y [m]", "\\F_g(x,y)"); g.setData(Fg, m, n); g.setRange(rx, ry, VecF2(minT,maxT)); g.surf(VecL3(2,3,3), false);
		g.setLabel("x [m]", "y [m]", "\\F_p(x,y)"); g.setData(Fp, m, n); g.setRange(rx, ry, VecF2(minT,maxT)); g.surf(VecL3(2,3,4), false);
		g.draw();

		delete[] vx;
		delete[] vy;
		delete[] T;
		delete[] Fg;
		delete[] Fp;
	}

	void txt(Setup &s, string name)
	{
		string s1 = s.run.name + "_pdeTxt_thr";
		string s2 = s.run.name + "_pdeTxt_fib";
		pdeT.txt(s.run.step, s1);
		pdeFg.txt(s.run.step, s2); 
	}

	void free()
	{
		pdeT.free();
		pdeFg.free();
	}

	void save(MyWriteFile &f)
	{
		pdeT.save(f);
		f.newLine();
		pdeFg.save(f);
		f.newLine();
		pdeFp.save(f);
		f.newLine();
		pdeTf.save(f);
		f.newLine();
		pdeTm.save(f);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(f.isTag(Const::PDE1)) pdeT.load(f);
		f.nextLine();
		if(f.isTag(Const::PDE2)) pdeFg.load(f);
		f.nextLine();
		if(f.isTag(Const::PDE3)) pdeFp.load(f);
		f.nextLine();
		if(f.isTag(Const::PDE4)) pdeFp.load(f);
		f.nextLine();
		if(f.isTag(Const::PDE5)) pdeFp.load(f);
		f.nextLine();
	}

	void paramHotChange(Setup &s)
	{
		critU = s.plt.critFib;
		pdeT.paramHotChange(s);
		pdeFg.paramHotChange(s);
		pdeFp.paramHotChange(s);
		pdeTf.paramHotChange(s);
		pdeTm.paramHotChange(s);
	}

	VecF3 getCon()
	{
		return VecF3(pdeT.getTotalCon(), pdeFg.getTotalCon(), pdeFp.getTotalCon());
	}


private:
};


#endif