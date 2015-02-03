#ifndef ANALYSIS3_H
#define	ANALYSIS3_H

class Measure2
{
public:
	DBL mass;
	DBL visc;
	DBL volume;
	DBL density;
	VecF2 vAvg;
	DBL vMax;
	DBL vCoeff;
	DBL vErr;
	DBL vMA;
	VecF2 clotPressure;

	Measure2() : mass(0), visc(0), volume(0), density(0), vAvg(VecF20), vMax(0), vCoeff(0), vErr(0), vMA(0), clotPressure(VecF20) {}

	void save(MyWriteFile &f)
	{
		f.write(mass);
		f.write(visc);
		f.write(volume);
		f.write(density);
		f.write(vAvg);
		f.write(vMax);
		f.write(vCoeff);
		f.write(vErr);
		f.write(vMA);
		//f.write(clotPressure);
	}

	void load(MyReadFile &f)
	{
		f.read(mass);
		f.read(visc);
		f.read(volume);
		f.read(density);
		f.read(vAvg);
		f.read(vMax);
		f.read(vCoeff);
		f.read(vErr);
		f.read(vMA);
		//f.read(clotPressure);
	}

private:
};

class Anl
{
public:
	Hist energy;	//energy evolution in time
	Hist plClot;	//total number of platelets in clot
	Hist plCore;	//nuber of platelets in clot core
	Hist heightClot;//height of clot
	Hist heightCore;//height of clot core
	Hist conT;		//total concentration of thrombin
	Hist conFg;		//total concentration of fibrinogen
	Hist conFp;		//total concentration of fibrin polymer
	Hist visc;		//viscosity
	Measure2 m;		//all temporary measurements
	DBL *v1;		//v_x(y)
	DBL *d1;		//d(y)
	DBL *v2x;		//v_x(x,y)
	DBL *v2y;		//v_y(x,y)
	DBL *d2;		//d(x,y)
	DBL *v1ns;		//v_x(y) Navier-Stokes(parabola)
	DBL *v2xPde;	//v_x(x,y) for PDE
	DBL *v2yPde;	//v_y(x,y) for PDE
	//DBL *px;		//pressure on platelets in clot (x direction)
	//DBL *py;		//pressure on platelets in clot (y direction)

	Anl()
	{
		v1		= NULL;
		d1		= NULL;
		v2x		= NULL;
		v2y		= NULL;
		d2		= NULL;
		v1ns	= NULL;
		v2xPde	= NULL;
		v2yPde	= NULL;
		//px		= NULL;
		//py		= NULL;
		anl		= Mat<Stat>();
		//anlPres	= Mat<Stat>();
		anlPde	= Mat<Stat>();
		m		= Measure2();
	}

	Anl(Setup &s)
	{
		energy	= Hist();
		plClot	= Hist();
		plCore	= Hist();
		visc	= Hist(true);		
		heightClot = Hist(true);
		heightCore = Hist(true);
		conT	= Hist(true);
		conFg	= Hist(true);
		conFp	= Hist(true);
		m		= Measure2();

		d		= s.anl.density;
		size	= s.spc.size;
		dx		= size.y / d;

		int nx	= size.x / size.y * d;
		int ny	= d;
		int n	= nx * ny;
		
		anl		= Mat<Stat>(nx, ny);
		anlPde	= Mat<Stat>(nx, ny);
		//anlPres	= Mat<Stat>(nx, ny);

		v1		= new DBL[ny];
		d1		= new DBL[ny];
		v1ns	= new DBL[ny];
		v2x		= new DBL[n];
		v2y		= new DBL[n];
		d2		= new DBL[n];
		v2xPde	= new DBL[n];
		v2yPde	= new DBL[n];
		//px		= new DBL[n];
		//py		= new DBL[n];
	}

	void update(Setup &s, AtomList &al)
	{
		DBL maxClotH = 0;
		DBL maxCoreH = 0;
		DBL coeff = pow(MU::vel, 2);
		for(int i=0; i<al.getSize(); i++)
		{
			if(!al.obj[i].free)
			{
				energy.update(pow(MODUL(al.obj[i].v), 2) * coeff);

				VecF2 ptmp = al.obj[i].pos / dx;
				VecL2 p = VecL2(ptmp.x, ptmp.y);

				anl(p).update(al.obj[i].v);
				anlPde(p).update(al.obj[i].v);
				//if(al.obj[i].inClot) anlPres(p).update(al.obj[i].F);
				
				if(s.plt.enabled && al.obj[i].isSPlt() && al.obj[i].inClot) 
				{
					maxClotH = max(maxClotH, al.obj[i].pos.y);
					plClot.update(1);
					if(al.obj[i].oldPlt)
					{
						plCore.update(1);
						maxCoreH = max(maxCoreH, al.obj[i].pos.y);  
					}
				}

				//TODO: analysis of clot growth for complex platelets (CPlat)
			}
		}

		energy.inc();
		if(s.plt.enabled) 
		{
			plClot.inc();
			plCore.inc();
			heightClot.update(maxClotH / s.spc.size.y);
			heightCore.update(maxCoreH / s.spc.size.y);
		}

		if(s.run.step % s.anl.anlStep == 0) 
			anlVD(s);

		if(s.run.step % s.anl.clotStep == 0) 
		{
			anlHist(s);
		}
	}

	//analize velocity and density profiles for PDE
	void anlPdeVD(Setup &s, PdeSys &ps, CProgressCtrl &p, CStatic &desc)
	{
		/*
		if(pdeTest)
		{
			ps.update(100);
			return;
		}
		*/

		desc.SetWindowTextA("Analyzing data.");
		for(int j=0; j<anlPde.n; j++)
		{
			for(int i=0; i<anlPde.m; i++)
				anlPde(i,j).avg();
		}
		p.StepIt();
	
		DBL maxVx = 0, maxVy = 0;

		desc.SetWindowTextA("Preparing data for graphs.");
		for(int j=0; j<anlPde.n; j++)
		{
			for(int i=0; i<anlPde.m; i++)
			{
				int k = anlPde.index(i,j);
				DBL vx = anlPde(i,j).v.x;
				DBL vy = anlPde(i,j).v.y;
				DBL dd = anlPde(i,j).n;
				v2xPde[k] = (dd!=0) ? vx * MU::vel / dd : 0;
				v2yPde[k] = (dd!=0) ? vy * MU::vel / dd : 0;
				maxVx = max(maxVx, v2xPde[k]);
				maxVy = max(maxVy, fabs(v2yPde[k]));
			}
		}
		p.StepIt();

		ps.update(anlPde, s.uni.u_out / 20, p, desc);

		VecF3 con = ps.getCon();

		conT.update(con.x);
		conFg.update(con.y);
		conFp.update(con.z);
		conT.store();
		conFg.store();
		conFp.store();

		//drawAnlPde(s, maxVx, maxVy);

		resetAnlPde();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::ANL);
		energy.save(f);
		plClot.save(f);
		plCore.save(f);
		visc.save(f);
		m.save(f);
		anl.save(f);
		anlPde.save(f);
		f.write(d);
		f.write(dx);
		f.write(size);
		//anlPres.save(f);
		heightClot.save(f);
		heightCore.save(f);
		conT.save(f);
		conFg.save(f);
		conFp.save(f);
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::ANL)) return;
		energy.load(f);
		plClot.load(f);
		plCore.load(f);
		visc.load(f);
		m.load(f);
		anl.load(f);
		anlPde.load(f);
		f.read(d);
		f.read(dx);
		f.read(size);
		//anlPres	= Mat<Stat>(anl.m, anl.n);
		heightClot.load(f);
		heightCore.load(f);
		conT.load(f);
		conFg.load(f);
		conFp.load(f);

		int nx	= size.x / size.y * d;
		int ny	= d;
		int n	= nx * ny;
		v1		= new DBL[ny];
		d1		= new DBL[ny];
		v1ns	= new DBL[ny];
		v2x		= new DBL[n];
		v2y		= new DBL[n];
		d2		= new DBL[n];
		v2xPde	= new DBL[n];
		v2yPde	= new DBL[n];
		//px		= new DBL[n];
		//py		= new DBL[n];
	}

	void clear()
	{
		free(v1);
		free(d1);
		free(v2x);
		free(v2y);
		free(v1ns);
		free(v2xPde);
		free(v2yPde);
		//free(px);
		//free(py);
		anl.clear();
		anlPde.clear();
		//anlPres.clear();
		energy.reset();
		plClot.reset();
		plCore.reset();
		heightClot.reset();
		heightCore.reset();
		conT.reset();
		conFg.reset();
		conFp.reset();
		visc.reset();
	}

	void txtHist(Setup &s) 
	{
		if(!s.plt.enabled) return;

		string str = toString(s.run.t * MU::time, 6);
		string s1 = s.run.name + "_clot";
		string s2 = s.run.name + "_core";
		string s5 = s.run.name + "_clotHeight";
		string s6 = s.run.name + "_coreHeight";
		string s7 = s.run.name + "_conT";
		string s8 = s.run.name + "_conFg";
		string s9 = s.run.name + "_conFp";
		MyWriteFile f1 = MyWriteFile(FS::getFN(FS::result, s1, "txt"), false);
		MyWriteFile f2 = MyWriteFile(FS::getFN(FS::result, s2, "txt"), false);
		MyWriteFile f5 = MyWriteFile(FS::getFN(FS::result, s5, "txt"), false);
		MyWriteFile f6 = MyWriteFile(FS::getFN(FS::result, s6, "txt"), false);
		MyWriteFile f7 = MyWriteFile(FS::getFN(FS::result, s7, "txt"), false);
		MyWriteFile f8 = MyWriteFile(FS::getFN(FS::result, s8, "txt"), false);
		MyWriteFile f9 = MyWriteFile(FS::getFN(FS::result, s9, "txt"), false);
		
		f1.setSep("");
		f2.setSep("");
		f5.setSep("");
		f6.setSep("");
		f7.setSep("");
		f8.setSep("");
		f9.setSep("");
		for(int i=0; i<plClot.hist.size(); i++) { f1.write(plClot.hist(i));	f1.newLine(); }
		for(int i=0; i<plCore.hist.size(); i++) { f2.write(plCore.hist(i));	f2.newLine(); }
		for(int i=0; i<heightClot.hist.size(); i++) { f5.write(heightClot.hist(i));	f5.newLine(); }
		for(int i=0; i<heightCore.hist.size(); i++) { f6.write(heightCore.hist(i));	f6.newLine(); }
		for(int i=0; i<conT.hist.size(); i++) { f7.write(conT.hist(i));	f7.newLine(); }
		for(int i=0; i<conFg.hist.size(); i++) { f8.write(conFg.hist(i));	f8.newLine(); }
		for(int i=0; i<conFp.hist.size(); i++) { f9.write(conFp.hist(i));	f9.newLine(); }

		f1.close();
		f2.close();
		f5.close();
		f6.close();
		f7.close();
		f8.close();
		f9.close();
	}

	void drawHist(Setup &s) 
	{
		VecF2 time = VecF2(0, s.run.t * MU::time);

		string s1 = s.run.name + "_energy";
		string s2 = s.run.name + "_clot";
		string s3 = s.run.name + "_clotHeight";
		string s4 = s.run.name + "_con";

		Graph2 g = Graph2();
		g.create(FS::result, s1, -1, 1, 1);
		g.setLabel("time [s]", "energy"); g.setData(energy.hist); g.setRange(time, VecF2(0, energy.maxVal * 1.1)); g.plot();
		g.draw();

		if(s.plt.enabled)
		{
			g.create(FS::result, s3, -1, 1, 1);
			g.setLabel("time [s]", "height"); g.setData(heightClot.hist, heightCore.hist); g.setRange(time, VecF2(0, 1)); g.plot();
			g.draw();

			g.create(FS::result, s2, -1, 1, 1);
			g.setLabel("time [s]", "platelets"); g.setData(plClot.hist, plCore.hist); g.setRange(time, VecF2(0, plClot.maxVal)); g.plot();
			g.draw();

			g.create(FS::result, s4, -1, 1, 1);
			DBL maxCon = max(max(conT.maxVal,conFg.maxVal),conFp.maxVal) ;
			g.setLabel("time [s]", "concentration"); g.setData(conT.hist, conFg.hist, conFp.hist); g.setRange(time, VecF2(0, maxCon)); g.plot();
			g.draw();
		}
	}
	
private:
	Mat<Stat>	anl;			//data for velocity and density, for general analysis
	//Mat<Stat>	anlPres;		//data for pressure
	Mat<Stat>	anlPde;			//data for velocity and density, for averaging in time period of DPD-PDE exchange
	DBL			d;				//analysis precission in y direction
	DBL			dx;				//calculated analysis step (size.y / d)
	VecF2		size;

	//analize velocity and density profiles
	void anlVD(Setup &s)
	{
		DBL maxD = 0, maxVx = 0;
		//DBL minPx = 0, minPy = 0, maxPx = 0, maxPy = 0;
		DBL vMax = 0;
		VecF2 vAvg = VecF20;
		//VecF2 pTot = VecF20;
		int count = 0;

		//density multipliers (for 1D and 2D)
		DBL cd1 = s.atom.m / (s.anl.anlStep * s.spc.size.x * dx * s.spc.size.y / s.atom.n.y) * MU::dens;
		DBL cd2 = s.atom.m / (s.anl.anlStep * pow(dx, 2) * s.spc.size.y / s.atom.n.y) * MU::dens;

		//DBL pUnit = MU::mass * s.atom.n.y / (MU::len * pow(MU::time, 2) * s.spc.size.y * dx);

		for(int j=0; j<anl.n; j++)
		{
			v1[j] = 0;
			d1[j] = 0;
			for(int i=0; i<anl.m; i++)
			{
				int k = anlPde.index(i,j);
				DBL vx = anl(i,j).v.x;
				DBL vy = anl(i,j).v.y;
				DBL dd = anl(i,j).n;
				v1[j] += vx;
				d1[j] += dd;
				d2[k] = dd * cd2;
				v2x[k] = dd == 0 ? 0 : vx * MU::vel / dd;
				v2y[k] = dd == 0 ? 0 : vy * MU::vel / dd;
				maxVx = max(maxVx, v2x[k]);
				maxD = max(maxD, d2[k]);
				vAvg += VecF2(vx, vy);
				count += dd;

				//k = anlPres.index(i,j);
				//vx = anlPres(i,j).v.x;
				//vy = anlPres(i,j).v.y;
				//dd = anlPres(i,j).n;
				//px[k] = dd == 0 ? 0 : vx * pUnit / s.anl.anlStep;
				//py[k] = dd == 0 ? 0 : vy * pUnit / s.anl.anlStep;
				//minPx = min(minPx, px[k]);
				//maxPx = max(maxPx, px[k]);
				//minPy = min(minPy, py[k]);
				//maxPy = max(maxPy, py[k]);

				//pTot += VecF2(px[k], py[k]);
			}
			if(d1[j] != 0) v1[j] /= d1[j];
			vMax = max(vMax, v1[j]);
			v1[j] *= MU::vel;
			d1[j] *= cd1;
		}
		if(count != 0) vAvg /= count;

		m.mass = count * s.atom.m * MU::mass / s.anl.anlStep;
		m.volume = s.spc.size.x * s.spc.size.y * s.spc.size.y * MU::vol / s.atom.n.y;
		m.density =  m.mass / m.volume;

		m.vAvg = vAvg * MU::vel;
		m.vMax = vMax * MU::vel;

		m.vMA = m.vMax / m.vAvg.x;
		m.vCoeff = 6. * m.vAvg.x / pow(s.spc.size.y * MU::len, 2);
		m.visc = m.density * s.met.G.x * MU::acc / (m.vCoeff * 2);

		visc.update(m.visc);
		visc.store();

		for(int j=0; j<d; j++)
		{
			DBL y = dx * (0.5 + j) * MU::len;
			v1ns[j] = m.vCoeff * y * (s.spc.size.y * MU::len - y); 
		}

		//drawAnl(s, maxD, maxVx, minPx, minPy, maxPx, maxPy);
		drawAnl(s, maxD, maxVx, 0, 0, 0, 0);
		txtAnl(s);
		
		resetAnl();
		resetAnlPres();
	}

	void anlHist(Setup &s)
	{
		energy.store();
		if(s.plt.enabled) 
		{
			int k = plClot.hist.size();

			plClot.store();
			plCore.store();
			heightClot.store();
			heightCore.store();
			
			if(k > 1)
			{
				//clotGrowthRate.update((plClot.hist(k-1) - plClot.hist(k-2)) / (s.run.dt * MU::time * s.anl.clotStep));
				//coreGrowthRate.update((plCore.hist(k-1) - plCore.hist(k-2)) / (s.run.dt * MU::time * s.anl.clotStep));
				if(plClot.hist(k-2) != plCore.hist(k-2) && plClot.hist(k-1) == plCore.hist(k-1)) s.run.stop = s.run.t + s.anl.clotStep * 10 * s.run.dt;
			}
			else
			{
				//clotGrowthRate.update(plClot.hist(k-1) - 10);
				//coreGrowthRate.update(plCore.hist(k-1));
			}

			//clotGrowthRate.store();
			//coreGrowthRate.store();
		}

		//drawHist(s);
	}

	void resetAnl()
	{
		anl.init(Stat());
	}

	void resetAnlPde()
	{
		anlPde.init(Stat());
	}

	void resetAnlPres()
	{
		//anlPres.init(Stat());
	}

	void drawAnl(Setup &s, DBL maxD, DBL maxVx, DBL minPx, DBL minPy, DBL maxPx, DBL maxPy) 
	{
		VecF2 rx	= VecF2(0, s.spc.size.x * MU::len);
		VecF2 ry	= VecF2(0, s.spc.size.y * MU::len);
		DBL tmp		= fabs(maxD - s.uni.rho_in);
		VecF2 rD	= VecF2(s.uni.rho_in - tmp, s.uni.rho_in + tmp);
		//rD	= VecF2(0, 2000);
		VecF2 rVx	= VecF2(0, maxVx * 1.1);
		DBL pxLim = max(fabs(minPx),fabs(maxPx));
		DBL pyLim = max(fabs(minPy),fabs(maxPy));
		VecF2 rPx	= (minPx < 0 && maxPx > 0) ? VecF2(-pxLim, pxLim) : (minPx < 0) ?  VecF2(-pxLim, 0) : VecF2(0, pxLim);
		VecF2 rPy	= (minPy < 0 && maxPy > 0) ? VecF2(-pyLim, pyLim) : (minPy < 0) ?  VecF2(-pyLim, 0) : VecF2(0, pyLim);

		string s1 = s.run.name + "_vel";
		string s2 = s.run.name + "_visc";

		Graph2 g = Graph2();

		g.create(FS::result, s1, s.run.step, 4, 4);
		g.setLabel("y [m]", "\\rho(y) [kg/m^3]"); g.setData(d1, anl.n, toString(s.uni.rho_in)); g.setRange(ry, rD); 			
		g.plot(VecL3(4,4,0));
		g.setLabel("y [m]", "\\v_x(y) [m/s]"); g.setData(v1, v1ns, anl.n);	g.setRange(ry, rVx); 					
		g.plot(VecL3(4,4,1));
		g.setLabel("x [m]", "y [m]"); g.setData(v2x, v2y, anl.m, anl.n); g.setRange(rx, ry);  
		g.vect(VecL3(1,4,2));
		g.setLabel("x [m]", "y [m]", "\\v_x(x,y) [m/s]"); g.setData(v2x, anl.m, anl.n); g.setRange(rx, ry, rVx); 
		g.dens(VecL3(1,4,1)); 
		g.setLabel("x [m]", "y [m]", "\\v_x(x,y) [m/s]"); g.setData(v2x, anl.m, anl.n); g.setRange(rx, ry, rVx); 
		g.surf(VecL3(4,4,3));
		g.setLabel("x [m]", "y [m]", "\\rho(x,y) [kg/m^3]"); g.setData(d2, anl.m, anl.n); g.setRange(rx, ry, rD); 
		g.dens(VecL3(1,4,3)); 
		g.setLabel("x [m]", "y [m]", "\\rho(x,y) [kg/m^3]"); g.setData(d2, anl.m, anl.n); g.setRange(rx, ry, rD); 
		g.surf(VecL3(4,4,2));
		g.draw();

		VecF2 rT = VecF2(0, s.run.t * MU::time);
		DBL pad = max((visc.maxVal - visc.minVal) * 0.1, 0.001);
		VecF2 rHV = VecF2(visc.minVal - pad, visc.maxVal + pad);

		g.create(FS::result, s2);
		g.setData(visc.hist, toString(s.uni.mu_in)); g.setRange(rT, rHV); g.setLabel("time", "\\mu"); g.plot();
		g.draw();
	}

	void txtAnl(Setup &s) 
	{
		string s1 = s.run.name + "_vel";
		bool append = (s.run.step / s.anl.anlStep > 1);
		string file = FS::getFN(FS::root, s1, "csv");
		MyWriteFile f(file, append);
		if(!append)
		{
			f.setSep("");
			f.writeTag("sep=;");
			f.newLine();
			f.write("time");
			f.setSep(";");
			f.write("step");
			f.write("energy");
			f.write("density");
			f.write("maxVx");
			f.write("avgVx");
			f.write("maxVx/avgVx");
			f.write("avgVy");
			f.write("current visc");
			f.write("overall visc");
			if(s.plt.enabled) 
			{
				f.write("clot");
				f.write("clot core");
				f.write("clot height");
				f.write("core height");
			}
			f.newLine();
		}
			
		f.setSep("");
		f.write(s.run.t);
		f.setSep(";");
		f.write(s.run.step);
		f.write(energy.getLast());
		f.write(m.density);
		f.write(m.vMax);
		f.write(m.vAvg.x);
		f.write(m.vMA);
		f.write(m.vAvg.y);
		f.write(m.visc);
		f.write(visc.getLast(10));
		if(s.plt.enabled) 
		{
			f.write(plClot.getLast());
			f.write(plCore.getLast());
			f.write(heightClot.getLast());
			f.write(heightCore.getLast());
		}
		f.newLine();

		f.close();
	}
};

#endif