#ifndef SAVE_H
#define	SAVE_H

namespace IO
{	
	static void saveDPDStp(Setup &s, string filePath)
	{
		MyWriteFile f = MyWriteFile(filePath, false);

		f.writeTag(Const::SAVT);
		f.write(Const::STP);
		f.newLine();

		bool tmp = s.run.newInstance;
		s.run.newInstance = true;
		s.save(f);			//Setup
		s.run.newInstance = tmp;

		f.close();
	}

	static void saveDPD(Setup &s, string filePath, AtomList &al, PPL &platPairs, Clock &c, PdeSys &ps, Anl &anl, CProgressCtrl &bar, CStatic &desc)
	{
		MyWriteFile f = MyWriteFile(filePath, false);

		f.writeTag(Const::SAVT);
		f.write(Const::S04);
		f.newLine();

		desc.SetWindowTextA("Setup");
		s.save(f);		//Setup
		c.save(f);		//Clock
		bar.StepIt();
		desc.SetWindowTextA("Particles");
		al.save(f);		//Atoms
		bar.StepIt();
		desc.SetWindowTextA("Complex particles");
		bar.StepIt();

		//Fibrin
		desc.SetWindowTextA("PDEs");
		if(s.fib.enabled) ps.save(f);
		bar.StepIt();
		desc.SetWindowTextA("Platelet connections");
		if(s.plt.enabled) platPairs.save(f);
		bar.StepIt();
		desc.SetWindowTextA("analysis");
		anl.save(f);
		bar.StepIt();
		
		f.close();
	}

	static bool loadDPD(Setup &s, string fileName, AtomList &al, PPL &platPairs, Clock &c, PdeSys &ps, Anl &anl, CProgressCtrl &bar, CStatic &desc)
	{
		MyReadFile f = MyReadFile(fileName);

		if(!f.isReadable()) return false;

		//header
		if(f.isTag(Const::SAVT)) 
		{
			string saveType = f.getStr();
			if(saveType.compare(Const::STP)==0)
			{
				f.nextLine();		
				s.load(f);		//Setup
				s.run.newInstance = true;
			}
			else if(saveType.compare(Const::S04)==0 || saveType.compare(Const::S03)==0) 
			{
				f.nextLine();
		
				desc.SetWindowTextA("Setup");
				s.load(f);		//Setup
				c.load(f);		//Clock		
				bar.StepIt();
				desc.SetWindowTextA("Particles");
				al.load(f);		//AtomList
				bar.StepIt();
			
				desc.SetWindowTextA("PDEs");
				if(s.fib.enabled)
				{	
					ps.free();
					ps = PdeSys(s);
					ps.load(f);
				}
				bar.StepIt();

				desc.SetWindowTextA("Platelet connections");
				if(s.plt.enabled) platPairs.load(f);
				bar.StepIt();

				desc.SetWindowTextA("Analysis");
				anl.load(f);
				bar.StepIt();
				s.run.newInstance = false;
			}
			else return false;
		}
		else return false;

		return true;
	}

	/*
	static void PaintBMP(Setup &s, CMyFrame *fr, bool all, string fileName)
	{
		CWaitCursor WaitCursor; 
		//MessageBeep(0xFFFFFFFF); 
		CRect rect; 
		fr->GetClientRect(rect); 

		int xBeg = 0, yBeg = 0;
		int x=0, y=0;
		if(all)
		{
			x = rect.right;
			y = rect.bottom;
		}

		CDC *WinDC;
		HDC CopyDC;
		CBitmap bm;

		// Use either NULL for the screen or hwnd for a window
		WinDC = fr->GetDC();
		CopyDC = CreateCompatibleDC(*WinDC);

		bm.CreateBitmap(x, y, 1, 32, NULL);
		SelectObject (CopyDC, bm);
		            
		BitBlt (CopyDC,
			0,0,
			x, y,
			*WinDC,
			xBeg,yBeg,
			SRCCOPY);
		            
		fr->ReleaseDC(WinDC);
		DeleteDC(CopyDC);

		CImage image;
		image.Attach(bm);
		image.Save(fileName.c_str());
	}

	//=====================================================================
	static void PaintBMP(Setup &s, CMyFrame *fr, bool all)
	{
		PaintBMP(s, fr, all, FS::getFN(FS::image, "img", "png", s.run.step));
	}
	*/

	//=====================================================================
	static void saveProps(Setup &s)
	{
		string fname = "config.txt";
		MyWriteFile f(fname, false);
		s.opt.save(f);
		f.close();
	}

	//=====================================================================
	static bool readProps(Setup &s)
	{
		string fname = "config.txt";
		MyReadFile f = MyReadFile(fname, false);
		s.opt.load(f);
		return true;
	}

	//=====================================================================
	static void txtParams(Setup &s, string folder = "")
	{
		string fol = folder.compare("")==0 ? FS::root : folder;
		string file = FS::getFN(fol, s.run.name + "_par", "csv");
		MyWriteFile f(file, false, 0, true);
		
		f.writeTag("sep="); f.newLine();

		//NAME
		f.writeTag("name"); f.write(s.run.name); f.newLine(); f.newLine();

		//UNITS
		f.writeTag("UNITS"); f.newLine();
		f.writeTag("length"); f.write(s.uni.mTo); f.newLine();
		f.writeTag("mass"); f.write(s.uni.gTo); f.newLine();
		f.writeTag("time"); f.write(s.uni.sTo); f.newLine(); f.newLine();

		//PHYSICAL
		f.writeTag("PHYSICAL"); f.newLine();
		f.writeTag("param"); f.write("SI units"); f.write("sim units"); f.newLine();
		f.writeTag("D"); f.write(s.uni.D_in); f.write(s.uni.D_out); f.newLine();
		f.writeTag("L"); f.write(s.uni.L_in); f.write(s.uni.L_out); f.newLine();
		f.writeTag("rho"); f.write(s.uni.rho_in); f.write(s.uni.rho_out); f.newLine();
		f.writeTag("mu"); f.write(s.uni.mu_in); f.write(s.uni.mu_out); f.newLine();
		f.writeTag("u_avg"); f.write(s.uni.u_in); f.write(s.uni.u_out); f.newLine();
		f.writeTag("Nx"); f.write(s.atom.n.x); f.newLine();
		f.writeTag("Ny"); f.write(s.atom.n.y); f.newLine();
		f.writeTag("r"); f.write(s.uni.r_in); f.write(s.uni.r_out); f.newLine();
		f.writeTag("m"); f.write(s.uni.m_in); f.write(s.uni.m_out); f.newLine();
		f.writeTag("G"); f.write(s.uni.G_in); f.write(s.uni.G_out); f.newLine(); 
		f.writeTag("gen area"); if(s.spc.prepArea) { f.write(s.spc.sizePrepX); } else { f.write("no"); } f.newLine(); f.newLine();

		//DPD
		f.writeTag("DPD"); f.newLine();
		f.writeTag("gamma (FD)"); f.write(s.met.gamma); f.newLine();
		f.writeTag("sigma (FR)"); f.write(s.met.sigma); f.newLine();
		f.writeTag("rc"); f.write(s.met.rc); f.newLine();
		f.writeTag("k"); f.write(s.met.p); f.newLine();
		f.writeTag("a (FC)"); f.write(s.met.A); f.newLine(); f.newLine();

		//ANALYSIS
		f.writeTag("ANALYSIS"); f.newLine();
		f.writeTag("density"); f.write(s.anl.density); f.newLine();
		f.writeTag("step anl"); f.write(s.anl.anlStep); f.newLine();
		f.writeTag("step clot"); f.write(s.anl.clotStep); f.newLine(); f.newLine();

		//RUN
		f.writeTag("RUN"); f.newLine();
		f.writeTag("dt"); f.write(s.run.dt); f.newLine();
		f.writeTag("stop time"); f.write(s.run.stop); f.newLine();
		f.writeTag("dTau"); f.write(s.run.excT); f.newLine(); f.newLine();
		f.writeTag("u_max"); f.write(s.atom.v_maxPF); f.newLine(); f.newLine();

		//OUTPUT
		f.writeTag("OUTPUT"); f.newLine();
		f.writeTag("save"); f.write(s.run.saveStep); f.newLine();
		f.writeTag("img"); f.write(s.run.bmpStep); f.newLine();
		f.writeTag("show"); f.write(s.run.showStep); f.newLine(); f.newLine();
		
		//PLATELET
		if(s.plt.enabled)
		{
			f.writeTag("PLATELET"); f.newLine();
			f.writeTag("amount coeff"); f.write(s.plt.percent); f.newLine();
			f.writeTag("dc"); f.write(s.plt.dc); f.newLine();
			f.writeTag("dd"); f.write(s.plt.dd); f.newLine();
			f.writeTag("const inflow"); f.write(s.plt.constInflow ? "Yes" : "No"); f.newLine();
			f.writeTag("F1"); f.write(s.plt.F2weak); f.newLine();
			f.writeTag("F2"); f.write(s.plt.F2strong); f.newLine();
			f.writeTag("F3"); f.write(s.plt.F2fib); f.newLine();
			f.writeTag("time"); f.write(s.plt.F2t); f.newLine();
			f.writeTag("fib crit"); f.write(s.plt.critFib); f.newLine(); f.newLine();
		}

		//FACTORS
		if(s.fib.enabled)
		{
			f.writeTag("FACTORS"); f.newLine();
			f.writeTag("dx"); f.write(s.fib.dx); f.newLine();
			f.writeTag("dt"); f.write(s.fib.dt); f.newLine();
			f.writeTag("wound beg"); f.write(s.fib.woundBeg); f.newLine();
			f.writeTag("wound end"); f.write(s.fib.woundEnd); f.newLine();
			f.writeTag("wound val (T)"); f.write(s.fib.woundVal); f.newLine();
			f.writeTag("initVal (Fg)"); f.write(s.fib.initFg); f.newLine();
			f.writeTag("alpha1 (Tf)"); f.write(s.fib.alpha1); f.newLine();
			f.writeTag("alpha2 (Tm)"); f.write(s.fib.alpha2); f.newLine();
			f.writeTag("alpha3 (T)"); f.write(s.fib.alpha3); f.newLine();
			f.writeTag("alpha4 (Fg)"); f.write(s.fib.alpha4); f.newLine();
			f.writeTag("gamma1  (Tf)"); f.write(s.fib.gamma1); f.newLine();
			f.writeTag("gamma2  (Tm)"); f.write(s.fib.gamma2); f.newLine();
			f.writeTag("gamma3  (T)"); f.write(s.fib.gamma3); f.newLine();
			f.writeTag("beta1  "); f.write(s.fib.beta1); f.newLine();
			f.writeTag("beta2  "); f.write(s.fib.beta2); f.newLine();
			f.writeTag("beta3  "); f.write(s.fib.beta3); f.newLine();
			f.writeTag("k1 (T)"); f.write(s.fib.k1); f.newLine();
			f.writeTag("k3 (Fg, Fp)"); f.write(s.fib.k3); f.newLine();
			f.writeTag("C0 (T)"); f.write(s.fib.C0); f.newLine();
			f.writeTag("T0 (T)"); f.write(s.fib.T0); f.newLine(); f.newLine();
		}

		f.close();
	}
}

#endif