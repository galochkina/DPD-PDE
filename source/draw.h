#ifndef DRAW_H
#define	DRAW_H

class Draw
{
public:
	Draw()
	{
		text = false;
		stp = NULL;
	}

	Draw(Setup &_stp)
	{
		text = false;
		stp = &_stp;
	}

	void SetPhysicalCoordinates(CMyFrame *fr)
	{
		VecF2 dr = VecF2(1,1);

		dr = stp->spc.size;
		if(stp->fib.enabled) dr.y*=2.05;
		dr *= 1.1;
		VecF2 rc = VecF2(dr.x / 2, dr.y / 2);
		fr->Getaabxby (rc, dr);
				
		CRect rect;
		fr->GetClientRect(rect);
		int I = rect.right;
		int J = rect.bottom;

		d = 30;
		textBeg = VecL2(I*4/5, d);
		textRow = 0;
		text = true;
		int W = (DBL)I*4/5 - 2*d;
		int H = (DBL)J - 2*d;

		VecF2 A = dr;
		if(A.y * W <= A.x * H)
		{
			DBL k = ((DBL) W) / A.x;
			L = VecF2(k, -k);
			beg = VecL2(d, d + ((DBL)H - A.y*L.y)/2);
		}
		else
		{
			DBL k = ((DBL) H) / A.y;
			L = VecF2(k, -k);
			beg = VecL2(d, d + H);
		}
	}
	
	void drawDPD(CMyFrame *fr, AtomList &al, PPL &platPairs, Pde3 &fib, Anl &anl, Mesh &mesh)
	{
		if(stp == NULL) return;

		CDC *MyDc = fr->MyDc;

		MyDc->SaveDC();

		drawOneDPD(fr, al, platPairs, VecF20, true, true, false, fib);

		if(stp->fib.enabled)
			drawOneDPD(fr, al, platPairs, VecF2(0, 1.05*stp->spc.size.y), false, true, true, fib);

		drawText(MyDc, anl);

		MyDc->RestoreDC(-1); 

	}

private:
	int d;
	bool text;
	VecL2 textBeg;
	int textRow;
	VecF2 L;
	VecL2 beg;
	Setup *stp;

	inline VecL2 getTextC() { return VecL2(textBeg.x, textBeg.y + stp->opt.textSpace * textRow); }

	inline CPoint getC(VecF2 point)
	{
		VecL2 p = beg + VecL2(point.x * L.x, point.y * L.y);
		return CPoint(p.x, p.y);
	}

	void textWrite(CDC *MyDc, string a)
	{
		VecL2 p = getTextC();
		MyDc->TextOutA(p.x, p.y, (LPCTSTR)CString(a.c_str()));
		textRow++;
	}

	void textWrite(CDC *MyDc, string desc, string val, string unit="", string val2="")
	{
		ostringstream ss;
		ss << left << setw(9) << desc << flush;
		ss << scientific << setprecision(3) << val << flush;
		if(unit != "") ss << " " << left << unit << flush;
		if(val2 != "") ss << "  (" << scientific << setprecision(3) << val2 << ")" << flush;
		textWrite(MyDc, ss.str().c_str());
	}

	void textWrite(CDC *MyDc, string desc, DBL val, string unit, DBL val2)
	{
		bool sci = val >= 1000 && val < 0.001;
		bool sci2 = val2 >= 1000 && val2 < 0.001;
		textWrite(MyDc, desc, toString(val, 3, sci), unit, toString(val2, 3, sci2));
	}

	void textWrite(CDC *MyDc, string desc, DBL val, string unit)
	{
		bool sci = val >= 1000 && val < 0.001;
		textWrite(MyDc, desc, toString(val, 3, sci), unit);
	}

	void textWrite(CDC *MyDc, string desc, DBL val)
	{
		bool sci = val >= 1000 && val < 0.001;
		textWrite(MyDc, desc, toString(val, 3, sci));
	}

	void textWrite(CDC *MyDc, string desc, DBL val, DBL val2)
	{
		bool sci = val >= 1000 && val < 0.001;
		bool sci2 = val2 >= 1000 && val2 < 0.001;
		textWrite(MyDc, desc, toString(val, 3, sci), "", toString(val2, 3, sci2));
	}

	void textNewRow()
	{
		textRow++;
	}

	void textReset()
	{
		textRow = 0;
	}

	void drawLine(CDC *MyDc, VecF2 p1, VecF2 p2, COLORREF c, int penSize)
	{
		CPen pen (PS_SOLID, penSize, c);
		MyDc->SelectObject(&pen);

		MyDc->MoveTo(getC(p1));
		MyDc->LineTo(getC(p2));
	}

	void drawAtom(CDC *MyDc, Atom &p, VecF2 tr, bool drPlU)
	{
		MyDc->SaveDC();
		VecF2 dR = (p.isSPlt()) ? stp->atom.r * VecF21 : stp->atom.rShow * VecF21;
		VecF2 R1 = tr + p.pos - dR;
		VecF2 R2 = tr + p.pos + dR;
		CRect rect(getC(R1),getC(R2));

		//CPen pen (PS_SOLID, 0, RGB(0,0,0));
		
		COLORREF c = RGB(0,0,0);
		if(p.isPlasma())
		{
			if(p.prep) c = RGB(155/2, 204/2, 204/2);
			else c = RGB(155, 204, 204);
		}
		else if(p.isSPlt())
		{
			if(drPlU)
			{
				if(p.oldPlt) c = RGB(0,100,0);
				else c = RGB(0,200,0);
			}
			else
				c = RGB(0,0,0);
		}
		else 
		{
			c = RGB(100, 204, 204);
		}

		CPen pen (PS_SOLID, 0, c);
		MyDc->SelectObject(&pen);

		CBrush brush(c);
		MyDc->SelectObject(&brush);

		MyDc->Ellipse(rect);

		MyDc->RestoreDC(-1); 
	}

	void drawConn(CDC *MyDc, PPL &plConns, AtomList &al, VecF2 tr)
	{
		for(vector<PP>::iterator iter = plConns.l.begin(); iter != plConns.l.end(); iter++)
		{
			int red = (iter->t >= stp->plt.F2t) ? 150 : 255;
			drawLine(MyDc, tr+al.obj[iter->a].pos, tr+al.obj[iter->b].pos, RGB(red, 0, 0), 2);
		}
	} 

	void drawGrid(CDC *MyDc, DBL dx, DBL dy)
	{
		for(DBL i = 0; i <= stp->spc.size.x; i+=dx)
			drawLine(MyDc, VecF2(i, 0), VecF2(i, stp->spc.size.y), RGB(0, 0, 0), 0);

		for(DBL i = 0; i <= stp->spc.size.y; i+=dy)
			drawLine(MyDc, VecF2(0, i), VecF2(stp->spc.size.x, i), RGB(0, 0, 0), 0);
	} 

	void drawMesh(CDC *MyDc, Mesh &mesh)
	{
		for(DBL x = 0; x <= stp->spc.sizePrepX; x+=mesh.getIntXP())
			drawLine(MyDc, VecF2(x, 0), VecF2(x, stp->spc.size.y), RGB(0, 0, 0), 2);

		for(DBL x = stp->spc.sizePrepX; x <= stp->spc.size.x; x+=mesh.getIntXS())
			drawLine(MyDc, VecF2(x, 0), VecF2(x, stp->spc.size.y), RGB(0, 0, 0), 2);

		for(DBL y = 0; y <= stp->spc.size.y; y+=mesh.getIntY())
			drawLine(MyDc, VecF2(0, y), VecF2(stp->spc.size.x, y), RGB(0, 0, 0), 2);
	} 

	void drawCurve(CDC *MyDc, MX &X, DBL dx)
	{
		for(int i = 0; i < X.dimR() - 1; i++)
		{
			DBL x = i * dx;
			drawLine(MyDc, VecF2(x, X(i,0)), VecF2(x+(dx), X(i+1,0)), RGB(0, 0, 0), 2);
		}
	}

	void drawBoxPairs(CDC *MyDc, Mesh &mesh)
	{
		DBL intX = mesh.getIntXS();
		DBL intY = mesh.getIntY();

		for(unsigned int i=0; i<mesh.boxPairsList.size(); i++)
		{
			for(unsigned int j=0; j<mesh.boxPairsList[i].size(); j++)
			{
				VecL2 posA = mesh.boxPairsList[i][j].A->pos;
				VecF2 a((0.5 + posA.x) * intX, (0.5 + posA.y) * intY);
				VecL2 posB = mesh.boxPairsList[i][j].B->pos;
				VecF2 b((0.5 + posB.x) * intX, (0.5 + posB.y) * intY);
				drawLine(MyDc, a, b, RGB(255, 0, 0), 2);
			}
		}
	}

	void drawText(CDC *MyDc, Anl &anl)
	{	
		CFont font;
		font.CreatePointFont(stp->opt.textSize, "Consolas", MyDc);
		MyDc->SelectObject(&font);

		if(text)
		{
			textReset();

			DBL cT = stp->uni.getTimeM();
			DBL cV = stp->uni.getVelocityM();
			DBL cMu = stp->uni.getViscM();
			DBL cL = stp->uni.getLengthM();
			DBL cM = stp->uni.getMassM();
			DBL cRho = stp->uni.getDensityM();
			DBL cAcc = stp->uni.getAccM();

			//non const
			textWrite(MyDc, "MEASURED");
			textWrite(MyDc, "time =", stp->run.t * cT, "s", stp->run.t);
			textWrite(MyDc, "energy =", anl.energy.getLast());
			textWrite(MyDc, "vMax =", anl.m.vMax, "m/s", anl.m.vMax / MU::vel);
			textWrite(MyDc, "vAvg =", anl.m.vAvg.x, "m/s", anl.m.vAvg.x / MU::vel);
			textWrite(MyDc, "vM/vA =", anl.m.vMax / anl.m.vAvg.x);
			textWrite(MyDc, "visc =", anl.m.visc, "kg/ms", anl.m.visc / MU::visc);
			if(stp->plt.enabled) 
			{
				textWrite(MyDc, "clot =", anl.plClot.getLast());
				textWrite(MyDc, "core =", anl.plCore.getLast());
			}
			textNewRow();

			//const
			textWrite(MyDc, "PHYSICAL");
			textWrite(MyDc, "D =", stp->spc.size.y * cL, "m", stp->spc.size.y);
			textWrite(MyDc, "L =", stp->spc.size.x * cL, "m", stp->spc.size.x);
			textWrite(MyDc, "rho=", anl.m.density, "kg/m3", anl.m.density / MU::dens);
			textWrite(MyDc, "G =", stp->met.G.x * cAcc, "m/s2", stp->met.G.x);
			textWrite(MyDc, "Nx =", stp->atom.n.x);
			textWrite(MyDc, "Ny =", stp->atom.n.y);
			textWrite(MyDc, "r =", stp->atom.r * cL, "m", stp->atom.r);
			textWrite(MyDc, "m =", stp->atom.m * cM, "kg", stp->atom.m);
			textNewRow();

			textWrite(MyDc, "DPD");
			textWrite(MyDc, "gamma =", stp->met.gamma);
			textWrite(MyDc, "sigma =", stp->met.sigma);
			textWrite(MyDc, "rc =", stp->met.rc);
			textWrite(MyDc, "k =", stp->met.p);
			textWrite(MyDc, "a =", stp->met.A);
			textWrite(MyDc, "dt =", stp->run.dt);
			textWrite(MyDc, "tau =", stp->run.pdeT);
			textNewRow();

			if(stp->plt.enabled)
			{
				textWrite(MyDc, "PLATELET");				
				textWrite(MyDc, "Fw =", stp->plt.F2weak);
				textWrite(MyDc, "Fs =", stp->plt.F2strong);
				textWrite(MyDc, "Ffib =", stp->plt.F2fib);
				textWrite(MyDc, "time =", stp->plt.F2t);
				textWrite(MyDc, "u crit =", stp->plt.critFib);
				textWrite(MyDc, "dc =", stp->plt.dc);
				textWrite(MyDc, "dd =", stp->plt.dd);
				textWrite(MyDc, "plts =", stp->plt.percent, "%");
				textNewRow();
			}

			if(stp->fib.enabled)
			{
				textWrite(MyDc, "PDEs");	
				textWrite(MyDc, "dx =", stp->fib.dx);
				textWrite(MyDc, "dt =", stp->fib.dt);
				textWrite(MyDc, "alpha1 =", stp->fib.alpha1);
				textWrite(MyDc, "alpha2 =", stp->fib.alpha2);
				textWrite(MyDc, "k1 =", stp->fib.k1);
				textWrite(MyDc, "k3 =", stp->fib.k3);
				textWrite(MyDc, "gamma3 =", stp->fib.gamma3);
				textWrite(MyDc, "C0 =", stp->fib.C0);
				textWrite(MyDc, "T0 =", stp->fib.T0);
				textWrite(MyDc, "init Fg =", stp->fib.initFg);
				textNewRow();
			}
		}

		font.DeleteObject();
	}

	void drawOneDPD(CMyFrame *fr, AtomList &al, PPL &platPairs, VecF2 tr, bool drPlasma, bool drPlConn, bool drPlU, Pde3 &fib)
	{
		//if(stp == NULL || stats == NULL) return;
		if(stp == NULL) return;

		CDC *MyDc = fr->MyDc;

		MyDc->SaveDC();

		drawLine(MyDc, VecF2(tr.x,tr.y), VecF2(tr.x+stp->spc.size.x, tr.y), RGB(0,0,0), 2);
		drawLine(MyDc, VecF2(tr.x, tr.y+stp->spc.size.y), VecF2(tr.x+stp->spc.size.x, tr.y+stp->spc.size.y), RGB(0,0,0), 2);
		
		//iterates over list of particles and draws them
		for(LA::iterator iter = al.obj.begin(); iter != al.obj.end(); iter++)
			if(!iter->free) 
				if(drPlasma || !iter->isPlasma())
					drawAtom(MyDc, *iter, tr, drPlU);

		if(stp->plt.enabled && drPlConn) drawConn(MyDc, platPairs, al, tr);

		if(stp->spc.prepArea)
		{
			VecF2 A = VecF2(tr.x+stp->spc.sizePrepX, tr.y);
			VecF2 B = VecF2(tr.x+stp->spc.sizePrepX, tr.y+stp->spc.size.y);

			MyDc->MoveTo(getC(A));
			MyDc->LineTo(getC(B));
		}

		MyDc->RestoreDC(-1); 
	}
};

#endif