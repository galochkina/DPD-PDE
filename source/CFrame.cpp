#include "stdafx.h"

Calculus *Calc = NULL;

//=====================================================================
BEGIN_MESSAGE_MAP(Calculus, CMyFrame)
	ON_WM_PAINT()
	ON_WM_SETCURSOR()
	ON_WM_KEYDOWN()
	ON_WM_LBUTTONDBLCLK()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_COMMAND(ID_FILE_DPD, OnFileNew)
	ON_COMMAND(ID_TOOLS_UNITS, OnTrans)
	ON_UPDATE_COMMAND_UI(ID_VIEW_RESET, OnUpdateViewReset)
	ON_COMMAND(ID_TOOLS_OPTIONS,		OnToolsOptions)
	ON_WM_KEYDOWN()
END_MESSAGE_MAP()

//=====================================================================
Calculus::Calculus() 
{    
	pDialog = NULL;
	pDialog2 = NULL;
	pDialogOpt = NULL;
	pDialogProg = NULL;

	s = Setup();
	IO::readProps(s);

	//application properties
	s.opt.textSize		= 160;
	s.opt.textSize		= 30;
	s.opt.proc			= 2;	
	//s.opt.outputPath	= "C:\Projets\bloodCoagulation\dpd\newVersionAlen\Sim04";
	IO::readProps(s);

	if(!FS::dirExists(s.opt.outputPath)) FS::dirCreate(s.opt.outputPath);

	//DEFAULT VALUES
	//SETUP RUN	
	s.run.dt			= 0.001;		//s.run.dt
	s.run.dt2			= 0.001;		//s.run.dt
	s.run.stop			= 1000;			//stop time

	s.run.showStep		= 100;			//show every ...th s.run.step
	s.run.saveStep		= 100;			//save every ...th s.run.step
	s.run.txtStep		= 100000;		//print every ...th s.run.step
	s.run.bmpStep		= 100;			//bmp every ...th s.run.step
	s.run.excT			= 0.1;

	s.run.newInstance	= true;
	s.run.step			= 0;			//current s.run.step
	s.run.t				= 0;			//current time
	
	s.run.preSimulation	= false;
	s.run.preAnimOver	= false;
	s.run.preAnimTime	= 0;

	s.run.name			= "sim";	

	//SETUP SPACE
	s.spc.size			= VecF2(30, 10);				//s.spc.size of real space
	s.spc.prepArea		= false;
	s.spc.maxR			= 0;						//maximum r	
	s.spc.sizePrepX		= 10;


	//SETUP METHOD
	s.met.G			= VecF2(10, 0);		//total external F
	s.met.maxV		= 0;				//maximum v
	s.met.gamma		= 4.5;				//F function gamma			(s.met.type == 1)
	s.met.sigma		= 3;				//F function sigma			(s.met.type == 1)
	s.met.A			= 18.75;			//F function coefficients	(s.met.type == 1)
	s.met.rc		= 1;//0.1;
	s.met.p			= 1;
	
	//SETUP ATOMS
	s.atom.n.x		= 90;
	s.atom.n.y		= 30;
	s.atom.r		= 0.1;
	s.atom.rShow	= 0.05;
	s.atom.m		= 1;
	s.atom.v_maxPF	= 0;

	//SETUP ANALYSIS
	s.anl.density	= 25;
	s.anl.anlStep	= 10000;
	s.anl.clotStep	= 200;

	//PLATELETS
	s.plt.enabled	= true;
	s.plt.percent	= 0.007;
	s.plt.F2weak	= 100000;
	s.plt.F2strong	= 1000000;
	s.plt.F2fib		= 10000000;
	s.plt.F2t		= 20;
	s.plt.dc		= 1;
	s.plt.dd		= 1.5;
	s.plt.constInflow = false;
	s.plt.critFib	= 0.8;


	//SETUP FIB
	s.fib.enabled	= true;
	s.fib.dx		= 0.1;
	s.fib.dt		= 0.01;
	s.fib.alpha1	= 1;
	s.fib.alpha2	= 1;
	s.fib.alpha3	= 1;
	s.fib.alpha4	= 1;
	s.fib.gamma1	= 0.1;
	s.fib.gamma2	= 0.1;
	s.fib.gamma3	= 0.1;
	s.fib.beta1		= 1;
	s.fib.beta2		= 0.1;
	s.fib.beta3		= 0.1;
	s.fib.woundBeg	= 30;
	s.fib.woundEnd	= 35;
	s.fib.woundVal	= 1;
	s.fib.k1		= 1;
	s.fib.k3		= 1;
	s.fib.C0		= 1;
	s.fib.T0		= 1;
	s.fib.initFg	= 1;

	s.uni = SetupUnits();

	//CALCULUS MEMBERS
	c	= Clock();
	bc	= BoundCond(s);	//boundary conditions
	ps  = PdeSys();

	draw = Draw();

	cleared = true;
}

//=====================================================================
void Calculus::init()
{
	cleared = false; 

	s.met.sigmaSqrtdt = s.met.sigma / sqrt(s.run.dt); 
	FS::setFS(s);
	FS::createFS();

	if(s.run.newInstance) anl = Anl(s);

	bc	= BoundCond(s);
	if(s.run.newInstance) 
	{
		setAtoms();
		s.run.imgCount = 0;
	}
	setMaxRadius();
	mesh = Mesh(s, bc);
	MU::init(s.uni.mTo - s.uni.mFrom, s.uni.gTo - s.uni.gFrom, s.uni.sTo - s.uni.sFrom);
	if(s.run.newInstance) 
	{
		if(s.fib.enabled) ps = PdeSys(s);
		else ps = PdeSys();
	}

	draw = Draw(s);
	draw.SetPhysicalCoordinates(this);

	save(true);
}

//=====================================================================
void Calculus::setAtoms()
{
	DBL kPF = 4 * s.atom.v_maxPF / pow(s.spc.size.y, 2);

	int	k = 0;
	int	plNum = 1. / s.plt.percent;
	bool plBool = false;

	VecF2 da = (s.atom.n.x > 0 && s.atom.n.y > 0) ? s.spc.size / VecF2(s.atom.n.x, s.atom.n.y) : VecF20;
	
	VecF2 pos = da * 0.5;
	for(DBL i=0.5; pos.x < s.spc.size.x; pos.x = da.x * (++i))
	{
		bool inPrepArea = pos.x <= s.spc.sizePrepX;
		pos.y = da.y * 0.5;
		for(DBL j=0.5; pos.y < s.spc.size.y; pos.y = da.y * (++j))
		{
			k++;
			plBool = (s.plt.enabled && k%plNum==0) ? true : false;

			VecF2 v = VecF20;
			v.x += kPF * pos.y * (s.spc.size.y - pos.y);

			if(s.plt.enabled && s.fib.enabled && pos.x >= s.fib.woundBeg && pos.x <= s.fib.woundEnd && pos.y < s.spc.size.y / s.atom.n.y)
			{
				al.addLater(Atom(Atom::TypeSPlt, pos, VecF20, inPrepArea, true, true));
			}
			else if(s.cplt.enabled && s.fib.enabled && pos.x >= s.fib.woundBeg && pos.x <= s.fib.woundEnd && pos.y < s.spc.size.y / s.atom.n.y)
			{
				//al.addLater(Atom(Atom::TypeCPlt, pos, VecF20, inPrepArea, true, true));
			}
			else
			{
				int type = plBool ? Atom::TypeSPlt : Atom::TypePlasma;
				al.addLater(Atom(type, pos, v, false, false, false));
			}
		}
	}

	al.update();
}

//=====================================================================
void Calculus::calculateBoxPair(BoxPair &bp)
{
	int n1 = bp.A->atoms.size();
	if(bp.singleBox)
	{
		for(int i=0; i<n1; i++)
			for(int j=i+1; j<n1; j++)
				Interact::atomAtom(s, bp.A->atoms[i], bp.A->atoms[j], bp.displacement, false, platPairs, ps);
	}
	else
	{
		int n2 = bp.B->atoms.size();
		for(int i=0; i<n1; i++)
			for(int j=0; j<n2; j++)
				Interact::atomAtom(s, bp.A->atoms[i], bp.B->atoms[j], bp.displacement, bp.oneWay, platPairs, ps);
	}
}

//=====================================================================
void Calculus::calcAtomF()
{
	
	int n = mesh.boxPairsList.size();
	int m = n / 2;
	#pragma omp parallel for
	for(int i=0; i<m; i++)
	{
		for(unsigned int j=0; j<mesh.boxPairsList[i].size(); j++)
			calculateBoxPair(mesh.boxPairsList[i][j]);	
	}

	#pragma omp parallel for
	for(int i=m; i<n; i++)
	{
		for(unsigned int j=0; j<mesh.boxPairsList[i].size(); j++)
			calculateBoxPair(mesh.boxPairsList[i][j]);	
	}
}

//=====================================================================
void Calculus::calcOtherF()
{
	//ADD PLATELET CONNECTION FORCES AND VERIFY CONNECTIONS
	DBL h0 = 2 * s.atom.r * s.plt.dd;
	int n = platPairs.l.size();
	//#pragma omp parallel for
	for(int i=0; i<n; i++)
	{		
		platPairs.l[i].t += s.run.dt2;
		Atom *a1 = &(al.obj[platPairs.l[i].a]);
		Atom *a2 = &(al.obj[platPairs.l[i].b]);
		if(!(a1->free || a2->free))
		{
			VecF2 r = a1->pos - a2->pos;
			DBL h = MODUL(r);
			if(h > h0) 
				platPairs.l[i].del = true;
			else if(h > s.atom.r * 2)
			{
				VecF2 tmpF = Interact::pltFnPlus(r, h, platPairs.l[i].d);

				if(a1->oldPlt && a2->oldPlt)
				{
					a1->FPL3 += tmpF;
					a2->FPL3 -= tmpF;
				}
				else if(platPairs.l[i].t < s.plt.F2t)
				{
					a1->FPL1 += tmpF;
					a2->FPL1 -= tmpF;
				}
				else
				{
					a1->FPL2 += tmpF;
					a2->FPL2 -= tmpF;
				}
			}
		}
		else
		{
			platPairs.l[i].del = true;
		}
	}

	platPairs.refresh();
}

//=====================================================================
void Calculus::setOldPlat()
{
	int n = al.obj.size();
	#pragma omp parallel for
	for(int i=0; i<n; i++)
	{
		if(!al.obj[i].free && ps.isOld(al.obj[i].pos)) al.obj[i].oldPlt = true;
	}
}

//=====================================================================
void Calculus::moveAtoms()
{
	int n = al.obj.size();
	#pragma omp parallel for
	for(int i=0; i<n; i++)
	{
		if(!al.obj[i].free) 
			Move::move(s, al.obj[i], al, false, s.run.dt);
	}
	al.update();
}

//=====================================================================
void Calculus::nextStep()
{
	int k = s.run.dt / s.run.dt2;

	mesh.setParticles(al.obj);	
	calcAtomF();
	
	Interact::selfMirror(s, al, ps);
	if((s.plt.enabled) && s.fib.enabled) setOldPlat();	

	for(int j = 0; j<k ; j++)
	{
		calcOtherF();
		
		int n = al.obj.size();
		#pragma omp parallel for
		for(int i=0; i<n; i++)
			if(!al.obj[i].free) Move::move(s, al.obj[i], al, true, s.run.dt2);
		al.update();
	}

	int n = al.obj.size();
	#pragma omp parallel for
	for(int i=0; i<n; i++)
		if(!al.obj[i].free) Move::move(s, al.obj[i], al, false, s.run.dt);
	al.update();

	anl.update(s, al);

	n = al.obj.size();
	#pragma omp parallel for
	for(int i=0; i<n; i++) 
	{
		al.obj[i].resetF();
	}
}//

//=====================================================================
void Calculus::setMaxRadius()
{
	s.spc.maxR = s.met.rc;
}

//=====================================================================
void Calculus::OnPaint()
{ 
	CPaintDC dc(this);

	CBitmap* pBitmap;   
	CRect rect;    
	GetClientRect(rect);   
	if(MyDc) delete MyDc;    
	MyDc = new CDC;    
	MyDc->CreateCompatibleDC(&dc);    
	CBitmap bitmap; bitmap.CreateCompatibleBitmap(&dc, rect.right, rect.bottom);    
	pBitmap = MyDc->SelectObject(&bitmap);   
	CBrush  BrushNew (RGB(255,255,255));    
	CBrush* BrushOld = MyDc->SelectObject(&BrushNew);   
	MyDc->FillRect(rect, &BrushNew);    
	MyDc->SelectObject(BrushOld);
	MyDc->SaveDC(); 
	CPen  PenNew (PS_SOLID, 2, RGB(0,0,255));
	MyDc->SelectObject(&PenNew); 

	draw.drawDPD(this, al, platPairs, ps.pdeFp, anl, mesh);

	dc.BitBlt(0, 0, rect.right, rect.bottom, MyDc, 0, 0, SRCCOPY); 
	MyDc->SelectObject(pBitmap);
	if(MyDc) delete MyDc;
	MyDc = new CClientDC(this); 
} // 

//=====================================================================
void Calculus::Run()
{ 
#if _DEBUG
	RunBodyDPD();
#else
	AfxBeginThread(ThreadProc, this);
#endif	
}

//=====================================================================
void Calculus::RunBodyDPD()
{	
	string fRes = "";

	if(s.run.newInstance)
	{
		clear();

		if(s.run.preSimulation) s.run.preAnimOver = false;

		if(s.fib.enabled) 
		{
			s.fib.time = 0;
			s.run.pdeT = 0;
		}

		s.run.step = 0;
		s.run.t = 0;

		init();

		IO::txtParams(s);

		c = Clock();
		c.start();
	}
	else
	{	
	}

	ostringstream ss;
	s.run.newInstance = false;

	bool closedWound = true;
	while(s.run.t<s.run.stop)
	{
		s.run.step++;
		//s.run.t+=s.run.dt;
		s.run.t = s.run.dt * s.run.step;
		//data exchange
		
		if(s.fib.enabled)
		{
			if(s.run.pdeT + s.run.excT <= s.run.t)
			{
				int st = 7 + s.run.excT / s.fib.dt;
				progressBeg("Calculating PDE...", st);
				anl.anlPdeVD(s, ps, pDialogProg->bar, pDialogProg->desc);
				
				pDialogProg->desc.SetWindowTextA("ADI method.");
				//bool closedWound = anl.plClot.getLast() > 20;
				closedWound = s.run.t * MU::time > 5;
				while(s.fib.time <= s.run.pdeT + s.run.excT)
				{
					ps.nextStep(s, closedWound);
					s.fib.time += s.fib.dt;
					pDialogProg->bar.StepIt();
				}
				s.run.pdeT += s.run.excT;

				pDialogProg->desc.SetWindowTextA("Making graphs.");
				ps.draw(s, "fib");
				pDialogProg->bar.StepIt();
				pDialogProg->desc.SetWindowTextA("Writing to txt files.");
				//ps->txt(s.run.step, "fib");
				pDialogProg->bar.StepIt();
				progressEnd();
			}
		}
		
		

		nextStep();

		c.addCounter();

		
		if(s.run.step % s.run.saveStep == 0) save(false);
		//if(s.run.step % s.run.txtStep == 0) ;

		//calls paint method
		if(s.run.step % s.run.showStep == 0)
		{
			c.stop();
			MySetWindowText();
			Invalidate(FALSE); 
			SendMessage(WM_PAINT);
			c.start();
		}

		if(s.run.step % s.run.bmpStep == 0) 
		{
			DrawImg img(s);
			Gdiplus::GdiplusStartupInput gdiplusStartupInput;
			ULONG_PTR gdiplusToken;
			Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
			string name = s.run.name + "_img";
			img.drawSimpleImage(FS::getFN(FS::image, name, "png", s.run.step), s, al, platPairs, ps);
			Gdiplus::GdiplusShutdown(gdiplusToken);
			s.run.imgCount++;
		}
		 
		if(STOPKEY) { break; }
	}

	MySetWindowText();
	Invalidate(FALSE); 
	SendMessage(WM_PAINT);
	anl.txtHist(s);
	anl.drawHist(s);

} // 

//=====================================================================
void Calculus::save(bool justStp, string filePath)
{
	if(justStp)
	{
		if(filePath == "")
			filePath = FS::getFN(FS::save, s.run.name, "stp");

		IO::saveDPDStp(s, filePath);
	}
	else
	{
		if(filePath == "")
			filePath = FS::getFN(FS::save, s.run.name, "s04", s.run.step);
		
		Calc->progressBeg("Saving...", 6);
		IO::saveDPD(s, filePath, al, platPairs, c, ps, anl, Calc->pDialogProg->bar, Calc->pDialogProg->desc);
		Calc->progressEnd();
	}
}
//=====================================================================
bool Calculus::load(string fileName)
{
	clear();
	Calc->progressBeg("Loading...", 6);
	if(IO::loadDPD(s, fileName, al, platPairs, c, ps, anl, Calc->pDialogProg->bar, Calc->pDialogProg->desc))
	{		
		init();
		MySetWindowText();
		Invalidate(FALSE); 
		SendMessage(WM_PAINT);
		Calc->progressEnd();
		return true;		
	}
	Calc->progressEnd();
	
	return false;
}

//=====================================================================
void Calculus::clear()
{
	if(cleared) return;
	free(al.obj);
	free(platPairs.l);
	mesh.destroy();
	anl.clear();
	ps.free();
	cleared = true;
}

//=====================================================================
void Calculus::OnToolsOptions()      
{ 
	if(!pDialogOpt) 
	{
		pDialogOpt = new COptions(); 
		pDialogOpt->Create(COptions::IDD);
	}
	pDialogOpt->cf1 = s.opt.cf1;
	pDialogOpt->cf2 = s.opt.cf2;
	pDialogOpt->cfib = s.opt.cfib;
	pDialogOpt->cpla = s.opt.cpla;
	pDialogOpt->cplt = s.opt.cplt;

	pDialogOpt->ShowWindow(SW_SHOW);
	pDialogOpt->UpdateData(false);
} //

//====================================================================
void Calculus::OnFileNew()      
{ 
	if(!pDialog)
	{
		pDialog = new CBeginNew(); 
		pDialog->Create(CBeginNew::IDD);
	}
	pDialog->ShowWindow(SW_SHOW);
	pDialog->UpdateData(false);
} //

//====================================================================
void Calculus::OnTrans()      
{ 
	if(!pDialog2) 
	{
		pDialog2 = new CTrans(); 
		pDialog2->Create(CTrans::IDD);
	}
	pDialog2->ShowWindow(SW_SHOW);
	pDialog2->UpdateData(false);
} //
//====================================================================
void Calculus::progressBeg(string title, int totalSteps)      
{ 
	if(!pDialogProg) 
	{
		pDialogProg = new CProgress(); 
		pDialogProg->Create(CProgress::IDD);
	}
	pDialogProg->bar.SetRange(0, totalSteps);
	pDialogProg->bar.SetStep(1);
	pDialogProg->bar.SetPos(0);
	pDialogProg->desc.SetWindowTextA("");
	pDialogProg->SetWindowTextA(title.c_str());
	pDialogProg->ShowWindow(SW_SHOW);
	pDialogProg->UpdateData(false);
} 

void Calculus::progressEnd()      
{ 
	pDialogProg->ShowWindow(SW_HIDE);
} 

//====================================================================
void Calculus::setupHotChange()
{
	ps.paramHotChange(s);
}