#include "stdafx.h"
#include "dialog.h"

extern Calculus*	Calc;

//======================================================================
BEGIN_MESSAGE_MAP(CBeginNew, CDialog)
	ON_BN_CLICKED(IDC_START,	OnBnClickedButton1)
	ON_BN_CLICKED(IDC_LOAD,		OnBnLoad)
	ON_BN_CLICKED(IDC_SAVE,		OnBnClickedButton3)
	ON_BN_CLICKED(IDC_RESUME,	OnBnClickedButton5)
	ON_BN_CLICKED(IDC_INIT,		OnBnClickedInit)
	ON_BN_CLICKED(IDC_PL_ENABLED, &CBeginNew::OnBnClickedPlEnabled)
	ON_BN_CLICKED(IDC_FAC_ENABLED, &CBeginNew::OnBnClickedFcEnabled)
END_MESSAGE_MAP()

//======================================================================
void CBeginNew::DoDataExchange(CDataExchange* pDX)
{
	//SETUP SPACE
	DDX_Text(pDX, IDC_SIZEX,			Calc->s.spc.size.x);
	DDX_Text(pDX, IDC_SIZEY,			Calc->s.spc.size.y);
	DDX_Text(pDX, IDC_PREP_AREA_X,		Calc->s.spc.sizePrepX);

	int tmp = boolToInt(Calc->s.spc.prepArea);
	DDX_Check(pDX, IDC_PREP_AREA, tmp);
	Calc->s.spc.prepArea = intToBool(tmp); 


	//****************************************************************************************************************
	//SETUP METHOD
	DDX_Text(pDX, IDC_EXT_X,			Calc->s.met.G.x);
	DDX_Text(pDX, IDC_EXT_Y,			Calc->s.met.G.y);
	DDX_Text(pDX, IDC_MAX_V,			Calc->s.met.maxV);

	DDX_Text(pDX, IDC_DPD_GAMMA,		Calc->s.met.gamma);
	DDX_Text(pDX, IDC_DPD_SIGMA,		Calc->s.met.sigma);
	DDX_Text(pDX, IDC_DPD_RC,			Calc->s.met.rc);
	DDX_Text(pDX, IDC_DPD_P,			Calc->s.met.p);
	DDX_Text(pDX, IDC_DPD_A,			Calc->s.met.A);

	//****************************************************************************************************************
	//SETUP ATOMS
	DDX_Text(pDX, IDC_COUNT_X,			Calc->s.atom.n.x);
	DDX_Text(pDX, IDC_COUNT_Y,			Calc->s.atom.n.y);
	DDX_Text(pDX, IDC_RADIUS,			Calc->s.atom.r);
	DDX_Text(pDX, IDC_RADIUS_SHOW,		Calc->s.atom.rShow);
	DDX_Text(pDX, IDC_MASS,				Calc->s.atom.m);
	DDX_Text(pDX, IDC_V_MAX_PF,			Calc->s.atom.v_maxPF);


	//****************************************************************************************************************
	//SETUP ANALYSIS
	DDX_Text(pDX, IDC_A_DENSITY,		Calc->s.anl.density);
	DDX_Text(pDX, IDC_A_STEP_ANL,		Calc->s.anl.anlStep);
	DDX_Text(pDX, IDC_A_STEP_CLOT,		Calc->s.anl.clotStep);		


	//****************************************************************************************************************
	//SETUP RUN
	DDX_Text(pDX, IDC_DT,				Calc->s.run.dt);
	DDX_Text(pDX, IDC_STOP,				Calc->s.run.stop);
	DDX_Text(pDX, IDC_SHOW,				Calc->s.run.showStep);
	DDX_Text(pDX, IDC_PRINT,			Calc->s.run.txtStep);
	DDX_Text(pDX, IDC_SAVE_STEP,		Calc->s.run.saveStep);
	DDX_Text(pDX, IDC_BMP_STEP,			Calc->s.run.bmpStep);
	DDX_Text(pDX, IDC_T,				Calc->s.run.t);
	DDX_Text(pDX, IDC_EXC,				Calc->s.run.excT);
	DDX_Text(pDX, IDC_DT2,				Calc->s.run.dt2);	

	CString csn (Calc->s.run.name.c_str());
	DDX_Text(pDX, IDC_NAME,				csn);
	Calc->s.run.name = string((LPCTSTR)csn);

	//****************************************************************************************************************
	//SETUP PLATELETS

	tmp = boolToInt(Calc->s.plt.enabled);
	DDX_Check(pDX, IDC_PL_ENABLED, tmp);
	Calc->s.plt.enabled = intToBool(tmp); 

	tmp = boolToInt(Calc->s.plt.constInflow);
	DDX_Check(pDX, IDC_PL_INFLOW, tmp);
	Calc->s.plt.constInflow = intToBool(tmp); 

	DDX_Text(pDX, IDC_PL_N,			Calc->s.plt.percent);
	DDX_Text(pDX, IDC_PL_FWEAK,		Calc->s.plt.F2weak);	
	DDX_Text(pDX, IDC_PL_FSTRONG,	Calc->s.plt.F2strong);
	DDX_Text(pDX, IDC_PL_FFIB,		Calc->s.plt.F2fib);
	DDX_Text(pDX, IDC_PL_FT,		Calc->s.plt.F2t);
	DDX_Text(pDX, IDC_PL_DC,		Calc->s.plt.dc);
	DDX_Text(pDX, IDC_PL_DD,		Calc->s.plt.dd);
	DDX_Text(pDX, IDC_PL_U,			Calc->s.plt.critFib);

	//****************************************************************************************************************
	//SETUP FACTORS
	tmp = boolToInt(Calc->s.fib.enabled);
	DDX_Check(pDX, IDC_FAC_ENABLED, tmp);
	Calc->s.fib.enabled = intToBool(tmp);

	DDX_Text(pDX, IDC_FAC_DX,		Calc->s.fib.dx);
	DDX_Text(pDX, IDC_FAC_DT,		Calc->s.fib.dt);
	DDX_Text(pDX, IDC_FAC_ALPHA1,	Calc->s.fib.alpha1);
	DDX_Text(pDX, IDC_FAC_ALPHA2,	Calc->s.fib.alpha2);
	DDX_Text(pDX, IDC_FAC_ALPHA3,	Calc->s.fib.alpha1);
	DDX_Text(pDX, IDC_FAC_ALPHA4,	Calc->s.fib.alpha2);
	DDX_Text(pDX, IDC_FAC_BETA1,	Calc->s.fib.beta1);
	DDX_Text(pDX, IDC_FAC_BETA2,	Calc->s.fib.beta2);
	DDX_Text(pDX, IDC_FAC_BETA3,	Calc->s.fib.beta3);
	DDX_Text(pDX, IDC_FAC_GAMMA1,	Calc->s.fib.gamma1);
	DDX_Text(pDX, IDC_FAC_GAMMA2,	Calc->s.fib.gamma2);
	DDX_Text(pDX, IDC_FAC_GAMMA3,	Calc->s.fib.gamma3);
		

	DDX_Text(pDX, IDC_FAC_W_BEG,	Calc->s.fib.woundBeg);
	DDX_Text(pDX, IDC_FAC_W_END,	Calc->s.fib.woundEnd);
	DDX_Text(pDX, IDC_FAC_W_VAL,	Calc->s.fib.woundVal);
	DDX_Text(pDX, IDC_FAC_K1,		Calc->s.fib.k1);
	DDX_Text(pDX, IDC_FAC_K3,		Calc->s.fib.k3);
	DDX_Text(pDX, IDC_FAC_C0,		Calc->s.fib.C0);
	DDX_Text(pDX, IDC_FAC_T0,		Calc->s.fib.T0);
	DDX_Text(pDX, IDC_FAC_INIT_FG,	Calc->s.fib.initFg);

	switchPl(Calc->s.plt.enabled);
	switchFc(Calc->s.fib.enabled);
} //

//======================================================================
//Start
void CBeginNew::OnBnClickedButton1()
{  	
	Calc->s.run.newInstance = true;
	CDialog::OnOK();
	Calc->Run(); 
}  //

//======================================================================
//Init
void CBeginNew::OnBnClickedInit()
{  	
	UpdateData(true);
	Calc->s.run.newInstance = true;
	Calc->init();
	Calc->s.run.newInstance = false;
}  //

//======================================================================
//Load
void CBeginNew::OnBnLoad()
{  	
	CFileDialog f(true); 
	if(f.DoModal() == IDOK)
	{
		if(Calc->load((string)f.GetPathName()))
			UpdateData(false);
	}
}  //

//======================================================================
//Save
void CBeginNew::OnBnClickedButton3()
{  	
	UpdateData(true);
	CFileDialog f(false); 
	if(f.DoModal() == IDOK)
		Calc->save(false, (string)f.GetPathName());
}  //

//======================================================================
//SaveSetup
void CBeginNew::OnBnClickedSavestp()
{
	UpdateData(true);
	CFileDialog f(false); 
	if(f.DoModal() == IDOK)
		Calc->save(true, (string)f.GetPathName());
}

//======================================================================
//Resume
void CBeginNew::OnBnClickedButton5()
{  	 
	CDialog::OnOK();
	Calc->bc = BoundCond(Calc->s);
	Calc->setupHotChange();
	Calc->Run(); 
}  //

//======================================================================
//Img
void CBeginNew::OnBnClickedImg()
{
	CFileDialog f(false); 
	if(f.DoModal() == IDOK)
	{
		string path = f.GetPathName();
		DrawImg img(Calc->s);
		Gdiplus::GdiplusStartupInput gdiplusStartupInput;
		ULONG_PTR gdiplusToken;
		Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
		img.drawSimpleImage(path, Calc->s, Calc->al, Calc->platPairs, Calc->ps);
		Gdiplus::GdiplusShutdown(gdiplusToken);
	}
}

//======================================================================
//Txt Params
void CBeginNew::OnBnClickedTxtParams()
{
	UpdateData(true);
	IO::txtParams(Calc->s, Calc->s.opt.outputPath + "\\");
}

//======================================================================
//Events
void CBeginNew::OnBnClickedPlEnabled()
{
	CButton *pEnabled	= (CButton *) GetDlgItem(IDC_PL_ENABLED);
	switchPl(pEnabled->GetCheck() == 1);
}

void CBeginNew::switchPl(bool enable)
{
	CButton *pN			= (CButton *) GetDlgItem(IDC_PL_N);
	CButton *pFweak		= (CButton *) GetDlgItem(IDC_PL_FWEAK);
	CButton *pFstrong	= (CButton *) GetDlgItem(IDC_PL_FSTRONG);
	CButton *pFfib		= (CButton *) GetDlgItem(IDC_PL_FFIB);
	CButton *pFt		= (CButton *) GetDlgItem(IDC_PL_FT);
	CButton *pDC		= (CButton *) GetDlgItem(IDC_PL_DC);
	CButton *pDD		= (CButton *) GetDlgItem(IDC_PL_DD);
	CButton *pInflow	= (CButton *) GetDlgItem(IDC_PL_INFLOW);
	CButton *pU			= (CButton *) GetDlgItem(IDC_PL_U);

	int a = (enable) ? 1 : 0;
	
	pN->EnableWindow(a);
	pFweak->EnableWindow(a);
	pFstrong->EnableWindow(a);
	pFfib->EnableWindow(a);
	pFt->EnableWindow(a);
	pDC->EnableWindow(a);
	pDD->EnableWindow(a);
	pInflow->EnableWindow(a);
	pU->EnableWindow(a);
}

void CBeginNew::OnBnClickedFcEnabled()
{
	CButton *pEnabled	= (CButton *) GetDlgItem(IDC_FAC_ENABLED);
	switchFc(pEnabled->GetCheck() == 1);
}

void CBeginNew::switchFc(bool enable)
{
	CButton *pDx		= (CButton *) GetDlgItem(IDC_FAC_DX);
	CButton *pDt		= (CButton *) GetDlgItem(IDC_FAC_DT);
	CButton *pAlpha1	= (CButton *) GetDlgItem(IDC_FAC_ALPHA1);
	CButton *pAlpha2	= (CButton *) GetDlgItem(IDC_FAC_ALPHA2);
	CButton *pAlpha3	= (CButton *) GetDlgItem(IDC_FAC_ALPHA3);
	CButton *pAlpha4	= (CButton *) GetDlgItem(IDC_FAC_ALPHA4);
	CButton *pGamma1	= (CButton *) GetDlgItem(IDC_FAC_GAMMA1);
	CButton *pGamma2	= (CButton *) GetDlgItem(IDC_FAC_GAMMA2);
	CButton *pGamma3	= (CButton *) GetDlgItem(IDC_FAC_GAMMA3);
	CButton *pBeta1		= (CButton *) GetDlgItem(IDC_FAC_BETA1);
	CButton *pBeta2		= (CButton *) GetDlgItem(IDC_FAC_BETA2);
	CButton *pBeta3		= (CButton *) GetDlgItem(IDC_FAC_BETA3);
	CButton *pWB		= (CButton *) GetDlgItem(IDC_FAC_W_BEG);
	CButton *pWE		= (CButton *) GetDlgItem(IDC_FAC_W_END);
	CButton *pWV		= (CButton *) GetDlgItem(IDC_FAC_W_VAL);
	CButton *pk1		= (CButton *) GetDlgItem(IDC_FAC_K1);
	CButton *pk3		= (CButton *) GetDlgItem(IDC_FAC_K3);
	CButton *pC0		= (CButton *) GetDlgItem(IDC_FAC_C0);
	CButton *pT0		= (CButton *) GetDlgItem(IDC_FAC_T0);
	CButton *initFg		= (CButton *) GetDlgItem(IDC_FAC_INIT_FG);

	int a = (enable) ? 1 : 0;
	
	pDx->EnableWindow(a);
	pDt->EnableWindow(a);
	pAlpha1->EnableWindow(a);
	pAlpha2->EnableWindow(a);
	pAlpha3->EnableWindow(a);
	pAlpha4->EnableWindow(a);
	pGamma3->EnableWindow(a);
	pGamma3->EnableWindow(a);
	pGamma3->EnableWindow(a);
	pBeta1->EnableWindow(a);
	pBeta1->EnableWindow(a);
	pBeta1->EnableWindow(a);
	pWB->EnableWindow(a);
	pWE->EnableWindow(a);
	pWV->EnableWindow(a);
	pk1->EnableWindow(a);
	pk3->EnableWindow(a);
	pC0->EnableWindow(a);
	pT0->EnableWindow(a);
	initFg->EnableWindow(a);
}

//======================================================================
BEGIN_MESSAGE_MAP(CTrans, CDialog)
	ON_BN_CLICKED(ID_TRANS_CALC, &CTrans::OnBnClickedTransCalc)
	ON_BN_CLICKED(ID_TRANS_TRANS, &CTrans::OnBnClickedTransTrans)
	ON_BN_CLICKED(ID_TRANS_GET, &CTrans::OnBnClickedTransGet)
	ON_BN_CLICKED(ID_TRANS_SET, &CTrans::OnBnClickedTransSet)
	ON_BN_CLICKED(ID_TRANS_CLOSE, &CTrans::OnBnClickedTransClose)
END_MESSAGE_MAP()

//======================================================================
void CTrans::DoDataExchange(CDataExchange* pDX)
{
	DDX_Text(pDX, IDC_IN_D,			Calc->s.uni.D_in);
	DDX_Text(pDX, IDC_IN_G,			Calc->s.uni.G_in);
	DDX_Text(pDX, IDC_IN_M,			Calc->s.uni.m_in);
	DDX_Text(pDX, IDC_IN_MU,		Calc->s.uni.mu_in);
	DDX_Text(pDX, IDC_IN_NY,		Calc->s.uni.Ny);
	DDX_Text(pDX, IDC_IN_NX,		Calc->s.uni.Nx);
	DDX_Text(pDX, IDC_IN_RHO,		Calc->s.uni.rho_in);
	DDX_Text(pDX, IDC_IN_U,			Calc->s.uni.u_in);
	DDX_Text(pDX, IDC_IN_L,			Calc->s.uni.L_in);
	DDX_Text(pDX, IDC_IN_R,			Calc->s.uni.r_in);

	DDX_Text(pDX, IDC_OUT_D,		Calc->s.uni.D_out);
	DDX_Text(pDX, IDC_OUT_G,		Calc->s.uni.G_out);
	DDX_Text(pDX, IDC_OUT_M,		Calc->s.uni.m_out);
	DDX_Text(pDX, IDC_OUT_MU,		Calc->s.uni.mu_out);
	DDX_Text(pDX, IDC_OUT_RHO,		Calc->s.uni.rho_out);
	DDX_Text(pDX, IDC_OUT_U,		Calc->s.uni.u_out);
	DDX_Text(pDX, IDC_OUT_L,		Calc->s.uni.L_out);
	DDX_Text(pDX, IDC_OUT_R,		Calc->s.uni.r_out);

	DDX_Text(pDX, IDC_FM_OTHER_V,	Calc->s.uni.mFromOth);
	int kilo = 0, none = 0, centi = 0, mili = 0, micro = 0, nano = 0, pico = 0, other = 0;
	if(Calc->s.uni.mFrom == NONE)		none = 1;
	else if(Calc->s.uni.mFrom == CENTI)	centi = 1;
	else if(Calc->s.uni.mFrom == MILI)	mili = 1;
	else if(Calc->s.uni.mFrom == MICRO)	micro = 1;
	else if(Calc->s.uni.mFrom == NANO)	nano = 1;
	else if(Calc->s.uni.mFrom == PICO)	pico = 1;
	else					other = 1;
	DDX_Check(pDX, IDC_FM,			none);
	DDX_Check(pDX, IDC_FM_CENTI,	centi);
	DDX_Check(pDX, IDC_FM_MILI,		mili);
	DDX_Check(pDX, IDC_FM_MICRO,	micro);
	DDX_Check(pDX, IDC_FM_NANO,		nano);
	DDX_Check(pDX, IDC_FM_PICO,		pico);
	DDX_Check(pDX, IDC_FM_OTHER,	other);
	Calc->s.uni.mFrom = none ? NONE : centi ? CENTI : mili ? MILI : micro ? MICRO : nano ? NANO : pico ? PICO : other ? Calc->s.uni.mFromOth : NONE;

	DDX_Text(pDX, IDC_TM_OTHER_V,	Calc->s.uni.mToOth);
	kilo = 0, none = 0, centi = 0, mili = 0, micro = 0, nano = 0, pico = 0, other = 0;
	if(Calc->s.uni.mTo == NONE)		none = 1;
	else if(Calc->s.uni.mTo == CENTI)	centi = 1;
	else if(Calc->s.uni.mTo == MILI)	mili = 1;
	else if(Calc->s.uni.mTo == MICRO)	micro = 1;
	else if(Calc->s.uni.mTo == NANO)	nano = 1;
	else if(Calc->s.uni.mTo == PICO)	pico = 1;
	else					other = 1;
	DDX_Check(pDX, IDC_TM,			none);
	DDX_Check(pDX, IDC_TM_CENTI,	centi);
	DDX_Check(pDX, IDC_TM_MILI,		mili);
	DDX_Check(pDX, IDC_TM_MICRO,	micro);
	DDX_Check(pDX, IDC_TM_NANO,		nano);
	DDX_Check(pDX, IDC_TM_PICO,		pico);
	DDX_Check(pDX, IDC_TM_OTHER,	other);
	Calc->s.uni.mTo = none ? NONE : centi ? CENTI : mili ? MILI : micro ? MICRO : nano ? NANO : pico ? PICO : other ? Calc->s.uni.mToOth : NONE;

	DDX_Text(pDX, IDC_FG_OTHER_V,	Calc->s.uni.gFromOth);
	kilo = 0, none = 0, centi = 0, mili = 0, micro = 0, nano = 0, pico = 0, other = 0;
	if(Calc->s.uni.gFrom == KILO)		kilo = 1;
	else if(Calc->s.uni.gFrom == NONE)	none = 1;
	else if(Calc->s.uni.gFrom == MILI)	mili = 1;
	else if(Calc->s.uni.gFrom == MICRO)	micro = 1;
	else if(Calc->s.uni.gFrom == NANO)	nano = 1;
	else if(Calc->s.uni.gFrom == PICO)	pico = 1;
	else					other = 1;
	DDX_Check(pDX, IDC_FG_KILO,		kilo);
	DDX_Check(pDX, IDC_FG,			none);
	DDX_Check(pDX, IDC_FG_MILI,		mili);
	DDX_Check(pDX, IDC_FG_MICRO,	micro);
	DDX_Check(pDX, IDC_FG_NANO,		nano);
	DDX_Check(pDX, IDC_FG_PICO,		pico);
	DDX_Check(pDX, IDC_FG_OTHER,	other);
	Calc->s.uni.gFrom = kilo ? KILO : none ? NONE : mili ? MILI : micro ? MICRO : nano ? NANO : pico ? PICO : other ? Calc->s.uni.gFromOth : NONE;

	DDX_Text(pDX, IDC_TG_OTHER_V,	Calc->s.uni.gToOth);
	kilo = 0, none = 0, centi = 0, mili = 0, micro = 0, nano = 0, pico = 0, other = 0;
	if(Calc->s.uni.gTo == KILO)		kilo = 1;
	else if(Calc->s.uni.gTo == NONE)	none = 1;
	else if(Calc->s.uni.gTo == MILI)	mili = 1;
	else if(Calc->s.uni.gTo == MICRO)	micro = 1;
	else if(Calc->s.uni.gTo == NANO)	nano = 1;
	else if(Calc->s.uni.gTo == PICO)	pico = 1;
	else					other = 1;
	DDX_Check(pDX, IDC_TG_KILO,		kilo);
	DDX_Check(pDX, IDC_TG,			none);
	DDX_Check(pDX, IDC_TG_MILI,		mili);
	DDX_Check(pDX, IDC_TG_MICRO,	micro);
	DDX_Check(pDX, IDC_TG_NANO,		nano);
	DDX_Check(pDX, IDC_TG_PICO,		pico);
	DDX_Check(pDX, IDC_TG_OTHER,	other);
	Calc->s.uni.gTo = kilo ? KILO : none ? NONE : mili ? MILI : micro ? MICRO : nano ? NANO : pico ? PICO : other ? Calc->s.uni.gToOth : NONE;

	DDX_Text(pDX, IDC_FS_OTHER_V,	Calc->s.uni.sFromOth);
	kilo = 0, none = 0, centi = 0, mili = 0, micro = 0, nano = 0, pico = 0, other = 0;
	if(Calc->s.uni.sFrom == NONE)	none = 1;
	else if(Calc->s.uni.sFrom == MILI)	mili = 1;
	else if(Calc->s.uni.sFrom == MICRO)	micro = 1;
	else					other = 1;
	DDX_Check(pDX, IDC_FS,			none);
	DDX_Check(pDX, IDC_FS_MILI,		mili);
	DDX_Check(pDX, IDC_FS_MICRO,	micro);
	DDX_Check(pDX, IDC_FS_OTHER,	other);
	Calc->s.uni.sFrom = kilo ? KILO : none ? NONE : mili ? MILI : micro ? MICRO : nano ? NANO : pico ? PICO : other ? Calc->s.uni.sFromOth : NONE;

	DDX_Text(pDX, IDC_TS_OTHER_V,	Calc->s.uni.sToOth);
	kilo = 0, none = 0, centi = 0, mili = 0, micro = 0, nano = 0, pico = 0, other = 0;
	if(Calc->s.uni.sTo == NONE)	none = 1;
	else if(Calc->s.uni.sTo == MILI)	mili = 1;
	else if(Calc->s.uni.sTo == MICRO)	micro = 1;
	else					other = 1;
	DDX_Check(pDX, IDC_TS,			none);
	DDX_Check(pDX, IDC_TS_MILI,		mili);
	DDX_Check(pDX, IDC_TS_MICRO,	micro);
	DDX_Check(pDX, IDC_TS_OTHER,	other);
	Calc->s.uni.sTo = kilo ? KILO : none ? NONE : mili ? MILI : micro ? MICRO : nano ? NANO : pico ? PICO : other ? Calc->s.uni.sToOth : NONE;
}
void CTrans::OnBnClickedTransCalc()
{
	UpdateData(true);
	
	DBL u = Calc->s.uni.u_in;
	DBL mu = Calc->s.uni.mu_in;
	DBL D = Calc->s.uni.D_in;
	DBL L = Calc->s.uni.L_in;
	int Ny = Calc->s.uni.Ny;
	int Nx = L * Ny / D;
	int N = Nx*Ny;
	DBL rho = Calc->s.uni.rho_in;
	DBL r = D * pow(3./(pi*4), 1./3) / Ny;
	DBL m = pow(D/Ny, 3)*rho;

	Calc->s.uni.Nx = Nx;
	Calc->s.uni.r_in = r;
	Calc->s.uni.m_in = m;
	Calc->s.uni.G_in = 12 * u * mu / (rho * pow(D, 2));
	
	UpdateData(false);
}

void CTrans::OnBnClickedTransTrans()
{
	UpdateData(true);

	DBL m = Calc->s.uni.mFrom - Calc->s.uni.mTo;
	DBL g = Calc->s.uni.gFrom - Calc->s.uni.gTo;
	DBL sec = Calc->s.uni.sFrom - Calc->s.uni.sTo;

	Calc->s.uni.m_out = Calc->s.uni.m_in * pow(10, g);
	Calc->s.uni.u_out = Calc->s.uni.u_in * pow(10, m - sec);
	Calc->s.uni.D_out = Calc->s.uni.D_in * pow(10, m);
	Calc->s.uni.rho_out = Calc->s.uni.rho_in * pow(10, g - 3 * m);
	Calc->s.uni.mu_out = Calc->s.uni.mu_in * pow(10, g - m - sec);
	Calc->s.uni.G_out = Calc->s.uni.G_in * pow(10, m - 2 * sec);
	Calc->s.uni.r_out = Calc->s.uni.r_in * pow(10, m);
	Calc->s.uni.L_out = Calc->s.uni.L_in * pow(10, m);
	
	UpdateData(false);
}

void CTrans::OnBnClickedTransGet()
{
}

void CTrans::OnBnClickedTransSet()
{
	UpdateData(false);

	Calc->s.spc.size.x = Calc->s.uni.L_out;
	Calc->s.spc.size.y = Calc->s.uni.D_out;
	Calc->s.atom.n.x = Calc->s.uni.Nx;
	Calc->s.atom.n.y = Calc->s.uni.Ny;
	Calc->s.atom.m = Calc->s.uni.m_out;
	Calc->s.met.G.x = Calc->s.uni.G_out;
	Calc->s.atom.r = Calc->s.uni.r_out;
	Calc->s.atom.rShow = Calc->s.uni.r_out / 2;
}

void CTrans::OnBnClickedTransClose()
{
	CDialog::OnCancel();
}

//======================================================================
BEGIN_MESSAGE_MAP(COptions, CDialog)
	ON_BN_CLICKED(ID_OPT_OK,		&COptions::OnBnClickedOptOK)
	ON_BN_CLICKED(ID_OPT_CANCEL,	&COptions::OnBnClickedOptCancel)
	ON_BN_CLICKED(IDC_OPT_COLPLA, &COptions::OnBnClickedOptColpla)
	ON_BN_CLICKED(IDC_OPT_COLPLT, &COptions::OnBnClickedOptColplt)
	ON_BN_CLICKED(IDC_OPT_COLFIB, &COptions::OnBnClickedOptColfib)
	ON_BN_CLICKED(IDC_OPT_COLF1, &COptions::OnBnClickedOptColf1)
	ON_BN_CLICKED(IDC_OPT_COLF2, &COptions::OnBnClickedOptColf2)
	ON_BN_CLICKED(ID_OPT_BROWSE, &COptions::OnBnClickedOptBrowse)
END_MESSAGE_MAP()

void COptions::DoDataExchange(CDataExchange* pDX)
{
	DDX_Text(pDX, IDC_OPT_TEXT_SIZE,	Calc->s.opt.textSize);
	DDX_Text(pDX, IDC_OPT_TEXT_SPACE,	Calc->s.opt.textSpace);
	DDX_Text(pDX, IDC_OPT_THREADS,		Calc->s.opt.proc);
	DDX_Text(pDX, IDC_OPT_IMGW,			Calc->s.opt.imgW);

	CString csn(Calc->s.opt.outputPath.c_str());
	DDX_Text(pDX, IDC_OPT_OUTPUT, csn);
	Calc->s.opt.outputPath = string((LPCTSTR)csn);
		
	int tmp = boolToInt(Calc->s.opt.showPla);
	DDX_Check(pDX, IDC_OPT_PLA, tmp);
	Calc->s.opt.showPla = intToBool(tmp); 

	tmp = boolToInt(Calc->s.opt.showPlt);
	DDX_Check(pDX, IDC_OPT_PLT, tmp);
	Calc->s.opt.showPlt = intToBool(tmp);

	tmp = boolToInt(Calc->s.opt.showFib);
	DDX_Check(pDX, IDC_OPT_FIB, tmp);
	Calc->s.opt.showFib = intToBool(tmp);

	tmp = boolToInt(Calc->s.opt.showF1);
	DDX_Check(pDX, IDC_OPT_F1, tmp);
	Calc->s.opt.showF1 = intToBool(tmp);

	tmp = boolToInt(Calc->s.opt.showF2);
	DDX_Check(pDX, IDC_OPT_F2, tmp);
	Calc->s.opt.showF2 = intToBool(tmp);

	setCol();
}

void COptions::setCol()
{
	CMFCButton *b1 = (CMFCButton *) GetDlgItem(IDC_OPT_COLPLA);
	b1->SetFaceColor(cpla);
	CMFCButton *b2 = (CMFCButton *) GetDlgItem(IDC_OPT_COLPLT);
	b2->SetFaceColor(cplt);
	CMFCButton *b3 = (CMFCButton *) GetDlgItem(IDC_OPT_COLFIB);
	b3->SetFaceColor(cfib);
	CMFCButton *b4 = (CMFCButton *) GetDlgItem(IDC_OPT_COLF1);
	b4->SetFaceColor(cf1);
	CMFCButton *b5 = (CMFCButton *) GetDlgItem(IDC_OPT_COLF2);
	b5->SetFaceColor(cf2);
}

void COptions::OnBnClickedOptOK()
{
	CString cs;
	GetDlgItemText(IDC_OPT_OUTPUT, cs);
	string s = (LPCTSTR)cs;
	
	bool good = true;
	
	if(!FS::dirExists(s)) 
	{
		good = FS::dirCreate(s);
	}

	if(good)
	{
		CDialog::OnOK();
		Calc->s.opt.outputPath = s;
		Calc->s.opt.cpla = cpla;
		Calc->s.opt.cplt = cplt;
		Calc->s.opt.cfib = cfib;
		Calc->s.opt.cf1 = cf1;
		Calc->s.opt.cf2 = cf2;

		IO::saveProps(Calc->s);
	}
}

void COptions::OnBnClickedOptCancel()
{
	CDialog::OnCancel();
}

void COptions::OnBnClickedOptColpla()
{
	CColorDialog d;
	if(d.DoModal() == IDOK)
	{
		cpla = d.GetColor();
		CMFCButton *b = (CMFCButton *) GetDlgItem(IDC_OPT_COLPLA);
		b->SetFaceColor(cpla);
	}
}


void COptions::OnBnClickedOptColplt()
{
	CColorDialog d;
	if(d.DoModal() == IDOK)
	{
		cplt = d.GetColor();
		CMFCButton *b = (CMFCButton *) GetDlgItem(IDC_OPT_COLPLT);
		b->SetFaceColor(cplt);
	}
}


void COptions::OnBnClickedOptColfib()
{
	CColorDialog d;
	if(d.DoModal() == IDOK)
	{
		cfib = d.GetColor();
		CMFCButton *b = (CMFCButton *) GetDlgItem(IDC_OPT_COLFIB);
		b->SetFaceColor(cfib);
	}
}


void COptions::OnBnClickedOptColf1()
{
	CColorDialog d;
	if(d.DoModal() == IDOK)
	{
		cf1 = d.GetColor();
		CMFCButton *b = (CMFCButton *) GetDlgItem(IDC_OPT_COLF1);
		b->SetFaceColor(cf1);
	}
}


void COptions::OnBnClickedOptColf2()
{
	CColorDialog d;
	if(d.DoModal() == IDOK)
	{
		cf2 = d.GetColor();
		CMFCButton *b = (CMFCButton *) GetDlgItem(IDC_OPT_COLF2);
		b->SetFaceColor(cf2);
	}
}


void COptions::OnBnClickedOptBrowse()
{
	CString cs;
	GetDlgItemText(IDC_OPT_OUTPUT, cs);
	string s = (LPCTSTR)cs;
	//CFileDialog f(true,); 

	CoInitialize(0);

	CFileDialog f(true, NULL, NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, NULL, NULL, 0, true);
	IFileOpenDialog *openDlgPtr = f.GetIFileOpenDialog();
	if(openDlgPtr != NULL)
	{
		openDlgPtr->SetOptions(FOS_PICKFOLDERS);
		openDlgPtr->Release();
		openDlgPtr->Show(NULL);
		IShellItem *pShellItemSelected;
		if(openDlgPtr->GetResult(&pShellItemSelected) == S_OK) s = f.GetPathName();
		SetDlgItemText(IDC_OPT_OUTPUT, s.c_str());
	}

	CoUninitialize();
		
}


//======================================================================
BEGIN_MESSAGE_MAP(CProgress, CDialog)
END_MESSAGE_MAP()

//======================================================================
void CProgress::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_DESC, desc);
	DDX_Control(pDX, IDC_PROGRESS1, bar);
}


