#ifndef CALCULUS_H
#define	CALCULUS_H

class  Calculus : public CMyFrame  
{ 
public:

	//*********************
	//***** VARIABLES *****
	//*********************
	
	//DIALOGS
	CBeginNew*		pDialog;
	CTrans*			pDialog2;
	COptions*		pDialogOpt;
	CProgress*		pDialogProg;

	//SIMULATION
	AtomList	al;			//atom list
	PPL			platPairs;		//platelet pair list
	Mesh		mesh;			//mesh
	Clock		c;				//clock
	BoundCond	bc;				//boundary conditions
	Setup		s;
	Draw		draw;
	Anl			anl;

	PdeSys		ps;

	bool		cleared;
	
	//*******************		
	//***** METHODS *****
	//*******************

	//constructors
	Calculus();		

	//animation methods
	void	Run();
	void	RunBodyDPD();
	void	setAtoms();
	void	progressBeg(string title, int totalSteps);
	void	progressEnd();
	void	setupHotChange();

	//system methods
	void	clear();
	void	init();
	void	calculateBoxPair(BoxPair &bp);
	void	calcAtomF();
	void	calcOtherF();
	void	moveAtoms();
	void	setOldPlat();
	void	nextStep();
	void	setMaxRadius();
	void	save(bool justStp, string filePath = "");
	bool	load(string fileName);

	//Nikolay
	virtual bool    ReadCalculus (FILE *file, bool INOUT) {file;INOUT; return true;};
	virtual void    WRITE(){};
	virtual bool    IsValidPaint() { return al.obj.size() != 0; };
	virtual void	ViewReset(){};

	//destructor
	~Calculus()
	 {
		if(pDialog) delete pDialog;
		if(pDialog2) delete pDialog2;
		if(pDialogOpt) delete pDialogOpt;
		if(pDialogProg) delete pDialogProg;
	 }

	void MySetWindowText() 
	{ 
		CString Text;
		vector<int> time = c.getTime();
		vector<int> timeR = c.getRemainingTime(s.run.stop / s.run.dt);
		string delim = "   ";
		string str = "             simulation time = %10.6f" + delim
			+ "real time = %02d:%02d:%02d ( %02d:%02d:%02d )" + delim
			+ "steps/sec = %.0f ( %.0f )";

		Text.Format(_T(str.c_str()), 
			s.run.t, 
			time[0], time[1], time[2], timeR[0], 
			timeR[1], timeR[2], 
			c.getRate(), c.getAvgRate());
		
		SetWindowText(s.run.name.c_str() + Text);                               
	} //

	static UINT Calculus::ThreadProc(LPVOID pArg)
	{
		Calculus *myPtr = reinterpret_cast<Calculus *>(pArg);	 
		myPtr->RunBodyDPD();
		return 0;
	}

protected:
	//animation methods
	afx_msg void	OnFileNew();
	afx_msg void	OnTrans();
	afx_msg void	OnToolsOptions();
	afx_msg void	OnPaint();
	DECLARE_MESSAGE_MAP()
};

#endif



