#ifndef MYCALC_H
#define	MYCALC_H

#include "matrix.h"
#include "helper.h"
#include "data.h"
#include "setup.h"
#include "filesystem.h"
#include "pde.h"
#include "objects.h"
#include "problem.h"
#include "analysis.h"
#include "print.h"

#include "interactions.h"
#include "move.h"
#include "save.h"
#include "draw.h"
#include "analysis2.h"
#include "MyCalc.h"

class  MyCalc : public Calculus  
{ 
public:

	//*********************
	//***** VARIABLES *****
	//*********************
	

	//SIMULATION
	AtomList	al;			//atom list
	PPL			platPairs;	//platelet pair list
	LC			cells;		//cell list
	Mesh		mesh;		//mesh
	Clock		c;			//clock
	BoundCond	bc;			//boundary conditions
	Stats		stats;		//analysis data
	Stats2		stats2;		//analysis data
	Draw		draw;
	FileSys		fs;

	PDE			*pde;
	CRDE		fib;
	
	//*******************		
	//***** METHODS *****
	//*******************

	//constructors
	MyCalc();		

	//animation methods
	void	setAtoms(bool keepCells = false);

	vector<int>	testCPU();
	void	testGPU(bool gpu);
	//void	applyForceCPU(Atom &a);

	//script
	void	runScript();
	void	setParamsScript(vector<string> &vec);

	//system methods
	void	clear();
	void	createDirectories();
	void	init();
	void	calculateBoxPair(BoxPair &bp);
	void	nextStep(bool preRun);
	void	setMaxRadius();
	void	save(string filePath = "");
	bool	load(string fileName, bool newInstance);
	bool	updateU();

	//Nikolay
	virtual bool    ReadCalculus (FILE *file, bool INOUT) {file;INOUT; return true;};
	virtual void    WRITE(){};
	virtual bool    IsValidPaint() { return al.atoms.size() != 0 || cells.size() != 0; };
	virtual void	ViewReset(){};

	//destructor
	~MyCalc()    
	 {
		if(pDialog) delete pDialog;
		if(pDialog2) delete pDialog2;
		if(pDialog3) delete pDialog3;
		if(pDialogOpt) delete pDialogOpt;
		if(pde) delete pde;
		stats.destroy();
	 }

	void MySetWindowText() 
	{ 
		CString Text;
		vector<int> time = c.getTime();
		vector<int> timeR = myStp.is_DPD() ? c.getRemainingTime(myStp.run.stop / myStp.run.dt) : 
			c.getRemainingTime(myStp.pde.stop / myStp.pde.dt);

		string delim = "   ";
		string str = "             simulation time = %10.6f" + delim
			+ "real time = %02d:%02d:%02d ( %02d:%02d:%02d )" + delim
			+ "steps/sec = %.0f ( %.0f )" + delim
			+ "mu = %10.3f" + delim
			+ "uErr = %.3f%%" + delim
			+ "uMax : uAvg = %.3f ( %.3f , %.3f )";

		Text.Format(_T(str.c_str()), 
			myStp.run.t, 
			time[0], time[1], time[2], timeR[0], 
			timeR[1], timeR[2], 
			c.getRate(), c.getAvgRate(), 
			stats.m.mu, 
			stats.m.uErr, 
			stats.m.uMA, 
			stats.m.uMax, 
			stats.m.uAvg); 

		SetWindowText(myStp.run.name.c_str() + Text);                           
	} //

	static UINT MyCalc::ThreadProc(LPVOID pArg)
	{
		MyCalc *myPtr = reinterpret_cast<MyCalc *>(pArg);	 
		if(myStp.is_DPD()) 
		{
			myPtr->RunBodyDPD();
		}
		else if(myStp.is_PDE()) 
		{
			myPtr->RunBodyPDE();
		}
		return 0;
	}

	static UINT MyCalc::ThreadProcScript(LPVOID pArg)
	{
		MyCalc *myPtr = reinterpret_cast<MyCalc *>(pArg);	 
		myPtr->RunScriptBody();
		return 0;
	}

protected:
	//animation methods
	afx_msg void	OnFileNew();
	afx_msg void	OnTrans();
	afx_msg void	OnPde();
	afx_msg void	OnToolsOptions();
	afx_msg void	OnPaint();
	DECLARE_MESSAGE_MAP()
};

#endif



