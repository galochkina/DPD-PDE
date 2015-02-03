#ifndef CONST_H
#define	CONST_H

namespace Const
{
	static const string SAVT	= "#SAVT";	//save type
	static const string SIMT	= "#SIMT";	//simulation type
	static const string SPDE	= "#SPDE";	//SetupPde
	static const string SSPC	= "#SSPC";	//SetupSpace	
	static const string SMET	= "#SMET";	//SetupMethod
	static const string SATM	= "#SATM";	//SetupAtoms
	static const string SCPL	= "#SCPL";	//SetupCPlt
	static const string SCER	= "#SCER";	//SetupCRBC
	static const string SCEL	= "#SCEL";	//SetupCells
	static const string SPLT	= "#SPLT";	//SetupPlatelets
	static const string MP	    = "#MP";	//SetupMP
	static const string SRUN	= "#SRUN";	//SetupRun
	static const string SANL	= "#SANL";	//SetupAnalysis
	static const string SDRW	= "#SDRW";	//SetupDraw
	static const string SUNI	= "#SUNI";	//SetupUnit
	static const string ALIS	= "#ALIS";	//atom list
	static const string SLIS	= "#SLIS";	//smart list
	static const string ATOM	= "#ATOM";	//atom
	static const string GOBJ	= "#GOBJ";	//generic object
	static const string CELL	= "#CELL";	//cell
	static const string FIBR	= "#FIBR";	//fibrin setup
	static const string FIBX	= "#FIBX";	//fibrin X (concentration matrix)
	static const string CLCK	= "#CLCK";	//Clock
	static const string MANL	= "#MANL";	//measurements
	static const string STAT	= "#STAT";	//stat
	static const string AVGF	= "#AVGF";	//pfSum
	static const string CLTH	= "#CLTH";	//clot size history
	static const string ENGH	= "#ENGH";	//energy history
	static const string PCON	= "#PCON";	//platelet pair list
	static const string PPAR	= "#PPAR";	//platelet pair
	static const string DVEL	= "#DVEL";	//DataVel
	static const string DVDE	= "#DVDE";	//DataVelDen
	static const string DHIS	= "#DHIS";	//DataHist
	static const string DFLW	= "#DFLW";	//DataFlow
	static const string DFDE	= "#DFDE";	//DataForceDensity
	static const string AVDE	= "#AVDE";	//AnlVelDen
	static const string AVDP	= "#AVDP";	//AnlVelDenPde
	static const string ACLT	= "#ACLT";	//AnlClotTotal
	static const string ACL2	= "#ACL2";	//AnlClotFib
	static const string AFLW	= "#AFLW";	//AnlFlow
	static const string AFDE	= "#AFDE";	//AnlForceDen
	static const string STA2	= "#STA2";	//Stats2
	static const string PDE1	= "#PDE1";	//T
	static const string PDE2	= "#PDE2";	//Fg
	static const string PDE3	= "#PDE3";	//Fp
	static const string PDE4	= "#PDE1";	//Tf
	static const string PDE5	= "#PDE2";	//Tm

	static const string ANL		= "#ANL";	//Fp


	static const string S02		= "S02";
	static const string S03		= "S03";
	static const string S04		= "S04";
	static const string STP		= "STP";

	static const int KILO		= 3;
	static const int NONE		= 0;
	static const int CENTI		= -2;
	static const int MILI		= -3;
	static const int MICRO		= -6;
	static const int NANO		= -9;
	static const int PICO		= -12;
}

#endif