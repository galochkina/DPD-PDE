#ifndef DIALOG_H
#define	DIALOG_H

#include "afxcmn.h"
#include "afxwin.h"

using namespace Const;

class CBeginNew : public CDialog 
{  
public:
	CBeginNew(): CDialog() {};

	enum { IDD = IDD_BEGIN_NEW };
	virtual void DoDataExchange(CDataExchange* pDX); 
	afx_msg void OnBnClickedButton1();
	afx_msg void OnBnLoad();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedSavestp();
	afx_msg void OnBnClickedButton4();
	afx_msg void OnBnClickedButton5();
	afx_msg void OnBnClickedPlEnabled();
	afx_msg void OnBnClickedFcEnabled();
	afx_msg void OnBnClickedInit();
	void switchPl(bool enable);
	void switchFc(bool enable);
	DECLARE_MESSAGE_MAP();	
	afx_msg void OnBnClickedImg();
	afx_msg void OnBnClickedTxtParams();
};

class CTrans : public CDialog 
{  
public:
	CTrans(): CDialog() {};

	enum { IDD = IDD_TRANS };
	virtual void DoDataExchange(CDataExchange* pDX); 
	DECLARE_MESSAGE_MAP();	
	afx_msg void OnBnClickedTransCalc();
	afx_msg void OnBnClickedTransTrans();
	afx_msg void OnBnClickedTransGet();
	afx_msg void OnBnClickedTransSet();
	afx_msg void OnBnClickedTransClose();
};

class COptions : public CDialog 
{  
public:
	COptions(): CDialog() {};

	COLORREF cpla;
	COLORREF cplt;
	COLORREF mp;
	COLORREF cfib;
	COLORREF cf1;
	COLORREF cf2;

	enum { IDD = IDD_OPT };
	virtual void DoDataExchange(CDataExchange* pDX); 
	afx_msg void OnBnClickedOptOK();
	afx_msg void OnBnClickedOptCancel();
	DECLARE_MESSAGE_MAP();	
	afx_msg void OnBnClickedOptColpla();
	afx_msg void OnBnClickedOptColplt();
	afx_msg void OnBnClickedOptColfib();
	afx_msg void OnBnClickedOptColf1();
	afx_msg void OnBnClickedOptColf2();
	void setCol();
	afx_msg void OnBnClickedOptBrowse();
};

class CProgress : public CDialog 
{  
public:
	CStatic desc;
	CProgressCtrl bar;

	CProgress(): CDialog() {};

	enum { IDD = IDD_PROGRESS };
	virtual void DoDataExchange(CDataExchange* pDX); 
	DECLARE_MESSAGE_MAP();	
};

#endif