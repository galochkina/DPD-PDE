#ifndef DRAW3_H
#define	DRAW3_H

class Graph2
{
public:
	Graph2(int uL = 800) { gr = NULL; unitL = uL; dataUsed = 0; }

	void setData(DBL *val1, int n)	{ v1.Set(val1, n); dataUsed = 1; }
	void setData(Vec<DBL> &val1)		{ setData(val1.arr(), val1.n); }
	void setData(DBL *val1, int m, int n)	{ v1.Set(val1, m, n); dataUsed = 1; }
	void setData(Mat<DBL> &val1)				{ setData(val1.arr(), val1.m, val1.n); }
	void setData(DBL *val1, DBL *val2, int n)	{ v1.Set(val1, n); v2.Set(val2, n); dataUsed = 2; }
	void setData(Vec<DBL> &val1, Vec<DBL> &val2)	{ if(val1.n != val2.n) return; setData(val1.arr(), val2.arr(), val1.n); }
	void setData(DBL *val1, DBL *val2, DBL *val3, int n)	{ v1.Set(val1, n); v2.Set(val2, n); v3.Set(val3, n); dataUsed = 3; }
	void setData(Vec<DBL> &val1, Vec<DBL> &val2, Vec<DBL> &val3)	{ if((val1.n != val2.n) || (val1.n != val3.n)) return; setData(val1.arr(), val2.arr(), val3.arr(), val1.n); }
	void setData(DBL *val1, DBL *val2, int m, int n)	{ v1.Set(val1, m, n); v2.Set(val2, m, n); dataUsed = 2; }
	void setData(Mat<DBL> &val1, Mat<DBL> &val2)		{ if(val1.m != val2.m || val1.n != val2.n) return; setData(val1.arr(), val2.arr(), val1.m, val1.n); }
	void setData(DBL *val1, int n, string fn, DBL from = 0, DBL to = 0)	{ v1.Set(val1, n); v2 = mglData(n); v2.Fill(from, to, 'x'); v2.Modify(fn.c_str()); dataUsed = 2; }
	void setData(Vec<DBL> &val1, string fn, DBL from = 0, DBL to = 0)	{ setData(val1.arr(), val1.n, fn, from, to); }

	void setLabel(string x, string y, string z = "") { lx = x; ly = y; lz = z; }
	void setRange(VecF2 x = VecF20, VecF2 y  = VecF20, VecF2 z = VecF20) { rx = x, ry = y; rz = z; }
	
	void create(string sFolder, string sName, int step = -1, DBL x = 1, DBL y = 1) { if(x >= y){ if(x > 3){ y = y / x * 3; x = 3; } }else if(y > x){ if(y > 3){ x = x / y * 3; y = 3; } } if(x*y<2)	size = VecL2(unitL * x * 2, unitL * y * 2); else size = VecL2(unitL * x, unitL * y);  name = FS::getFN(sFolder, sName, "png", step); gr = new mglGraph(0, size.x, size.y);}

	void plot(VecL3 sp = VecL30) { prepare(sp); if(dataUsed >= 1) gr->Plot(v1, "b-2"); if(dataUsed >= 2) gr->Plot(v2, "r-2"); if(dataUsed >= 3) gr->Plot(v3, "g-2"); }
	void vect(VecL3 sp = VecL30) { prepare(sp); gr->Vect(v1, v2); }
	void dens(VecL3 sp = VecL30) { prepare(sp); gr->Dens(v1); gr->Colorbar(); }
	void surf(VecL3 sp = VecL30, bool mesh = true) { prepare(sp, true); if(mesh) gr->Surf(v1, "BbcyrR#"); else gr->Surf(v1, "BbcyrR");}

	void draw() { gr->WritePNG(name.c_str()); free(); }

private:
	mglGraph *gr;
	mglData v1, v2, v3;
	int dataUsed;
	int unitL;
	string name;
	VecL2 size;
	string lx, ly, lz;
	VecF2 rx, ry, rz;

	void free()	{ if(gr != NULL) delete gr;	dataUsed = 0; }

	void prepare(VecL3 sp = VecL30, bool rotate = false)
	{
		if(sp != VecL30) gr->SubPlot(sp.x, sp.y, sp.z);
		if(rx != VecF20 && ry != VecF20) gr->SetRanges(rx.x, rx.y, ry.x, ry.y, rz.x, rz.y);
		if(rx != VecF20) gr->SetTicks('x', (rx.y - rx.x) / 4);  
		if(ry != VecF20) gr->SetTicks('y', (ry.y - ry.x) / 4);  
		if(rz != VecF20) gr->SetTicks('z', (rz.y - rz.x) / 4);

		if(rotate) gr->Rotate(60,40); else gr->Rotate(0,0);
		gr->Box();
		gr->Label('x', lx.c_str());
		gr->Label('y', ly.c_str());
		gr->Label('z', lz.c_str());
		gr->Axis();
	}

};

#endif
