#ifndef DRAW2_H
#define	DRAW2_H
/*
class Graph
{
public:
	Graph()
	{
		gr = NULL;
		//unitL = 100;
		//unitL = 300;
		unitL = 800;
		//unitL = 1200;
		size = VecL2(unitL, unitL);
		rx = ry = rz = VecF20;
		lx = ly = lz = "";
		name = "graph.png";
		type = -1;
		nLines = 1;
		points = false;
	}

	void create()
	{
		gr = new mglGraph(0, size.x, size.y);
	}

	void put(VecL3 subplot = VecL30)
	{
		if(subplot != VecL30) gr->SubPlot(subplot.x, subplot.y, subplot.z);
		if(rx != VecF20 && ry != VecF20) gr->SetRanges(rx.x, rx.y, ry.x, ry.y, rz.x, rz.y);
		if(rx != VecF20) gr->SetTicks('x', (rx.y - rx.x) / 4);  
		if(ry != VecF20) gr->SetTicks('y', (ry.y - ry.x) / 4);  
		if(rz != VecF20) gr->SetTicks('z', (rz.y - rz.x) / 4);

		if(type == DRAW_SURF) gr->Rotate(60,40);
		gr->Box();
		gr->Label('x', lx.c_str());
		gr->Label('y', ly.c_str());
		if(type == DRAW_SURF) gr->Label('z', lz.c_str());
		gr->Axis();
		//gr->Title(title.c_str());
		
		if(type == DRAW_PLOT) 
		{
			if(nLines == 1)
				gr->Plot(val1, "b-2");
			else if(nLines == 2)
				gr->Plot(val1, val2, "r-2");
		}
		else if(type == DRAW_SURF) 
		{
			mgl_set_meshnum(gr, 10);
			gr->Surf(val1, "BbcyrR#");
		}
		else if(type == DRAW_DENS) 
		{
			gr->Dens(val1);
			gr->Colorbar(); 
		}
		else if(type == DRAW_VECT) 
		{
			gr->Vect(val1, val2); 
		}		
	}

	void draw()
	{
		gr->WritePNG(name.c_str());
		free();
	}

	void free()
	{
		if(gr != NULL) delete gr;
	}

	//setters
	void setLabel(string x, string y = "", string z = "") { lx = x; ly = y; lz = z; }
	void setName(string s) { name = s; }
	void setName(string sFolder, string sName, int step = -1) { name = FS::getFN(sFolder, sName, "png", step); }
	void setRange(VecF2 x, VecF2 y, VecF2 z = VecF20) { if(x != VecF20) { rx = x; ry = y; rz = z; } }
	//void setSize(int x, int y, int unitLength = 1200) { unitL = unitLength; int m = 1, n = 1; if(x > y) m = x / y; else if(x < y) n = y / x; size = VecL2(m * unitL, n * unitL); }
	//void setSize(int x, int y, int imgX, int imgY, int unitLength = 1200) { unitL = unitLength; int m = 1, n = 1; if(x > y) m = x / y; else if(x < y) n = y / x; size = VecL2(m * imgX * unitL, n * imgY * unitL); }
	//void setSize(VecL2 s) { size = s; }
	void setType(int t) { type = t; }
	void setTitle(string s = "") { title = s; }

	void setSize(DBL x, DBL y) 
	{ 
		if(x > y && x > 3) 
		{ 
			if(x>3)
			{
				y = y / x * 3;
				x = 3;
			}
		}
		else if(y > x)
		{
			if(y>3)
			{
				x = x / y * 3;
				y = 3;
			}
		}

		size = VecL2(unitL * x, unitL * y);
	}

	void setData(int m, int n, DBL* valA, DBL* valB, DBL* ptsA, DBL* ptsB, bool freeMemory) 
	{ 
		if(type == DRAW_PLOT)
		{
			if(valA != NULL) 
			{
				val1.Set(valA, m);
				if(ptsA != NULL) 
				{
					pts1.Set(ptsA, m);
					points = true;
				}
				if(valB != NULL) 
				{
					nLines = 2;
					val2.Set(valB, m);
					if(ptsB != NULL) pts2.Set(ptsB, n);
					else if(points) pts2.Set(ptsA, m);
				}
			}			
		}
		else
		{
			if(valA != NULL) val1.Set(valA, m, n);
			if(valB != NULL) val2.Set(valB, m, n);
		}

		if(freeMemory)
		{
			if(valA != NULL) delete[] valA;
			if(valB != NULL) delete[] valB;
			if(ptsA != NULL) delete[] ptsA;
			if(ptsB != NULL) delete[] ptsB;
		}
	}


	//static draw calls
	//PLOT

	static void drawPlot(string sTitle, string folder, string name, int step, DBL coeff, Vec<DBL> &val, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", val.n, 0, val.arr(coeff), NULL, NULL, NULL); } 

	static void drawPlot(string sTitle, string folder, string name, int step, DBL *val, int n, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", n, 0, val, NULL, NULL, NULL); } 

	static void drawPlot(string sTitle, string folder, string name, int step, DBL coeff, Vec<DBL> &val, Vec<DBL> &pts, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", val.n, 0, val.arr(coeff), NULL, pts.arr(), NULL); } 

	static void drawPlot(string sTitle, string folder, string name, int step, DBL coeff, Vec<DBL> &val1, Vec<DBL> &val2, Vec<DBL> &pts, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", val1.n, 0, val1.arr(coeff), val2.arr(coeff), pts.arr(), NULL); } 

	static void drawPlot(string sTitle, string folder, string name, int step, DBL coeff, Vec<DBL> &val1, Vec<DBL> &val2, Vec<DBL> &pts1, Vec<DBL> &pts2, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", val1.n, val2.n, val1.arr(coeff), val2.arr(coeff), pts1.arr(), pts2.arr()); } 

	//SURF
	static void drawSurf(string sTitle, string folder, string name, int step, DBL coeff, Mat<DBL> &val, VecF2 rx = VecF20, VecF2 ry = VecF20, VecF2 rz = VecF20, string lx = "", string ly = "",  string lz = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_SURF, rx, ry, rz, lx, ly, lz, val.m, val.n, val.arr(coeff), NULL, NULL, NULL); } 

	//DENS
	static void drawDens(string sTitle, string folder, string name, int step, DBL coeff, Mat<DBL> &val, VecF2 rx = VecF20, VecF2 ry = VecF20, VecF2 rz = VecF20, string lx = "", string ly = "",  string lz = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_DENS, rx, ry, rz, lx, ly, lz, val.m, val.n, val.arr(coeff), NULL, NULL, NULL); } 

	//VECT
	static void drawVect(string sTitle, string folder, string name, int step, DBL coeff, Mat<DBL> &val1, Mat<DBL> &val2, VecF2 rx = VecF20, VecF2 ry = VecF20, VecF2 rz = VecF20, string lx = "", string ly = "",  string lz = "") 
	{ drawToFile(sTitle, folder, name, step, DRAW_VECT, rx, ry, rz, lx, ly, lz, val1.m, val1.n, val1.arr(coeff), val2.arr(coeff), NULL, NULL); } 


	//static calls for images composed of more graphs
	void compSet(string folder, string name, int step, int m, int n)
	{ setName(folder, name, step); setSize(m, n); create(); }

	//PLOT
	void compPlot(string sTitle, VecL3 subplot, int m, DBL* val, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToBuff(sTitle, subplot, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", m, 0, val, NULL, NULL, NULL); } 

	void compPlot(string sTitle, VecL3 subplot, int m, DBL* val, DBL* pts, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToBuff(sTitle, subplot, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", m, 0, val, NULL, pts, NULL); } 

	void compPlot(string sTitle, VecL3 subplot, int m, DBL* val1, DBL* val2, DBL* pts, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToBuff(sTitle, subplot, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", m, 0, val1, val2, pts, NULL); } 

	void compPlot(string sTitle, VecL3 subplot, int m, int n, DBL* val1, DBL* val2,DBL* pts1, DBL* pts2, VecF2 rx = VecF20, VecF2 ry = VecF20, string lx = "", string ly = "") 
	{ drawToBuff(sTitle, subplot, DRAW_PLOT, rx, ry, VecF20, lx, ly, "", m, n, val1, val2, pts1, pts2); } 

	//SURF
	void compSurf(string sTitle, VecL3 subplot, int m, int n, DBL* val, VecF2 rx = VecF20, VecF2 ry = VecF20, VecF2 rz = VecF20, string lx = "", string ly = "",  string lz = "") 
	{ drawToBuff(sTitle, subplot, DRAW_SURF, rx, ry, rz, lx, ly, lz, m, n, val, NULL, NULL, NULL); } 

	//DENS
	void compDens(string sTitle, VecL3 subplot, int m, int n, DBL* val, VecF2 rx = VecF20, VecF2 ry = VecF20, VecF2 rz = VecF20, string lx = "", string ly = "",  string lz = "") 
	{ drawToBuff(sTitle, subplot, DRAW_DENS, rx, ry, rz, lx, ly, lz, m, n, val, NULL, NULL, NULL); } 

	//VECT
	void compVect(string sTitle, VecL3 subplot, int m, int n, DBL* val1, DBL* val2, VecF2 rx = VecF20, VecF2 ry = VecF20, VecF2 rz = VecF20, string lx = "", string ly = "",  string lz = "") 
	{ drawToBuff(sTitle, subplot, DRAW_VECT, rx, ry, rz, lx, ly, lz, m, n, val1, val2, NULL, NULL); } 

	
private:
	mglGraph *gr;
	mglData pts1, pts2; //for ticks  (domain)
	mglData val1, val2; //for values (codomain)
	string name;
	VecF2 rx, ry, rz;
	string lx, ly, lz;
	VecL2 size;
	int unitL;
	int type;
	int nLines;
	bool points;
	string title;

	static void drawToFile(string sTitle, string folder, string name, int step, int type, VecF2 rx, VecF2 ry, VecF2 rz, string labx, string laby, string labz, int m, int n, DBL* val1, DBL* val2, DBL* pts1, DBL* pts2)
	{
		Graph gr = Graph();
		gr.setTitle(sTitle);
		gr.setName(folder, name, step);
		gr.setType(type);
		gr.setRange(rx, ry, rz);
		gr.setLabel(labx, laby, labz);
		if(type == DRAW_DENS || type == DRAW_VECT) gr.setSize(m, n);
		gr.setData(m, n, val1, val2, pts1, pts2, true);
		gr.create();
		gr.put();	
		gr.draw();	
	}

	void drawToBuff(string sTitle, VecL3 subplot, int type, VecF2 rx, VecF2 ry, VecF2 rz, string labx, string laby, string labz, int m, int n, DBL* val1, DBL* val2, DBL* pts1, DBL* pts2)
	{
		setTitle(sTitle);
		setType(type);
		setRange(rx, ry, rz);
		setLabel(labx, laby, labz);
		if(type == DRAW_DENS || type == DRAW_VECT) setSize(m, n);
		setData(m, n, val1, val2, pts1, pts2, false);
		put(subplot);	
	}
};
*/

#endif
