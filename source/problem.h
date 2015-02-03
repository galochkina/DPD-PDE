#ifndef PROBLEM_H
#define	PROBLEM_H


//*******************************************************************************************************************************************************
//PlateletPair
class PP
{
public:
	int a;
	int b;
	DBL d;
	DBL t;
	bool del;

	PP() : a(-1), b(-1), d(1) { t = 0; del = false; }

	PP(Atom &a1, Atom &a2, DBL dist)
	{
		a = a1.id;
		b = a2.id;
		d = dist;
		del = false;
		t = 0;
	}

	PP(MyReadFile &f)
	{
		if(f.isTag(Const::PPAR)) 
		{
			f.read(a);
			f.read(b);
			f.read(d);
			f.read(t);
			del = false;
		}
		f.nextLine();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PPAR);
		f.write(a);
		f.write(b);
		f.write(d);
		f.write(t);
		f.newLine();
	}
};


//*******************************************************************************************************************************************************
//PlateletPairList
class PPL
{
public:
	vector<PP> l;

	PPL()
	{
		l = vector<PP>();
	}

	void add(Atom &a1, Atom &a2, DBL d)
	{
		PP conn = (a1.id < a2.id) ? PP(a1,a2,d) : PP(a2,a1,d);
		vector<PP>::iterator i = l.begin();

		while(i != l.end() && i->a < conn.a) i++;
		if(i != l.end() && i->a == conn.a)
		{
			while(i != l.end() && i->a == conn.a && i->b < conn.b) i++;
			if(i != l.end() && i->a == conn.a && i->b == conn.b) return;
		}
		l.insert(i, conn);
		if(a1.inClot || a2.inClot) 
		{
			a1.inClot = true;
			a2.inClot = true;
		}
	}

	int find(int x, int y)
	{
		if(x == y) return -1;
		int a = x; int b = y;
		if(y < x)
		{
			a = y;
			b = x;
		}
		for(unsigned int i=0; i<l.size(); i++)
		{
			if(l[i].a == a && l[i].b == b) return i;
			i++;
		}
		return -1;
	}

	bool contains(int x, int y)
	{
		if(find(x, y) > 0) return true;
		else return false;
	}

	void remove(int i)
	{
		l[i].del = true;
	}

	void refresh()
	{
		vector<PP> lNew = vector<PP>();

		for(unsigned int i=0; i<l.size(); i++)
			if(!l[i].del) lNew.push_back(l[i]);

		l = lNew;
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PCON);
		f.write(l.size());
		f.newLine();

		for(vector<PP>::iterator p=l.begin(); p!=l.end(); p++)
			p->save(f);
	}

	void load(MyReadFile &f)
	{
		free(l);
		int lSize = 0;
		if(f.isTag(Const::PCON)) 
		{
			f.read(lSize);
			l.reserve(lSize);
		}
		f.nextLine();	

		//PP
		for(int k=0; k<lSize; k++)
			l.push_back(PP(f));
	}
};

//*******************************************************************************************************************************************************
//boundary conditions
class BoundCond
{
public:
	bool	periodicTB;			//top and bottom periodic boundary condition
	bool	periodicLR;			//left and right periodic boundary condition

	BoundCond(){};

	BoundCond(Setup &s)
	{
		periodicTB = false;
		//periodicLR = !s.spc.prepArea;
		periodicLR = true;
	}

private:
	
};

//*******************************************************************************************************************************************************
class Box
{
public:
	LAP atoms;
	VecL2 pos;

	Box() {};

	Box(int x, int y) 
	{
		pos = VecL2(x,y);
	}

	void add(Atom *a)
	{
		//#pragma omp critical
		{
			atoms.push_back(a);
		}
	}

	void clear()
	{
		atoms.clear();
	}

private:
	
};

//*******************************************************************************************************************************************************
class BoxPair
{
public:
	Box *A;
	Box *B;
	VecF2 displacement;
	bool singleBox;
	bool oneWay;

	BoxPair() : A(NULL), B(NULL), displacement(VecF20), singleBox(true), oneWay(false) {}

	BoxPair(Box *A) : A(A), B(A), displacement(VecF20), singleBox(true), oneWay(false) {}

	BoxPair(Box *A, Box *B, VecF2 displacement, bool oneWay) : A(A), B(B), displacement(displacement), singleBox(false), oneWay(oneWay) {}

private:
	
};

//*******************************************************************************************************************************************************
typedef vector<BoxPair> LBP;
typedef vector<LBP>		LLBP;

//*******************************************************************************************************************************************************
class Mesh
{
public:
	
	LBP boxPairs;
	LLBP boxPairsList;
	//default constructor: not for use
	Mesh(){ boxes = NULL; }

	Mesh(Setup &s, BoundCond &bc) 
	{
		prepArea = s.spc.prepArea;
		if(prepArea)
		{
			spaceSizeXP = s.spc.sizePrepX;
			mSizeXP = spaceSizeXP / s.spc.maxR;
			intXP = spaceSizeXP / mSizeXP;
		}
		else
		{
			spaceSizeXP = 0;
			mSizeXP = 0;
			intXP = 0;
		}

		spaceSize = VecF2(s.spc.size.x, s.spc.size.y);
		mSizeXS = (spaceSize.x - spaceSizeXP) / s.spc.maxR;
		intXS = (spaceSize.x - spaceSizeXP) / mSizeXS;

		mSize = VecL2(mSizeXP + mSizeXS, spaceSize.y / s.spc.maxR);
		intY = spaceSize.y / mSize.y;

		boxes = new Box[mSize.x*mSize.y];
		initBoxes();
		proc = s.opt.proc;
		setBoxPairs(bc);

		printBoxPairs();
	}

	void printBoxPairs()
	{
		string str = "boxes.txt";
		MyWriteFile f(str, false);
		for(unsigned int i = 0; i<boxPairsList.size(); i++)
		{
			f.write("ARRAY");
			f.write(i);
			f.newLine();
			for(unsigned int j = 0; j<boxPairsList[i].size(); j++)
			{
				f.write(boxPairsList[i][j].A->pos.x);
				f.write(boxPairsList[i][j].A->pos.y);
				f.write("   ");
				f.write(boxPairsList[i][j].B->pos.x);
				f.write(boxPairsList[i][j].B->pos.y);
				if(boxPairsList[i][j].oneWay) f.write(" ***");
				f.newLine();
			}
		}
		f.close();
	}

	void setParticles(LA &atoms)
	{
		clearBoxes();

		int n = atoms.size();
		//#pragma omp parallel for
		for(int i=0; i<n; i++)
		{
			if(!atoms[i].free)
			{
				Atom *a = &(atoms[i]);
				VecL2 c = getBoxCoords(a->pos);
				add(c, a);
			}
		}
	}

	//adds the particle pointer to the list of paricles which is element of matrix given by its coordinates
	void add(VecL2 pos, Atom *a)
	{
		boxes[pos.x * mSize.y + pos.y].add(a);
	}

	void clearBoxes()
	{		
		int n = mSize.x * mSize.y;
		for(int i=0; i<n; i++)
			boxes[i].clear();
	}

	void destroy() 
	{ 
		if(boxes != NULL) delete [] boxes;
		boxes = NULL;
		free(boxPairs); 
	}

	//returns the number of columns of the matrix
	int sizeX() { return mSize.x; }
	//returns the number of rows of the matrix
	int sizeY() { return mSize.y; }
	//returns the number of prep area rows of the matrix
	DBL getIntXP() { return intXP; }
	//returns the number of rows of the matrix
	DBL getIntXS() { return intXS; }
	//returns the number of columns of the matrix
	DBL getIntY() { return intY; }

private:
	bool	prepArea;		//is there a preparation area
	int		mSizeXP;		//mesh number of boxes in x direction for preparation area
	int		mSizeXS;		//mesh number of boxes in x direction for simulation area
	VecL2	mSize;			//mesh number of boxes	
	Box*	boxes;			//pointer to an array of pointers to lists
	DBL		intXP;			//interval for a box in x direction for prep. area
	DBL		intXS;			//interval for a box in x direction for sim. area
	DBL		intY;			//interval for a box in y direction for sim. area
	DBL		spaceSizeXP;	//space size.x for prep.area
	VecF2	spaceSize;		//space size
	int		proc;

	VecL2 getBoxCoords(VecF2 &pos)
	{
		int x = 0;
		if(pos.x < spaceSizeXP)
			x = pos.x / intXP;
		else
			x = mSizeXP + (pos.x-spaceSizeXP) / intXS;
		if(x == mSize.x) x--;
		int y = pos.y / intY;
		if(y == mSize.y) y--;
		return VecL2(x, y);
	}

	//returns pointer to element of matrix given by its coordinates
	Box* get(int x, int y)
	{
		return &(boxes[x * mSize.y + y]);
	}

	//init boxes positions
	void initBoxes()
	{
		for(int x=0; x<mSize.x; x++)
			for(int y=0; y<mSize.y; y++)
				get(x,y)->pos = VecL2(x,y);
	}

	//sets the list of pair of boxes that are to be calculated
	void setBoxPairs(BoundCond &bc)
	{
		//SET d
		int dy = mSize.y / proc;
		if(mSize.y % proc > 0) dy++;
		int dn = mSize.y / dy;
		if(mSize.y % dy > 0) dn++;
		
		//ALLOCATE VECTORS
		boxPairsList = vector<LBP>(dn * 2);
		for(int j=0; j<mSize.y; j++)
			for(int i=0; i<mSize.x; i++)
				setPointUp(dy, i, j, bc.periodicLR);
		
		for(int j=0; j<mSize.y; j+=dy)
			for(int i=0; i<mSize.x; i++)
				setPointDown(dn, dy, i, j, bc.periodicLR, bc.periodicTB);
	}
	/*
	void setPointUp(int dy, int x, int y, bool xPeriodicBC)
	{
		int n = y / dy;		

		bool up		= (y+1)/dy == n && y+1 < mSize.y;
		bool left	= x-1 >= 0 || xPeriodicBC || prepArea;
		bool right	= x+1 < mSize.x || xPeriodicBC;

		bool leftPA_0	= prepArea && x == 0;
		bool leftPA_XP	= prepArea && x == mSizeXP;
		bool rightPA	= prepArea && x == mSizeXP-1;

		// o
		boxPairsList[n].push_back(BoxPair(get(x,y)));
		
		VecF2 dis = VecF20;

		// o.
		if(right) 
		{
			if(rightPA)
			{
				dis.x = 0;
				boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), y), dis, true));
				dis.x = spaceSizeXP;
				boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSizeXP), y), dis, false));
			}
			else
			{
				dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
				boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), y), dis, false));
			}		
		}

		if(up)
		{
			//. 
			// o
			if(left) 
			{
				if(leftPA_0)
				{
					dis.x = -spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSizeXP), modInt(y+1, mSize.y)), dis, false));
				}
				else if(leftPA_XP)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(modInt(x-1, mSize.x), modInt(y+1, mSize.y)), get(x,y), dis, true));
				}
				else 
				{
					dis.x = (x-1 < 0) ? -spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSize.x), modInt(y+1, mSize.y)), dis, false));
				}	
			}

			// .
			// o
			dis.x = 0;
			boxPairsList[n].push_back(BoxPair(get(x,y), get(x,  modInt(y+1, mSize.y)), dis, false));
			
			//  .
			// o 
			if(right) 
			{
				if(rightPA)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y+1, mSize.y)), dis, true));
					dis.x = spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSizeXP), modInt(y+1, mSize.y)), dis, false));
				}
				else 
				{
					dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y+1, mSize.y)), dis, false));
				}	
			}
		}
	}

	void setPointDown(int dn, int dy, int x, int y, bool xPeriodicBC, bool yPeriodicBC)
	{
		int n = y / dy + dn;	

		bool down	= y-1 >= 0 || yPeriodicBC;
		if(!down) return;

		bool left	= x-1 >= 0 || xPeriodicBC || prepArea;
		bool right	= x+1 < mSize.x || xPeriodicBC;

		bool leftPA_0	= prepArea && x == 0;
		bool leftPA_XP	= prepArea && x == mSizeXP;
		bool rightPA	= prepArea && x == mSizeXP-1;
		
		VecF2 dis = VecF20;
		if(y-1 < 0) dis.y = -spaceSize.y;

		if(down)
		{
			// o 
			//.
			if(left) 
			{
				if(leftPA_0)
				{
					dis.x = -spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSizeXP), modInt(y-1, mSize.y)), dis, false));
				}
				else if(leftPA_XP)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(modInt(x-1, mSize.x), modInt(y-1, mSize.y)), get(x,y), dis, true));
				}
				else 
				{
					dis.x = (x-1 < 0) ? -spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSize.x), modInt(y-1, mSize.y)), dis, false));
				}					
			}

			// o
			// .
			dis.x = 0;
			boxPairsList[n].push_back(BoxPair(get(x,y), get(x, modInt(y-1, mSize.y)), dis, false));

			// o
			//  .
			if(right) 
			{
				if(rightPA)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y-1, mSize.y)), dis, true));
					dis.x = spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSizeXP), modInt(y-1, mSize.y)), dis, false));
				}
				else 
				{
					dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y-1, mSize.y)), dis, false));
				}
			}
		}
	}
	*/

	void setPointUp(int dy, int x, int y, bool xPeriodicBC)
	{
		int n = y / dy;		

		bool up		= (y+1)/dy == n && y+1 < mSize.y;
		bool left	= x-1 >= 0 || xPeriodicBC || prepArea;
		bool right	= x < mSize.x || xPeriodicBC;

		bool leftPA_0	= prepArea && x == 0;
		bool leftPA_XP	= prepArea && x == mSizeXP;
		bool rightPA	= prepArea && x == mSizeXP-1;
		bool rightPAT	= prepArea && x == mSize.x-1;

		// o
		boxPairsList[n].push_back(BoxPair(get(x,y)));
		
		VecF2 dis = VecF20;

		// o.
		if(right) 
		{
			if(rightPA)
			{
				dis.x = 0;
				boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), y), dis, true));
				dis.x = spaceSizeXP;
				boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSizeXP), y), dis, false));
			}
			else if(rightPAT)
			{
				dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
				boxPairsList[n].push_back(BoxPair(get(modInt(x+1, mSize.x), y), get(x,y), -dis, true));
			}
			else
			{
				dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
				boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), y), dis, false));
			}		
		}

		if(up)
		{
			//. 
			// o
			if(left) 
			{
				if(leftPA_0)
				{
					dis.x = -spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSizeXP), modInt(y+1, mSize.y)), dis, false));

					dis.x = -spaceSize.x;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSize.x), modInt(y+1, mSize.y)), dis, true));
				}
				else if(leftPA_XP)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(modInt(x-1, mSize.x), modInt(y+1, mSize.y)), get(x,y), dis, true));
				}
				else 
				{
					dis.x = (x-1 < 0) ? -spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSize.x), modInt(y+1, mSize.y)), dis, false));
				}	
			}

			// .
			// o
			dis.x = 0;
			boxPairsList[n].push_back(BoxPair(get(x,y), get(x,  modInt(y+1, mSize.y)), dis, false));
			
			//  .
			// o 
			if(right) 
			{
				if(rightPA)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y+1, mSize.y)), dis, true));
					dis.x = spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSizeXP), modInt(y+1, mSize.y)), dis, false));
				}
				else if(rightPAT)
				{
					dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(modInt(x+1, mSize.x), modInt(y+1, mSize.y)), get(x,y), -dis, true));
				}
				else 
				{
					dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y+1, mSize.y)), dis, false));
				}	
			}
		}
	}

	void setPointDown(int dn, int dy, int x, int y, bool xPeriodicBC, bool yPeriodicBC)
	{
		int n = y / dy + dn;	

		bool down	= y-1 >= 0 || yPeriodicBC;
		if(!down) return;

		bool left	= x-1 >= 0 || xPeriodicBC || prepArea;
		bool right	= x < mSize.x || xPeriodicBC;

		bool leftPA_0	= prepArea && x == 0;
		bool leftPA_XP	= prepArea && x == mSizeXP;
		bool rightPA	= prepArea && x == mSizeXP-1;
		bool rightPAT	= prepArea && x == mSize.x-1;
		
		VecF2 dis = VecF20;
		if(y-1 < 0) dis.y = -spaceSize.y;

		if(down)
		{
			// o 
			//.
			if(left) 
			{
				if(leftPA_0)
				{
					dis.x = -spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSizeXP), modInt(y-1, mSize.y)), dis, false));

					dis.x = -spaceSize.x;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSize.x), modInt(y-1, mSize.y)), dis, true));
				}
				else if(leftPA_XP)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(modInt(x-1, mSize.x), modInt(y-1, mSize.y)), get(x,y), dis, true));
				}
				else 
				{
					dis.x = (x-1 < 0) ? -spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x-1, mSize.x), modInt(y-1, mSize.y)), dis, false));
				}					
			}

			// o
			// .
			dis.x = 0;
			boxPairsList[n].push_back(BoxPair(get(x,y), get(x, modInt(y-1, mSize.y)), dis, false));

			// o
			//  .
			if(right) 
			{
				if(rightPA)
				{
					dis.x = 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y-1, mSize.y)), dis, true));
					dis.x = spaceSizeXP;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSizeXP), modInt(y-1, mSize.y)), dis, false));
				}
				else if(rightPAT)
				{
					dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(modInt(x+1, mSize.x), modInt(y-1, mSize.y)), get(x,y), -dis, true));
				}
				else 
				{
					dis.x = (x+1 >= mSize.x) ? spaceSize.x : 0;
					boxPairsList[n].push_back(BoxPair(get(x,y), get(modInt(x+1, mSize.x), modInt(y-1, mSize.y)), dis, false));
				}
			}
		}
	}
};


#endif
