#ifndef OBJECTS_H
#define	OBJECTS_H

class Atom
{
public:
	static const int TypePlasma	= 0;
	static const int TypeSPlt	= 1;
	
	int		id;				//index of the atom in the attom array
	int		type;			//type of the atom (plasma, simple plt,...)
	VecF2	pos;			//position
	VecF2	v;				//velocity
	VecF2	F;				//force (DPD)
	bool	prep;			//in preparation area or not
	bool	fix;			//fixed position
	bool	free;
	bool	inClot;
	bool	copied;
	bool	oldPlt;
	VecF2	acc;
	VecF2	vOld;

	VecF2 Fprep;

	VecF2 FC;
	VecF2 FD1;
	VecF4 FD2;
	VecF2 FR;
	VecF2 FPL1M;
	VecF2 FPL1;
	VecF2 FPL2;
	VecF2 FPL3;

	Atom() : id(-1), type(TypePlasma), pos(0,0), v(0,0), F(0,0), prep(false), fix(false), free(false), 
		inClot(false), copied(false), oldPlt(false), acc(0,0), vOld(0,0),
		FC(VecF20), FD1(VecF20), FD2(VecF40), FR(VecF20), FPL1M(VecF20), FPL1(VecF20), FPL2(VecF20), FPL3(VecF20) {}

	Atom(int type, VecF2 pos, VecF2 v, bool prep, bool fix, bool inClot) 
		: id(-1), type(type), pos(pos), v(v), prep(prep), fix(fix), inClot(inClot), acc(VecF20), F(VecF20), free(false), copied(false), oldPlt(false), vOld(v), 
		FC(VecF20), FD1(VecF20), FD2(VecF40), FR(VecF20), FPL1M(VecF20), FPL1(VecF20), FPL2(VecF20), FPL3(VecF20) {}

	Atom(MyReadFile &f) : FC(VecF20), FD1(VecF20), FD2(VecF40), FR(VecF20), FPL1M(VecF20), FPL1(VecF20), FPL2(VecF20), FPL3(VecF20)
	{
		f.read(id);
		f.read(type);
		f.read(pos.x);
		f.read(pos.y);
		f.read(v.x);
		f.read(v.y);
		f.read(prep);
		f.read(fix);
		f.read(free);
		f.read(inClot);
		f.read(copied);
		f.read(oldPlt);
		f.read(acc.x);
		f.read(acc.y);
		F = VecF20;
		vOld = VecF20;
	}

	Atom(const Atom &a)
	{
		id		= a.id;
		type	= a.type;
		pos		= a.pos;
		v		= a.v;
		F		= a.F;

		FC = VecF20;
		FD1 = VecF20;
		FD2 = VecF40;
		FR = VecF20;
		FPL1 = VecF20;
		FPL1M = VecF20;
		FPL2 = VecF20;
		FPL3 = VecF20;

		prep	= a.prep;
		fix		= a.fix;
		free	= a.free;
		inClot  = a.inClot;
		copied	= a.copied;
		oldPlt	= a.oldPlt;
		acc		= a.acc;
		vOld	= a.vOld;
	}

	inline void resetF() 
	{ 
		F = VecF20;  

		FC = VecF20;
		FD1 = VecF20;
		FD2 = VecF40;
		FR = VecF20;
		FPL1 = VecF20;
		FPL1M = VecF20;
		FPL2 = VecF20;
		FPL3 = VecF20;
	}

	inline void resetFdpd() 
	{ 
		F = VecF20; 
				
		FC = VecF20;
		FD1 = VecF20;
		FD2 = VecF40;
		FR = VecF20;
		FPL1M = VecF20;
	}

	inline void resetFplt() 
	{ 
		FPL1 = VecF20;
		FPL2 = VecF20;
		FPL3 = VecF20;	
	}

	inline bool isPlasma() 
	{
		return type == TypePlasma;
	}

	inline bool isSPlt() 
	{
		return type == TypeSPlt;
	}



	void save(MyWriteFile &f)
	{		
		f.write(id);
		f.write(type);
		f.write(pos.x);
		f.write(pos.y);
		f.write(v.x); 
		f.write(v.y);
		f.write(prep);
		f.write(fix);
		f.write(free);
		f.write(inClot);
		f.write(copied);
		f.write(oldPlt);
		f.write(acc.x);
		f.write(acc.y);
	}
};

typedef vector<Atom*> LAP;
typedef vector<Atom> LA;

//*******************************************************************************************************************************************************
//SmartList
template <class T>
class SmartList
{
public:
	vector<T> obj; //to

	SmartList() {};

	SmartList(int objCount, DBL reserve)
	{
		obj = vector<T>(0);
		newObj = vector<T>(0);
		nRes = objCount * (1. + reserve);
		obj.reserve(nRes);

		objCount = 0;
		freeCount = 0;
		free = list<int>();
	}

	SmartList(int reserve)
	{
		obj = vector<T>(0);
		newObj = vector<T>(0);
		nRes = reserve;
		obj.reserve(nRes);

		objCount = 0;
		freeCount = 0;
		free = list<int>();
	}


	//puts object to temporary list and updates the main list after call of method update (this combination is used for parallel processing)
	void addLater(const T &a)
	{
		#pragma omp critical
		newObj.push_back(a);
	}

	//add objects from temp. list (newObj) to the main object list (obj)
	void update()
	{
		for(vector<T>::iterator i = newObj.begin(); i != newObj.end(); i++)
			addNow(*i);
		newObj.clear();
	}

	//adds the new object directly to the main list
	int addNow(const T &a)
	{
		int k = -1;
		if(freeCount <= 0)
		{
			k = obj.size();
			obj.push_back(a);
			obj[k].id = k;
			obj[k].free = false;
		}
		else
		{
			k = *(free.begin());
			obj[k] = a;
			obj[k].id = k;
			obj[k].free = false;
			free.pop_front();
			freeCount--;
		}
		objCount++;
		return k;
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::SLIS);
		f.write(getSize());
		f.write(count());
		f.newLine();

		for(vector<T>::iterator p=obj.begin(); p!=obj.end(); p++)
		{
			f.writeTag(Const::GOBJ);
			p->save(f);
			f.newLine();
		}
	}

	void load(MyReadFile &f)
	{
		clear();
		int objSize = 0;
			
		if(f.isTag(Const::SLIS)) 
		{
			f.read(objSize);
			int count = f.getInt();
			obj.reserve(objSize);
			f.nextLine();	

			//Objects
			for(int k=0; k<objSize; k++)
			{
				if(f.isTag(Const::GOBJ)) 
				{
					load(T(f));
				}
				f.nextLine();	
			}
		}
	}

	void remove(int index)
	{
		if(index < (int)obj.size() && !obj[index].free)
		{
			obj[index].free = true;
			#pragma omp critical
			{
				free.push_front(index);
				freeCount++;
				objCount--;
			}
		}
	}

	int count()
	{
		return objCount;
	}

	int getSize()
	{
		return obj.size();
	}

	void clear()
	{
		obj.clear();
		newObj.clear();
		free.clear();
		objCount = 0;
		freeCount = 0;
	}

private:
	LA newObj;
	int freeCount;
	list<int> free;
	int objCount;
	int nRes;	

	void load(T a)
	{
		obj.push_back(a);
		if(a.free) 
		{
			freeCount++;
			free.push_back(a.id);
		}
		else objCount++;
	}
};

//*******************************************************************************************************************************************************
typedef SmartList<Atom> AtomList;

#endif