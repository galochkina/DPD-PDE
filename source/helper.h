#ifndef HELPER_H
#define	HELPER_H

//*******************************************************************************************************************************************************
template <class T>
inline std::string toString (const T& t, int prec = 0, bool sci = false)
{
	std::stringstream ss;
	if(prec >= 0) ss.precision(prec);
	if(sci) ss << scientific;
	ss << t << flush;
	return ss.str();
}

//*******************************************************************************************************************************************************
//HELPER METHOD
//m mod n
inline DBL mod (DBL x, DBL n)
{
	DBL ret = fmod(x, n);
	if(ret < 0) ret += n;
	return ret;
}

inline int modInt (int x, int n)
{
	int ret = x % n;
	if(ret < 0) ret += n;
	return ret;
}
//*******************************************************************************************************************************************************
//FREE VECTOR MEMORY
template <typename T>
static void free(vector<T> &t) 
{
	t.clear();
    vector<T> tmp = vector<T>();
    t.swap(tmp);
}

template <typename T>
static void free(T *t) 
{
	if(t!=NULL) delete[] t;
	t = NULL;
}

//*******************************************************************************************************************************************************
//HELPER METHOD FOR DIALOGS
//conversion from bool to int and vice versa
inline int boolToInt(bool b)
{
	if(b) return 1;
	else return 0;
}

inline bool intToBool(int i)
{
	if(i > 0) return true;
	else return false;
}

//*******************************************************************************************************************************************************
class Stat
{
public:
	DBL n;
	VecF2 v;

	Stat() : n(0), v(VecF20) {};

	inline void update(VecF2 vel)
	{
		v += vel;
		n++;
	}

	inline void avg(DBL coeff = 1)
	{
		v = (n!=0) ? v * coeff / n : VecF20;
	}	
};


//*******************************************************************************************************************************************************
//splits given string by any char in delimeter and returns vector of splited parts
inline vector<string> splitString(string str, string delim)
{
	vector<string> results = vector<string>();
	unsigned int cutAt;
	while( (cutAt = str.find_first_of(delim)) != str.npos )
	{
		if(cutAt > 0)
		{
			results.push_back(str.substr(0,cutAt));
		}
		str = str.substr(cutAt+1);
	}
	if(str.length() > 0)
		results.push_back(str);

	return results;
}

//*******************************************************************************************************************************************************
//HELPER METHOD FOR IO
inline static vector<string> parseLine(string &line, string separator = ";")
{
	return splitString(line, separator);
}

class MyReadFile
{
public:
	MyReadFile(string &fileName, bool withTagsFlag = true)
	{
		withTags = withTagsFlag;
		sep = ";";
		file = vector<string>();
		ifstream f(fileName.c_str());
		if(f.is_open())
		{
			string line = "";
			while(getline(f, line))
				file.push_back(line);
			f.close();
			readable = true;

			iFile = 0;
			fileSize = file.size();

			nextLine();
		}
		else readable = false;
	}

	inline bool isReadable()	{ return readable; }
	inline bool nextLine()		{ if(iFile >= fileSize) { line.clear(); return false; } line = parseLine(file[iFile], sep); iLine = 0; lineSize = line.size(); iFile++; if(withTags) read(curTag); return true; }
	
	inline void read(string &s)	{ s = check() ? line[iLine++] : ""; }
	inline void read(bool &b)	{ b = check() ? intToBool(atoi(line[iLine++].c_str())) : false; }
	inline void read(int &i)	{ i = check() ? atoi(line[iLine++].c_str()) : 0; }
	inline void read(unsigned int &i)	{ i = check() ? atoi(line[iLine++].c_str()) : 0; }
	inline void read(long &l)	{ l = check() ? atol(line[iLine++].c_str()) : 0; }
	inline void read(unsigned long &l)	{ l = check() ? atol(line[iLine++].c_str()) : 0; }
	inline void read(DBL &d)	{ d = check() ? atof(line[iLine++].c_str()) : 0; }
	void read(VecF2 &v)			{ v.x = check() ? atof(line[iLine++].c_str()) : 0; v.y = check() ? atof(line[iLine++].c_str()) : 0;	}
	void read(VecL2 &v)			{ v.x = check() ? atoi(line[iLine++].c_str()) : 0; v.y = check() ? atoi(line[iLine++].c_str()) : 0;	}
	void read(VecF3 &v)			{ v.x = check() ? atof(line[iLine++].c_str()) : 0; v.y = check() ? atof(line[iLine++].c_str()) : 0; v.z = check() ? atof(line[iLine++].c_str()) : 0; }
	void read(VecL3 &v)			{ v.x = check() ? atoi(line[iLine++].c_str()) : 0; v.y = check() ? atoi(line[iLine++].c_str()) : 0;	v.z = check() ? atoi(line[iLine++].c_str()) : 0; }
	void read(Stat &s)			{ s.n = check() ? atoi(line[iLine++].c_str()) : 0; s.v.x = check() ? atoi(line[iLine++].c_str()) : 0; s.v.y = check() ? atoi(line[iLine++].c_str()) : 0; }


	inline string getStr()		{ return check() ? line[iLine++] : ""; }
	inline bool getBool()		{ return check() ? intToBool(atoi(line[iLine++].c_str())) : false; }
	inline int getInt()			{ return check() ? atoi(line[iLine++].c_str()) : 0; }
	inline long getLong()		{ return check() ? atol(line[iLine++].c_str()) : 0; }
	inline DBL getDbl()			{ return check() ? atof(line[iLine++].c_str()) : 0; }

	inline bool isTag(const string &s) { if(!withTags) return false; else return s.compare(curTag) == 0; }
	inline void setWithTags(bool tf) { withTags = tf; }

private:
	vector<string> file;
	vector<string> line;
	string curTag;
	unsigned int iFile;
	unsigned int iLine;
	unsigned int lineSize;
	unsigned int fileSize;
	bool readable;
	string sep;
	bool withTags;

	bool check()
	{
		return iLine < line.size();
	}
};

class MyWriteFile
{
public:
	MyWriteFile(string &fileName, bool append, int width = 0, bool autoFormat = false) : w(width), autoFormat(autoFormat)
	{
		sep = ";";
		if(append) os.open(fileName.c_str(), ios_base::app);
		else os.open(fileName.c_str());
		os << setprecision(6) << flush;
	}

	inline void writeTag(const string &tag)	{ os << tag; }
	inline void write(const string &s)	{ os << sep << setw(w) << s; }
	inline void write(const char* c)	{ string str(c); os << sep << setw(w) << str; }
	inline void write(const int i)		{ if(autoFormat) writeNum(i); else os << sep << setw(w) << scientific << i; }
	inline void write(const long l)		{  if(autoFormat) writeNum(l); else os << sep << setw(w) << scientific << l; }
	inline void write(const unsigned int i)	{  if(autoFormat) writeNum(i); else os << sep << setw(w) << scientific << i; }
	inline void write(const unsigned long l) {  if(autoFormat) writeNum(l); else os << sep << setw(w) << scientific << l; }
	inline void write(const DBL d)		{  if(autoFormat) writeNum(d); else os << sep << setw(w) << scientific << d; }
	inline void write(const bool b)		{ os << sep << setw(w) << b; }
	inline void write(const VecF2 &v)	{ if(autoFormat){ writeNum(v.x); writeNum(v.y); } else {os << sep << setw(w) << scientific << v.x; os << sep << setw(w) << scientific << v.y;} }
	inline void write(const VecL2 &v)	{ if(autoFormat){ writeNum(v.x); writeNum(v.y); } else {os << sep << setw(w) << scientific << v.x; os << sep << setw(w) << scientific << v.y;} }
	inline void write(const VecF3 &v)	{ if(autoFormat){ writeNum(v.x); writeNum(v.y); writeNum(v.z); } else {os << sep << setw(w) << scientific << v.x; os << sep << setw(w) << scientific << v.y; os << sep << setw(w) << scientific << v.z;} }
	inline void write(const VecL3 &v)	{ if(autoFormat){ writeNum(v.x); writeNum(v.y); writeNum(v.z); } else {os << sep << setw(w) << scientific << v.x; os << sep << setw(w) << scientific << v.y; os << sep << setw(w) << scientific << v.z;} }
	inline void write(const Stat &s)	{ if(autoFormat){ writeNum(s.n); writeNum(s.v.x); writeNum(s.v.y); } else {os << sep << setw(w) << scientific << s.n; os << sep << setw(w) << scientific << s.v.x; os << sep << setw(w) << scientific << s.v.y;} }
	inline void newLine()				{ os << sep << endl; }
	inline void close()					{ os.close(); }
	inline void setSep(string separator)	{ sep = separator; }

private:
	ofstream os;
	int w;
	string sep;
	bool autoFormat;

	void writeNum(DBL f)
	{
		if(!autoFormat) 
		{
			os << sep << setw(w) << scientific << f;
			return;
		}

		os << sep << setw(w);
		if(f == 0)
		{
			os << setprecision(0) << fixed << f;
			return;
		}

		DBL d = fabs(f);
		int p = 6;
		stringstream ss;
		ss << setprecision(p) << scientific << d << flush;
		string s = ss.str();
		string tmp;
		int nAbs = atoi(s.substr(10, 3).c_str());
		tmp = s.at(9);
		int n = (tmp.compare("-")==0) ? -nAbs : nAbs;

		//0.000000e+000
		int i=0;
		for(i=7; i>1; i--) 
		{
			tmp = s.at(i);
			if(atoi(tmp.c_str()) != 0) break;
		}
		p=i-1;

		if(n==nAbs)
		{
			if(p<=3 && n<=3)
			{
				p = (p<=n) ? 0 : p - n;
				os << fixed;
			}
			else 
			{
				if(p==0) p=1;
				os << scientific;
			}
		}
		else
		{
			if(p+nAbs<=4) 
			{
				p = p+nAbs;
				os << fixed;
			}
			else 
			{
				if(p==0) p=1;
				os << scientific;
			}
		}		

		os << setprecision(p) << f;
	}
};

//*******************************************************************************************************************************************************
template <class T>
static void saveVector(vector<T> &v, MyWriteFile &f)
{
	f.write(v.size());
	for(int i=0; i<v.size(); i++)
		f.write(v[i]);
}

template <class T>
static void loadVector(vector<T> &v, MyReadFile &f)
{
	free(v);
	int n = f.getInt();
	v = vector<T>(n);
	for(int i=0; i<v.size(); i++)
		f.read(v[i]);
}

//*******************************************************************************************************************************************************
class Clock
{
public:
	//time
	clock_t beg;
	//average over time for more start/stop-s
	clock_t sumTime;
	DBL avgRate;
	DBL sumCounter;


	Clock() : beg(clock()), a(0), b(0), counter(0), time(0), rate(0), sumTime(0), avgRate(0), sumCounter(0) {};

	void start() 
	{ 
		a = clock(); 
		counter = 0; 
	}
	
	void stop() 
	{ 
		b = clock(); 
		time = b-a; 
		rate = counter / (time / (DBL) CLOCKS_PER_SEC); 
		sumTime += time; 
		sumCounter += counter; 
		avgRate = sumCounter / (sumTime / (DBL) CLOCKS_PER_SEC);
	}

	void addCounter() { counter++; }

	vector<int> getTime() { return formatTime(b - beg); }
	DBL getRate() { return rate; }
	DBL getAvgRate() { return avgRate; }

	vector<int> getRemainingTime(int steps) 
	{ 
		clock_t elapsed = b - beg;
		clock_t tr = (clock_t)(elapsed * (steps / (DBL)sumCounter - 1));
		return formatTime(tr);
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::CLCK);
		f.write(beg);
		f.write(sumTime);
		f.write(avgRate);
		f.write(sumCounter);
		f.newLine();
	}

	void load(MyReadFile &f)
	{
		if(!f.isTag(Const::CLCK)) return;
		f.read(beg);
		f.read(sumTime);
		f.read(avgRate);
		f.read(sumCounter);
		f.nextLine();
	}
	
private:
	//time
	clock_t a; 
	clock_t b;
	//for one start/stop
	int counter;
	clock_t time;
	DBL rate;	

	vector<int> formatTime(clock_t c)
	{
		//h:m:s
		vector<int> t = vector<int>(3);
		long s = c / CLOCKS_PER_SEC;
		long m = s / 60;
		t[2] = s % 60;
		t[1] = m % 60;
		t[0] = m / 60;
		return t;
	}
};

//*******************************************************************************************************************************************************
template <class T> class Mat;

template <class T>
class Vec
{
public:		
	int n;
	vector<T> v;

	Vec()
	{
		n = 0;
		v = vector<T>();
	}
	
	Vec(int size)
	{
		n = size;
		v = vector<T>(n);
	}

	Vec(int size, T val)
	{
		n = size;
		v = vector<T>(n, val);
	}

	Vec(vector<T> &vec)
	{
		n = vec.size();
		v = vector<T>(n);
		for(int i=0; i<n; i++)
			v[i] = vec[i];
	}

	Vec(MyReadFile &f)
	{
		f.read(n);
		v = vector<T>(n);
		for(int i=0; i<n; i++) f.read(v[i]);
	}

	void pushBack(T &t)
	{
		v.push_back(t);
		n++;
	}

	T& operator ()(const int i)
	{
		return v[i];
	}

	void init(T val)
	{
		for(int i=0; i<n; i++)
			v[i] = val;
	}

	T* arr(DBL coeff = 1)
	{
		T* arr = new T[n];
		for(int i=0; i<n; i++)
			arr[i] = v[i] * coeff;
		return arr;
	}

	void multiply(DBL k)
	{
		for(int i=0; i<n; i++)
			v[i] *= k; 
	}

	void copy(T* arr, int size, DBL coeff = 1)
	{
		n = size;
		free(v);
		v = vector<T>(n);
		for(int i=0; i<n; i++)
			v[i] = arr[i] * coeff;
	}

	void copy(vector<T> &vec, DBL coeff = 1)
	{
		n = vec.size();
		free(v);
		v = vector<T>(n);
		for(int i=0; i<n; i++)
			v[i] = vec[i] * coeff;
	}

	void copy(Vec<T> &V, DBL coeff = 1)
	{
		n = V.n;
		free(v);
		v = vector<T>(n);
		for(int i=0; i<n; i++)
			v[i] = V(i) * coeff;
	}

	void copyMat(Mat<T> &M, DBL coeff = 1)
	{
		n = M.m * M.n;
		free(v);
		v = vector<T>(n);
		for(int i=0; i<n; i++)
			v[i] = M.v[i] * coeff;
	}


	VecF2 range(bool add = false)
	{
		if(add)
		{
			VecF2 r = getRange(v);
			DBL t = (r.y - r.x) / 10;
			r.x -= t;
			r.y += t;
			return r;
		}
		return getRange(v);
	}

	void save(MyWriteFile &f)
	{
		f.write(n);
		for(int i=0; i<n; i++) f.write(v[i]);
	}

	void load(MyReadFile &f)
	{
		free(v);
		f.read(n);
		v = vector<T>(n);
		for(int i=0; i<n; i++) f.read(v[i]);
	}

	void clear()
	{
		n = 0;
		free(v);
	}
	inline int size()
	{
		return n;
	}

	inline bool empty()
	{
		return n == 0;
	}
};

template <class T>
class Mat
{
public:	

	int m;
	int n;
	vector<T> v;

	Mat()
	{
		m = 0;
		n = 0;
		v = vector<T>();
	}
	
	Mat(int sizeX, int sizeY)
	{
		m = sizeX;
		n = sizeY;
		v = vector<T>(m*n);
	}

	Mat(MyReadFile &f)
	{
		f.read(m);
		f.read(n);
		v = vector<T>(m*n);
		for(int i=0; i<m*n; i++) f.read(v[i]);
	}

	T& operator ()(const int i, const int j)
	{
		return v[j*m+i];
	}

	T& operator ()(const VecL2 pos)
	{
		return v[pos.y*m+pos.x];
	}

	void init(T val)
	{
		for(int i=0; i<m*n; i++)
			v[i] = val;
	}

	T* arr(DBL coeff = 1)
	{
		T* arr = new T[m*n];
		for(int i=0; i<m*n; i++)
			arr[i] = v[i] * coeff;
		return arr;
	}

	void copy(Mat<T> &M, DBL coeff = 1)
	{
		m = M.m;
		n = M.n;
		free(v);
		v = vector<T>(m*n);
		if(coeff!=1)
		{
			for(int i=0; i<m*n; i++)
				v[i] = M.v[i] * coeff;
		}
		else
		{
			for(int i=0; i<m*n; i++)
				v[i] = M.v[i];
		}
	}

	void copyVec(Vec<T> &V, int dm, int dn, bool byRows = true)
	{
		m = dm;
		n = dn;
		free(v);
		v = vector<T>(m*n);
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				if(byRows) v[j*m+i] = V(i);	//stores vector by filling rows
				else v[i*n+j] = V(i);
	}

	VecF2 range()
	{
		return getRange(v);
	}

	void save(MyWriteFile &f)
	{
		f.write(m);
		f.write(n);
		for(int i=0; i<m*n; i++) f.write(v[i]);
	}

	void load(MyReadFile &f)
	{
		f.read(m);
		f.read(n);
		free(v);
		v = vector<T>(m*n);
		for(int i=0; i<m*n; i++) f.read(v[i]);
	}

	void clear()
	{
		m = 0;
		n = 0;
		free(v);
	}

	//returns the index in array v for position (m,n) in the matrix
	inline int index(int i, int j)
	{
		return j*m+i;
	}
};

//*******************************************************************************************************************************************************
class Hist
{
public:
	long n;			// counter for averaging
	DBL sum;		// current sum
	bool autoInc;	// automatic increment option
	Vec<DBL> hist;	// history
	DBL minVal;
	DBL maxVal;

	Hist(bool autoIncrement = false)
	{
		n = 0;
		sum = 0;
		autoInc = autoIncrement;
		hist = Vec<DBL>();
		maxVal = 0;
		minVal = 0;
	}

	void update(DBL value)
	{
		sum += value;
		if(autoInc) n++;
	}

	void inc()
	{
		if(!autoInc) n++;
	}

	void store(bool reset = true)
	{
		DBL sum2 = (sum == 0) ? 0 : sum / n;
		hist.pushBack(sum2);
		if(hist.size() < 2) 
		{
			minVal = getLast();
			maxVal = getLast();
		}
		else
		{
			minVal = min(minVal, sum2);
			maxVal = max(maxVal, sum2);
		}
		if(reset)
		{
			sum = 0;
			n = 0;
		}
	}

	void reset()
	{
		hist.clear();
		maxVal = 0;
		sum = 0;
		n = 0;
	}

	DBL getLast(int k=1)
	{
		if(k<1 || hist.empty()) return 0;
		if(k == 1) return hist(hist.size()-1);
		else 
		{
			int m = hist.n;
			DBL sum = 0;
			if(k>m) k = m;
			for(int i = m-k; i<m; i++) sum += hist(i);
			return sum / k;
		}
	}

	VecF2 getRange(bool pad = false)
	{
		if(pad)
		{
			DBL d = (maxVal - minVal) * 0.1;
			return VecF2(minVal - d, maxVal + d);
		}
		else return VecF2(minVal, maxVal);
	}

	void save(MyWriteFile &f)
	{
		f.write(n);
		f.write(sum);
		f.write(autoInc);
		f.write(minVal);
		f.write(maxVal);
		hist.save(f);
	}

	void load(MyReadFile &f)
	{
		f.read(n);
		f.read(sum);
		f.read(autoInc);
		f.read(minVal);
		f.read(maxVal);
		hist.load(f);
	}
};

//*******************************************************************************************************************************************************
template <class T> static VecF2 getRange(T* v, int n)
{
	VecF2 r(v[0], v[0]);
	for(int i=1; i<n; i++)
	{
		r.x = min(v[i], r.x);
		r.y = max(v[i], r.y);
	}
	return r;
}

template <class T> static VecF2 getRange(vector<T> &v)
{
	if(v.size()==0) return VecF20;
	VecF2 r(v[0], v[0]);
	for(int i=1; i<v.size(); i++)
	{
		r.x = min(v[i], r.x);
		r.y = max(v[i], r.y);
	}
	return r;
}

//*******************************************************************************************************************************************************
//d-data, f-filter
static DBL conv(Vec<DBL> &d, Vec<DBL> &f, int pos)
{
	int k = f.n/2;
	DBL sum = 0;
	for(int i=0, j=pos-k; i<f.n; i++, j=pos-k+i)
		if(j>-1 && j<d.n) 
			sum += d(j) * f(i);
		else
		{
			if(j<0)
				sum -= d(-j) * f(i);
			else if(j>d.n-1)
			{
				int k = d.n-1 - j;
				sum -= d(d.n+k) * f(i);
			}
		}
	return sum;
}

static VecF2 convX(Mat<VecF2> &M, int x, int y, Vec<DBL> &f, bool mirrorEdge)
{
	VecF2 meanL = VecF20;
	VecF2 meanR = VecF20;
	int k = f.n/2;
	if(!mirrorEdge)
	{
		for(int i=0; i<k; i++)
		{
			meanL += M(i,y); 
			meanR += M(M.m-1-i,y); 
		}
		meanL/=k;
		meanR/=k;
	}
	VecF2 sum = VecF20;
	for(int i=0, j=x-k; i<f.n; i++, j=x-k+i)
	{
		if(j>-1 && j<M.m) 
			sum += M(j,y) * f(i);
		else
		{
			DBL sig = mirrorEdge ? -1 : 1;

			if(j<0)
				if(mirrorEdge)
					sum += sig * M(-j,y) * f(i);
				else
					sum += meanL * f(i);
			else if(j>M.m-1)
			{
				if(mirrorEdge)
				{
					int k = M.m-1 - j;
					sum += sig * M(M.m-1+k,y) * f(i);
				}
				else
					sum += meanR * f(i);
			}
		}
	}
	return sum;
}

static VecF2 convY(Mat<VecF2> &M, int x, int y, Vec<DBL> &f, bool mirrorEdge)
{
	VecF2 meanL = VecF20;
	VecF2 meanR = VecF20;
	int k = f.n/2;
	if(!mirrorEdge)
	{
		for(int i=0; i<k; i++)
		{
			meanL += M(x,i); 
			meanR += M(x,M.n-1-i); 
		}
		meanL/=k;
		meanR/=k;
	}
	VecF2 sum = VecF20;
	for(int i=0, j=y-k; i<f.n; i++, j=y-k+i)
	{
		if(j>-1 && j<M.n) 
			sum += M(x,j) * f(i);
		else
		{
			DBL sig = mirrorEdge ? -1 : 1;
			
			if(j<0)
				if(mirrorEdge)
					sum +=sig * M(x,-j) * f(i);
				else
					sum += meanR * f(i);
			else if(j>M.n-1)
			{
				if(mirrorEdge)
				{
					int k = M.n-1 - j;
					sum += sig * M(x,M.n-1+k) * f(i);
				}
				else
					sum += meanR * f(i);
			}
		}
	}
	return sum;
}

//d-data, f-filter
static Vec<DBL> convAll(Vec<DBL> &d, Vec<DBL> &f)
{
	Vec<DBL> r = Vec<DBL>(d.n);
	for(int i=0; i<d.n; i++)
		r(i) = conv(d, f, i);
	return r;
}

static void convAll(Mat<VecF2> &M, Vec<DBL> &f, CProgressCtrl &p)
{
	Mat<VecF2> R = Mat<VecF2>(M.m, M.n);
	//smoothing in x direction
	int m = R.m;
	int n = R.n;
	#pragma omp parallel for
	for(int j=0; j<n; j++)
	{
		for(int i=0; i<m; i++)
			R(i,j) = convX(M, i, j, f, false);
	}
	p.StepIt();
	M.copy(R);
	//smoothing in y direction
	#pragma omp parallel for
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<n; j++)
			R(i,j) = convY(M, i, j, f, true);
	}
	p.StepIt();
	M.copy(R);
}

static Vec<DBL> gaussFilt(int n, DBL sigma)
{
	int k = n/2;
	Vec<DBL> f = Vec<DBL>(2*k+1);
	DBL sum = 0;
	//filter construction
	for(int i=0, j=-k; i<f.n; i++, j=i-k)
	{
		f(i) = exp((-pow(j,2.)/(2*pow(sigma,2))));
		sum += f(i);
	}
	//filter normalization
	for(int i=0; i<n; i++)
		f(i) /= sum;

	return f;
}

static inline DBL bilin(DBL x, DBL y, DBL x1, DBL x2, DBL y1, DBL y2, DBL Q11, DBL Q12, DBL Q21, DBL Q22)
{
	return ( (Q11 * (y2-y) + Q12 * (y-y1)) * (x2-x) + (Q21 * (y2-y) + Q22 * (y-y1)) * (x-x1) ) / ((x2-x1)*(y2-y1));
}

static inline DBL circlePartMult(DBL r, DBL x, DBL y)
{
	DBL r2 = r*r;
	DBL K = x*y;
	DBL A = r2*acos(y/r)-y*sqrt(r2-y*y);
	DBL B = r2*acos(x/r)-x*sqrt(r2-x*x);
	DBL C = r2*pi;
	DBL res = (C/4-A/2-K)/B + 0.5;
	return res;
}

//*******************************************************************************************************************************************************
static inline DBL dist(VecF2 &a, VecF2 &b)
{
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

static DBL angle(VecF2 &ap, VecF2 &a, VecF2 &an)
{
	VecF2 vp = UNIT(ap - a);
	VecF2 vn = UNIT(an - a);

	DBL vpc = acos(vp.x);
	DBL vps = asin(vp.y);
	DBL vpa = (vps < 0) ? 2.*pi-vpc : vpc;

	DBL vnc = acos(vn.x);
	DBL vns = asin(vn.y);
	DBL vna = (vns < 0) ? 2.*pi-vnc : vnc;

	return (vpa > vna) ? vpa-vna : 2.*pi+vpa-vna;
}

#endif
