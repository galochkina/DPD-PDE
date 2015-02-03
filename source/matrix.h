#ifndef MX_H
#define	MX_H

//MX interface
interface MXI
{
public:
	DBL& el(int row, int column);
	int dim();
};

//three-diagonal MX
class TriDiagM : MXI
{
public:
	TriDiagM()
	{
		n = 0;
	}

	TriDiagM(int dim)
	{
		n = dim;
		vec = vector<DBL>(dim*3-2, 0);
	}

	//TODO: test this if it really returns the reference on the asked element
	DBL& el(int row, int column)
	{
		int k = row - column;
		if(k == -1) return vec[row];
		else if(k == 0) return vec[n-1 + row];
		else if(k == 1) return vec[2*n-1 + column];
	}

	int dim()
	{
		return n;
	}

private:
	int n;			//dimension of the MX
	vector<DBL> vec; //elements occupying three non-trivial diagonals; first go elements from upper , then main, and at last, lower diagonal
};

//three-diagonal MX
class TDM
{
public:
	vector<DBL> a, b, c;

	TDM()
	{
		n = 0;
	}

	TDM(int dim)
	{
		n = dim;
		a = vector<DBL>(n,0);
		b = vector<DBL>(n,0);
		c = vector<DBL>(n,0);
	}

	int dim()
	{
		return n;
	}

private:
	int n;			//dimension of the MX
};

//three-diagonal MX
class MX : MXI
{
public:
	static const int NONE	= 0;	//any MX
	static const int D		= 1;	//diagonal MX
	static const int L		= 2;	//L from LU decomposition
	static const int U		= 3;	//U from LU decomposition
	static const int D3		= 4;	//tri-diagonal MX

	MX()
	{
		m = 0;
		n = 0;
	}

	MX(int rows, int columns)
	{
		m = rows;
		n = columns;
		vec = vector<DBL>(rows * columns, 0);
	}

	MX(const MX& A)
	{
		m = A.m;
		n = A.n;
		vec = vector<DBL>(m * n, 0);
		for(unsigned int i=0; i<vec.size(); i++) vec[i] = A.vec[i];
	}

	DBL& operator ()(const int row, const int column)
	{
		return vec[row * n + column];
	}

	const MX operator - ()                                            
	{
		MX B = MX(*this);
		for(unsigned int i=0; i<vec.size(); i++) B.vec[i] = -B.vec[i];
		return B;
	}

	const friend MX operator + (const MX& A, const MX& B)                                            
	{
		MX C = MX(A);
		for(unsigned int i=0; i<C.vec.size(); i++) C.vec[i] = C.vec[i] + B.vec[i];
		return C;
	}

	const friend MX operator - (const MX& A, const MX& B)                                            
	{
		MX C = MX(A);
		for(unsigned int i=0; i<C.vec.size(); i++) C.vec[i] = C.vec[i] - B.vec[i];
		return C;
	}

	const friend MX operator * (MX& A, MX& B)                                            
	{
		MX C = MX(A.m, B.n);
		for(int i=0; i<A.m; i++)
			for(int j=0; j<B.n; j++)
				for(int k=0; k<A.n; k++)
					C(i,j) += A(i,k) * B(k,j);
		return C;
	}

	int dimR()
	{
		return m;
	}

	int dimC()
	{
		return n;
	}

	int isNone() { return type == NONE; }
	int isD() { return type == D; }
	int isL() { return type == L; }
	int isU() { return type == U; }
	int isD3() { return type == D3; }

private:
	int m;				//dimension of the MX (rows)
	int n;				//dimension of the MX (columns)
	vector<DBL> vec;	//elements stored row by row
	int type;
};


#endif

