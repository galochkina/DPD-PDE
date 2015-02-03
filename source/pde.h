#ifndef PDE_H
#define	PDE_H

//thrombin
// dT/dt = alpha (d^2T/dx^2) - V.(Del T) + f(T)(C0-T) - gamma T
// f(T) = k1 (T^2)/(T0+T)
class Pde1
{
public:
	DBL dx;
	DBL dy;
	DBL dt;
	Mat<DBL> *X;
	Mat<DBL> *Xnew;
	int m;
	int n;
	DBL alpha3;
	DBL k1;
	DBL gamma3;
	DBL beta1;
	DBL beta2;
	DBL beta3;
	DBL C0;
	DBL T0;
	Mat<DBL> B;
	DBL t2, t2x, at2x2, atx2, t2y, at2y2, aty2;
	VecF2 size;
	vector<DBL> posX;
	Mat<VecF2> V;
	Mat<DBL> lastX;
	VecF2 wound;

	Pde1()
	{
		X		= NULL;
		Xnew	= NULL;
		V		= Mat<VecF2>();
		lastX	= Mat<DBL>();
		B		= Mat<DBL>();
	}

	Pde1(Setup &s)
	{
		if(!s.fib.enabled) return;
		size = s.spc.size;
		alpha3 = s.fib.alpha3;
		k1 = s.fib.k1;
		gamma3 = s.fib.gamma3;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dx = s.fib.dx;
		dy = s.fib.dx;
		dt = s.fib.dt;
		wound = VecF2(s.fib.woundBeg, s.fib.woundEnd);
		
		m = size.x / dx + 1;
		n = size.y / dx + 1;
		X = new Mat<DBL>(m, n);
		B = Mat<DBL>(m, n);
		Xnew = new Mat<DBL>(m, n);
		Mat<DBL> B(m, n);
		V = Mat<VecF2>(m,n);
		lastX = Mat<DBL>(m,n);

		V.init(VecF20);
		//initial conditions
		X->init(0);
		Xnew->init(0);
		openedWound();

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha3 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha3 * t2y / dy;
		aty2 = at2y2 * 2;

		lastX.copy(*X);
	}

	void paramHotChange(Setup &s)
	{
		alpha3 = s.fib.alpha3;
		k1 = s.fib.k1;
		gamma3 = s.fib.gamma3;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dt = s.fib.dt;
		wound = VecF2(s.fib.woundBeg, s.fib.woundEnd);

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha3 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha3 * t2y / dy;
		aty2 = at2y2 * 2;
	}

	void nextStep(bool closedWound, Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		halfStepX(closedWound,  Tm, Tf);
		halfStepY(closedWound,  Tm, Tf);
	}

	void free()
	{
		if(X != NULL) 
		{
			delete X;
			X = NULL;
		}
		if(Xnew != NULL) 
		{
			delete Xnew;
			Xnew = NULL;
		}
	}

	void saveWave()
	{
		int y = n / 2;
		for(int i=0; i<m; i++)
		{
			if(((*X)(i,y)-0.5)*((*X)(i+1,y)-0.5) < 0)
			{
				DBL aaa = (*X)(i,y)-(*X)(i+1,y);
				posX.push_back(((*X)(i,y)/aaa-0.5/aaa+i)*dx);
				break;
			}
		}
	}


	//returns the concentration on given coordinates p
	inline DBL getU(VecF2 &p)
	{	
		int i = p.x / dx;
		int j = p.y / dy;
		return bilin(p.x, p.y, dx*i, dx*(i+1), dy*j, dy*(j+1), lastX(i, j), lastX(i+1, j), lastX(i, j+1), lastX(i+1, j+1));
		//return (lastX(i, j) + lastX(i+1, j) + lastX(i, j+1) + lastX(i+1, j+1)) / 4;
		
	}

	void update(DBL coeff)
	{
		for(int j=0; j<n; j++)
		{
			DBL y = dy * j;
			for(int i=0; i<m; i++)
			{
				V(i,j) = VecF2(coeff * y*(size.y - y) * 4. / pow(size.y, 2), 0);
				//V(i,j) = VecF2(coeff,0);
			}
		}

		lastX.copy(*X);
	}

	void update(VecF2 v)
	{
		for(int j=0; j<n; j++)
			for(int i=0; i<m; i++)
				V(i,j) = v;
		lastX.copy(*X);
	}

	void update(Mat<Stat> &M, DBL ignoreVal, CProgressCtrl &p, CStatic &desc)
	{
		DBL da = size.x / M.m;
		DBL db = size.y / M.n;
		DBL dah = da / 2;
		DBL dbh = db / 2;
		DBL daL = size.x - da;
		DBL dbL = size.y - db;

		int i1=0, i2=0, j1=0, j2=0;
		DBL dxL=0, dx1=0, dx2=0, dyL=0, dy1=0, dy2=0;

		desc.SetWindowTextA("Adjusting velocity profile for PDE.");
		for(int i=0; i<m; i++)
		{
			bool x0 = false;
			DBL x = dx * i - dah;

			if(x < 0)
			{
				x0 = true;
				i1 = 0;
			}
			else if(x >= daL)
			{
				x0 = true;
				i1 = M.m-1;
			}
			else
			{
				i1 = x / da;
				i2 = i1 + 1;
				dxL = da;
				dx1 = dx * i - (i1 * da + dah);
				dx2 = i2 * da + dah - dx * i;
			} 

			for(int j=0; j<n; j++)
			{
				DBL y = dy * j - dbh;

				if(y < 0)
				{
					dy1 = dy * j;
					if(x0) 
						V(i,j) = M(i1,0).v * dy1 / dbh;
					else 
						V(i,j) = ((M(i1,0).v * dx2 + M(i2,0).v * dx1)) * dy1 / (da * dbh);
				}
				else if(y >= dbL)
				{
					dy2 = size.y - dy * j;
					if(x0)
						V(i,j) = M(i1, M.n-1).v * dy2 / dbh;
					else
						V(i,j) = (M(i1, M.n-1).v * dx2 + M(i2, M.n-1).v * dx1) * dy2 / (da * dbh);
				}
				else
				{
					j1 = y / db;
					j2 = j1 + 1;
					dyL = db;
					dy1 = dy * j - (j1 * db + dbh);
					dy2 = j2 * db + dbh - dy * j;

					if(x0)
						V(i,j) = (M(i1, j1).v * dy2 + M(i1, j1).v * dy1) / db;
					else
						V(i,j) = (M(i1, j1).v * dx2 * dy2 
							+ M(i2, j1).v * dx1 * dy2
							+ M(i1, j2).v * dx2 * dy1
							+ M(i2, j2).v * dx1 * dy1)
							/ (da * db);
				} 
			}
		}

		p.StepIt();

		DBL sigma = 0.1;
		int N = min(V.m, V.n);
		Vec<DBL> f = gaussFilt(N, sigma * N/2);

		desc.SetWindowTextA("Smoothing velocity profile.");
		convAll(V, f, p);

		//ignoring small velocities
		for(int i=0; i<m; i++)
		{
			for(int j=0; j<n; j++)
			{
				//if(fabs(V(i,j).x) < ignoreVal) V(i,j).x = 0;
				if(fabs(V(i,j).y) < ignoreVal) V(i,j).y = 0;
			}
		}

		lastX.copy(*X);
	}

	void drawSmoothed(Mat<DBL>& Vx, Mat<DBL>& Vy)
	{
		VecF2 rx = VecF2(0, size.x * MU::len);
		VecF2 ry = VecF2(0, size.y * MU::len);

		Graph2 g = Graph2();
		g.create(FS::result, "vSmooth", -1, 1, 2);
		g.setLabel("y [m]", "x [m]", "\\v_x [m/s]"); g.setData(Vx); g.setRange(rx, ry); g.plot(VecL3(1,2,0));
		g.setLabel("y [m]", "x [m]", "\\v_y [m/s]"); g.setData(Vy); g.setRange(rx, ry); g.plot(VecL3(1,2,1));
		g.draw();
	}

	void txt(int step, string name)
	{
		MyWriteFile f = MyWriteFile(FS::getFN(FS::result, name, "txt", step), false);
		f.setSep("");
		for(int i=0; i<X->m; i++)
		{
			int k = X->n;
			for(int j=0; j<k; j++)
			{
				f.write((*X)(i,j));
				if(j < k-1) f.write(";");
			}
			f.newLine();
		}
		f.close();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PDE1);
		for(int i=0; i<m*n; i++)
			f.write(X->v[i]);
	}

	void load(MyReadFile &f)
	{		
		for(int i=0; i<m*n; i++) 
			f.read(X->v[i]);
		lastX.copy(*X);
	}

	DBL getTotalCon()
	{
		DBL con = 0;
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				con += (*X)(i,j);
		return con * dx * dx;
	}

private:
	void copyX(bool closedWound)
	{
		Mat<DBL> *p = X;
		X = Xnew;
		Xnew = p;
		Xnew->init(0);
		if(!closedWound) openedWound();
	}

	//implicit in x direction
	void setAx(TDM &Ax, int j, Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		//left D bc
		//Ax.a[0] = 0;
		//Ax.b[0] = 1;
		//Ax.c[0] = 0;

		//right D bc
		//Ax.a[m-1] = 0;
		//Ax.b[m-1] = 1;
		//Ax.c[m-1] = 0;

		//left N bc
		int i=0;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) - t2x * V(i,j).x + t2*gamma3;
			Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ax.c[i] = - at2x2;
		}
		else
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ax.c[i] = - at2x2;
		}
		//Ax.a[i] = 0;
		//Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ax.c[i] = - at2x2;
		
		//right N bc
		i = m-1;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ax.c[i] = 0;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
			Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2x * V(i,j).x + t2*gamma3;
			Ax.c[i] = 0;
		}
		else
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ax.c[i] = 0;
		}
		//Ax.a[i] = - at2x2;
		//Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ax.c[i] = 0;

		for(int i=1; i<m-1; i++)
		{
			if(V(i,j).x < 0)
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) - t2x * V(i,j).x + t2*gamma3;
				Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
			}
			else if(V(i,j).x > 0)
			{
				Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
				Ax.b[i] = atx2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2x * V(i,j).x + t2*gamma3;
				Ax.c[i] = - at2x2;
			}
			else
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
				Ax.c[i] = - at2x2;
			}
		}
	}

	//implicit in y direction
	void setAy(TDM &Ay, int i, Mat<DBL> &Tm, Mat<DBL> &Tf)
	{	
		//bottom N bc
		int j=0;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) - t2y * V(i,j).y + t2*gamma3;
			Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ay.c[j] = - at2y2;
		}
		else
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ay.c[j] = - at2y2;
		}
		//Ay.a[j] = 0;
		//Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ay.c[j] = - at2y2;

		//top N bc
		j = n-1;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ay.c[j] = 0;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
			Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2y * V(i,j).y + t2*gamma3;
			Ay.c[j] = 0;
		}
		else
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
			Ay.c[j] = 0;
		}
		//Ay.a[j] = - at2y2;
		//Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ay.c[j] = 0;
			
		for(int j=1; j<n-1; j++)
		{
			if(V(i,j).y < 0)
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) - t2y * V(i,j).y + t2*gamma3;
				Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
			}
			else if(V(i,j).y > 0)
			{
				Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
				Ay.b[j] = aty2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2y * V(i,j).y + t2*gamma3;
				Ay.c[j] = - at2y2;
			}
			else
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. + t2 * f((*X)(i,j), Tm(i,j), Tf(i,j)) + t2*gamma3;
				Ay.c[j] = - at2y2;
			}
		}
	}

	//explicit in y direction
	void setBx(Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		//left D bc
		//B(0,j) = 0;
		//int i0 = 1;

		//right D bc
		//B(m-1,j) = 0;
		//int iLim = m-1;		

		//left N bc
		int i0 = 0;
		//right N bc
		int iLim = m;

		//bottom N bc
		int j = 0;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j).y) * (*X)(i,j) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2y2 * (*X)(i,j+1);
			else
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2y2 * (*X)(i,j+1);
		}

		//top N bc
		j = n-1;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ (at2y2 - t2y * V(i,j).y) * (*X)(i,j); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2y2 * (*X)(i,j);
			else
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2y2 * (*X)(i,j);
		}
		
		#pragma omp parallel for
		for(int j=1; j<n-1; j++)
		{
			for(int i=i0; i<iLim; i++)
			{
				if(V(i,j).y < 0)
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
						+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
				else if(V(i,j).y > 0)
					B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
						+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
						+ at2y2 * (*X)(i,j+1);
				else
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
						+ at2y2 * (*X)(i,j+1);
			}
		}
		
	}

	//explicit in x direction
	void setBy(Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		//bottom D bc
		//B(i,0) = 0;
		//int j0 = 1;
		//top D bc
		//B(i,n-1) = 0;		
		//int jLim = n-1;
			
		//bottom N bc
		int j0 = 0;
		//top N bc
		int jLim = n;
		
		//left N bc
		int i = 0;
		for(int j=j0; j<jLim; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i,j).x) * (*X)(i,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2x2 * (*X)(i+1,j);
			else
			{
				B(i,j) = at2x2 * (*X)(i,j) 
					+ (1. - atx2) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2x2 * (*X)(i+1,j);
			}
		}

		//right N bc
		i = m-1;
		for(int j=j0; j<jLim; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ (at2x2 - t2x * V(i,j).x) * (*X)(i,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2x2 * (*X)(i,j);
			else
			{
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
					+ at2x2 * (*X)(i,j);
			}
		}

		#pragma omp parallel for
		for(int i=1; i<m-1; i++)
		{
			for(int j=j0; j<jLim; j++)
			{
				if(V(i,j).x < 0)
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
						+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
				else if(V(i,j).x > 0)
					B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
						+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
						+ at2x2 * (*X)(i+1,j);
				else
				{
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2) * (*X)(i,j) + t2 * C0 * f((*X)(i,j), Tm(i,j), Tf(i,j))
						+ at2x2 * (*X)(i+1,j);
				}
			}
				
		}
	}

	//implicit in x, explicit in y direction
	void thomasX(Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		//0 <= j < n-1
		#pragma omp parallel for
		for(int j=0; j<n; j++)
		{
			vector<DBL> c = vector<DBL>(m-1, 0);
			vector<DBL> d = vector<DBL>(m, 0);
			
			TDM Ax = TDM(m);
			setAx(Ax, j, Tm, Tf);

			c[0] = Ax.c[0] / Ax.b[0];
			for(int i=1; i<m-1; i++)
				c[i] = Ax.c[i] / (Ax.b[i] - c[i-1] * Ax.a[i]);

			d[0] = B(0,j) / Ax.b[0];
			for(int i=1; i<m; i++)
			{
				d[i] = (B(i,j) - d[i-1] * Ax.a[i]) / (Ax.b[i] - c[i-1] * Ax.a[i]);
			}

			(*Xnew)(m-1, j) = d[m-1];
			for(int i=m-2; i>=0; i--)
			{
				(*Xnew)(i,j) = d[i] - c[i]*(*Xnew)(i+1,j);
			}

		}	
	}

	//implicit in y, explicit in x direction
	void thomasY(Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		//0 < i < m-1
		#pragma omp parallel for
		for(int i=0; i<m; i++)	//Neumann
		//for(int i=1; i<m-1; i++)	//Dirichlet
		{
			vector<DBL> c = vector<DBL>(n-1, 0);
			vector<DBL> d = vector<DBL>(n, 0);

			TDM Ay = TDM(n);
			setAy(Ay, i, Tm, Tf);
			
			c[0] = Ay.c[0] / Ay.b[0];
			for(int j=1; j<n-1; j++)
				c[j] = Ay.c[j] / (Ay.b[j] - c[j-1] * Ay.a[j]);

			d[0] = B(i,0) / Ay.b[0];
			for(int j=1; j<n; j++)
			{
				d[j] = (B(i,j) - d[j-1] * Ay.a[j]) / (Ay.b[j] - c[j-1] * Ay.a[j]);
			}

			(*Xnew)(i, n-1) = d[n-1];
			for(int j=n-2; j>=0; j--)
			{
				(*Xnew)(i,j) = d[j] - c[j]*(*Xnew)(i,j+1);	
			}
		}	
	}

	//implicit in x, explicit in y direction
	void halfStepX(bool closedWound,  Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		setBx(Tm, Tf);
		thomasX(Tm, Tf);
		copyX(closedWound);
	}

	//implicit in y, explicit in x direction
	void halfStepY(bool closedWound, Mat<DBL> &Tm, Mat<DBL> &Tf)
	{
		setBy(Tm, Tf);
		thomasY(Tm, Tf);
		copyX(closedWound);
	}

	inline DBL f(DBL T, DBL Tm, DBL Tf)
	{
		return (beta2 * pow(T, 2) / (T0 + T)) + (beta1 * Tf) + (beta3 * Tm * T);
	}

	void openedWound()
	{
		for(int i=0; i<m; i++)	
		{
			DBL x = i*dx;		

			if(x >= wound.x && x <= wound.y)
				(*X)(i,0) = 1;
		}
	}
};



//Tissue factor
class Pde4
{
public:
	DBL dx;
	DBL dy;
	DBL dt;
	Mat<DBL> *X;
	Mat<DBL> *Xnew;
	int m;
	int n;
	DBL alpha1;
	DBL k1;
	DBL gamma1;
	DBL C0;
	DBL T0;
	Mat<DBL> B;
	DBL t2, t2x, at2x2, atx2, t2y, at2y2, aty2;
	VecF2 size;
	vector<DBL> posX;
	Mat<VecF2> V;
	Mat<DBL> lastX;
	VecF2 wound;

	Pde4()
	{
		X		= NULL;
		Xnew	= NULL;
		V		= Mat<VecF2>();
		lastX	= Mat<DBL>();
		B		= Mat<DBL>();
	}

	Pde4(Setup &s)
	{
		if(!s.fib.enabled) return;
		size = s.spc.size;
		alpha1 = s.fib.alpha1;
		k1 = s.fib.k1;
		gamma1 = s.fib.gamma1;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dx = s.fib.dx;
		dy = s.fib.dx;
		dt = s.fib.dt;
		wound = VecF2(s.fib.woundBeg, s.fib.woundEnd);
		
		m = size.x / dx + 1;
		n = size.y / dx + 1;
		X = new Mat<DBL>(m, n);
		B = Mat<DBL>(m, n);
		Xnew = new Mat<DBL>(m, n);
		Mat<DBL> B(m, n);
		V = Mat<VecF2>(m,n);
		lastX = Mat<DBL>(m,n);

		V.init(VecF20);
		//initial conditions
		X->init(0);
		Xnew->init(0);
		openedWound();

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha1 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha1 * t2y / dy;
		aty2 = at2y2 * 2;

		lastX.copy(*X);
	}

	void paramHotChange(Setup &s)
	{
		alpha1 = s.fib.alpha1;
		k1 = s.fib.k1;
		gamma1 = s.fib.gamma1;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dt = s.fib.dt;
		wound = VecF2(s.fib.woundBeg, s.fib.woundEnd);

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha1 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha1 * t2y / dy;
		aty2 = at2y2 * 2;
	}

	void nextStep(bool closedWound)
	{
		halfStepX(closedWound);
		halfStepY(closedWound);
	}

	void free()
	{
		if(X != NULL) 
		{
			delete X;
			X = NULL;
		}
		if(Xnew != NULL) 
		{
			delete Xnew;
			Xnew = NULL;
		}
	}

	void saveWave()
	{
		int y = n / 2;
		for(int i=0; i<m; i++)
		{
			if(((*X)(i,y)-0.5)*((*X)(i+1,y)-0.5) < 0)
			{
				DBL aaa = (*X)(i,y)-(*X)(i+1,y);
				posX.push_back(((*X)(i,y)/aaa-0.5/aaa+i)*dx);
				break;
			}
		}
	}


	//returns the concentration on given coordinates p
	inline DBL getU(VecF2 &p)
	{	
		int i = p.x / dx;
		int j = p.y / dy;
		return bilin(p.x, p.y, dx*i, dx*(i+1), dy*j, dy*(j+1), lastX(i, j), lastX(i+1, j), lastX(i, j+1), lastX(i+1, j+1));
		//return (lastX(i, j) + lastX(i+1, j) + lastX(i, j+1) + lastX(i+1, j+1)) / 4;
		
	}

	void update(DBL coeff)
	{
		for(int j=0; j<n; j++)
		{
			DBL y = dy * j;
			for(int i=0; i<m; i++)
			{
				V(i,j) = VecF2(coeff * y*(size.y - y) * 4. / pow(size.y, 2), 0);
				//V(i,j) = VecF2(coeff,0);
			}
		}

		lastX.copy(*X);
	}

	void update(VecF2 v)
	{
		for(int j=0; j<n; j++)
			for(int i=0; i<m; i++)
				V(i,j) = v;
		lastX.copy(*X);
	}

	void update(Mat<Stat> &M, DBL ignoreVal, CProgressCtrl &p, CStatic &desc)
	{
		DBL da = size.x / M.m;
		DBL db = size.y / M.n;
		DBL dah = da / 2;
		DBL dbh = db / 2;
		DBL daL = size.x - da;
		DBL dbL = size.y - db;

		int i1=0, i2=0, j1=0, j2=0;
		DBL dxL=0, dx1=0, dx2=0, dyL=0, dy1=0, dy2=0;

		desc.SetWindowTextA("Adjusting velocity profile for PDE.");
		for(int i=0; i<m; i++)
		{
			bool x0 = false;
			DBL x = dx * i - dah;

			if(x < 0)
			{
				x0 = true;
				i1 = 0;
			}
			else if(x >= daL)
			{
				x0 = true;
				i1 = M.m-1;
			}
			else
			{
				i1 = x / da;
				i2 = i1 + 1;
				dxL = da;
				dx1 = dx * i - (i1 * da + dah);
				dx2 = i2 * da + dah - dx * i;
			} 

			for(int j=0; j<n; j++)
			{
				DBL y = dy * j - dbh;

				if(y < 0)
				{
					dy1 = dy * j;
					if(x0) 
						V(i,j) = M(i1,0).v * dy1 / dbh;
					else 
						V(i,j) = ((M(i1,0).v * dx2 + M(i2,0).v * dx1)) * dy1 / (da * dbh);
				}
				else if(y >= dbL)
				{
					dy2 = size.y - dy * j;
					if(x0)
						V(i,j) = M(i1, M.n-1).v * dy2 / dbh;
					else
						V(i,j) = (M(i1, M.n-1).v * dx2 + M(i2, M.n-1).v * dx1) * dy2 / (da * dbh);
				}
				else
				{
					j1 = y / db;
					j2 = j1 + 1;
					dyL = db;
					dy1 = dy * j - (j1 * db + dbh);
					dy2 = j2 * db + dbh - dy * j;

					if(x0)
						V(i,j) = (M(i1, j1).v * dy2 + M(i1, j1).v * dy1) / db;
					else
						V(i,j) = (M(i1, j1).v * dx2 * dy2 
							+ M(i2, j1).v * dx1 * dy2
							+ M(i1, j2).v * dx2 * dy1
							+ M(i2, j2).v * dx1 * dy1)
							/ (da * db);
				} 
			}
		}

		p.StepIt();

		DBL sigma = 0.1;
		int N = min(V.m, V.n);
		Vec<DBL> f = gaussFilt(N, sigma * N/2);

		desc.SetWindowTextA("Smoothing velocity profile.");
		convAll(V, f, p);

		//ignoring small velocities
		for(int i=0; i<m; i++)
		{
			for(int j=0; j<n; j++)
			{
				//if(fabs(V(i,j).x) < ignoreVal) V(i,j).x = 0;
				if(fabs(V(i,j).y) < ignoreVal) V(i,j).y = 0;
			}
		}

		lastX.copy(*X);
	}

	void drawSmoothed(Mat<DBL>& Vx, Mat<DBL>& Vy)
	{
		VecF2 rx = VecF2(0, size.x * MU::len);
		VecF2 ry = VecF2(0, size.y * MU::len);

		Graph2 g = Graph2();
		g.create(FS::result, "vSmooth", -1, 1, 2);
		g.setLabel("y [m]", "x [m]", "\\v_x [m/s]"); g.setData(Vx); g.setRange(rx, ry); g.plot(VecL3(1,2,0));
		g.setLabel("y [m]", "x [m]", "\\v_y [m/s]"); g.setData(Vy); g.setRange(rx, ry); g.plot(VecL3(1,2,1));
		g.draw();
	}

	void txt(int step, string name)
	{
		MyWriteFile f = MyWriteFile(FS::getFN(FS::result, name, "txt", step), false);
		f.setSep("");
		for(int i=0; i<X->m; i++)
		{
			int k = X->n;
			for(int j=0; j<k; j++)
			{
				f.write((*X)(i,j));
				if(j < k-1) f.write(";");
			}
			f.newLine();
		}
		f.close();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PDE4);
		for(int i=0; i<m*n; i++)
			f.write(X->v[i]);
	}

	void load(MyReadFile &f)
	{		
		for(int i=0; i<m*n; i++) 
			f.read(X->v[i]);
		lastX.copy(*X);
	}

	DBL getTotalCon()
	{
		DBL con = 0;
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				con += (*X)(i,j);
		return con * dx * dx;
	}

private:
	void copyX(bool closedWound)
	{
		Mat<DBL> *p = X;
		X = Xnew;
		Xnew = p;
		Xnew->init(0);
		if(!closedWound) openedWound();
	}

	//implicit in x direction
	void setAx(TDM &Ax, int j)
	{
		//left D bc
		//Ax.a[0] = 0;
		//Ax.b[0] = 1;
		//Ax.c[0] = 0;

		//right D bc
		//Ax.a[m-1] = 0;
		//Ax.b[m-1] = 1;
		//Ax.c[m-1] = 0;

		//left N bc
		int i=0;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. - t2x * V(i,j).x + t2*gamma1;
			Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2*gamma1;
			Ax.c[i] = - at2x2;
		}
		else
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2*gamma1;
			Ax.c[i] = - at2x2;
		}
		//Ax.a[i] = 0;
		//Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ax.c[i] = - at2x2;
		
		//right N bc
		i = m-1;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2*gamma1;
			Ax.c[i] = 0;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
			Ax.b[i] = at2x2 + 1. + t2x * V(i,j).x + t2*gamma1;
			Ax.c[i] = 0;
		}
		else
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2*gamma1;
			Ax.c[i] = 0;
		}
		//Ax.a[i] = - at2x2;
		//Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ax.c[i] = 0;

		for(int i=1; i<m-1; i++)
		{
			if(V(i,j).x < 0)
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. - t2x * V(i,j).x + t2*gamma1;
				Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
			}
			else if(V(i,j).x > 0)
			{
				Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
				Ax.b[i] = atx2 + 1. + t2x * V(i,j).x + t2*gamma1;
				Ax.c[i] = - at2x2;
			}
			else
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. + t2*gamma1;
				Ax.c[i] = - at2x2;
			}
		}
	}

	//implicit in y direction
	void setAy(TDM &Ay, int i)
	{	
		//bottom N bc
		int j=0;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. - t2y * V(i,j).y + t2*gamma1;
			Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2*gamma1;
			Ay.c[j] = - at2y2;
		}
		else
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2*gamma1;
			Ay.c[j] = - at2y2;
		}
		//Ay.a[j] = 0;
		//Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ay.c[j] = - at2y2;

		//top N bc
		j = n-1;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2*gamma1;
			Ay.c[j] = 0;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
			Ay.b[j] = at2y2 + 1. + t2y * V(i,j).y + t2*gamma1;
			Ay.c[j] = 0;
		}
		else
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2*gamma1;
			Ay.c[j] = 0;
		}
		//Ay.a[j] = - at2y2;
		//Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ay.c[j] = 0;
			
		for(int j=1; j<n-1; j++)
		{
			if(V(i,j).y < 0)
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. - t2y * V(i,j).y + t2*gamma1;
				Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
			}
			else if(V(i,j).y > 0)
			{
				Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
				Ay.b[j] = aty2 + 1. + t2y * V(i,j).y + t2*gamma1;
				Ay.c[j] = - at2y2;
			}
			else
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. + t2*gamma1;
				Ay.c[j] = - at2y2;
			}
		}
	}

	//explicit in y direction
	void setBx()
	{
		//left D bc
		//B(0,j) = 0;
		//int i0 = 1;

		//right D bc
		//B(m-1,j) = 0;
		//int iLim = m-1;		

		//left N bc
		int i0 = 0;
		//right N bc
		int iLim = m;

		//bottom N bc
		int j = 0;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j) 
					+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j).y) * (*X)(i,j) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) 
					+ at2y2 * (*X)(i,j+1);
			else
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2) * (*X)(i,j) +
					+ at2y2 * (*X)(i,j+1);
		}

		//top N bc
		j = n-1;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)  
					+ (at2y2 - t2y * V(i,j).y) * (*X)(i,j); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) 
					+ at2y2 * (*X)(i,j);
			else
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2) * (*X)(i,j) + 
					+ at2y2 * (*X)(i,j);
		}
		
		#pragma omp parallel for
		for(int j=1; j<n-1; j++)
		{
			for(int i=i0; i<iLim; i++)
			{
				if(V(i,j).y < 0)
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)  
						+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
				else if(V(i,j).y > 0)
					B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
						+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j)  
						+ at2y2 * (*X)(i,j+1);
				else
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2) * (*X)(i,j)  
						+ at2y2 * (*X)(i,j+1);
			}
		}
		
	}

	//explicit in x direction
	void setBy()
	{
		//bottom D bc
		//B(i,0) = 0;
		//int j0 = 1;
		//top D bc
		//B(i,n-1) = 0;		
		//int jLim = n-1;
			
		//bottom N bc
		int j0 = 0;
		//top N bc
		int jLim = n;
		
		//left N bc
		int i = 0;
		for(int j=j0; j<jLim; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) 
					+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i,j).x) * (*X)(i,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j)
					+ at2x2 * (*X)(i+1,j);
			else
			{
				B(i,j) = at2x2 * (*X)(i,j) 
					+ (1. - atx2) * (*X)(i,j)
					+ at2x2 * (*X)(i+1,j);
			}
		}

		//right N bc
		i = m-1;
		for(int j=j0; j<jLim; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j)
					+ (at2x2 - t2x * V(i,j).x) * (*X)(i,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) 
					+ at2x2 * (*X)(i,j);
			else
			{
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2) * (*X)(i,j)
					+ at2x2 * (*X)(i,j);
			}
		}

		#pragma omp parallel for
		for(int i=1; i<m-1; i++)
		{
			for(int j=j0; j<jLim; j++)
			{
				if(V(i,j).x < 0)
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) 
						+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
				else if(V(i,j).x > 0)
					B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
						+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) 
						+ at2x2 * (*X)(i+1,j);
				else
				{
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2) * (*X)(i,j)
						+ at2x2 * (*X)(i+1,j);
				}
			}
				
		}
	}

	//implicit in x, explicit in y direction
	void thomasX()
	{
		//0 <= j < n-1
		#pragma omp parallel for
		for(int j=0; j<n; j++)
		{
			vector<DBL> c = vector<DBL>(m-1, 0);
			vector<DBL> d = vector<DBL>(m, 0);
			
			TDM Ax = TDM(m);
			setAx(Ax, j);

			c[0] = Ax.c[0] / Ax.b[0];
			for(int i=1; i<m-1; i++)
				c[i] = Ax.c[i] / (Ax.b[i] - c[i-1] * Ax.a[i]);

			d[0] = B(0,j) / Ax.b[0];
			for(int i=1; i<m; i++)
			{
				d[i] = (B(i,j) - d[i-1] * Ax.a[i]) / (Ax.b[i] - c[i-1] * Ax.a[i]);
			}

			(*Xnew)(m-1, j) = d[m-1];
			for(int i=m-2; i>=0; i--)
			{
				(*Xnew)(i,j) = d[i] - c[i]*(*Xnew)(i+1,j);
			}

		}	
	}

	//implicit in y, explicit in x direction
	void thomasY()
	{
		//0 < i < m-1
		#pragma omp parallel for
		for(int i=0; i<m; i++)	//Neumann
		//for(int i=1; i<m-1; i++)	//Dirichlet
		{
			vector<DBL> c = vector<DBL>(n-1, 0);
			vector<DBL> d = vector<DBL>(n, 0);

			TDM Ay = TDM(n);
			setAy(Ay, i);
			
			c[0] = Ay.c[0] / Ay.b[0];
			for(int j=1; j<n-1; j++)
				c[j] = Ay.c[j] / (Ay.b[j] - c[j-1] * Ay.a[j]);

			d[0] = B(i,0) / Ay.b[0];
			for(int j=1; j<n; j++)
			{
				d[j] = (B(i,j) - d[j-1] * Ay.a[j]) / (Ay.b[j] - c[j-1] * Ay.a[j]);
			}

			(*Xnew)(i, n-1) = d[n-1];
			for(int j=n-2; j>=0; j--)
			{
				(*Xnew)(i,j) = d[j] - c[j]*(*Xnew)(i,j+1);	
			}
		}	
	}

	//implicit in x, explicit in y direction
	void halfStepX(bool closedWound)
	{
		setBx();
		thomasX();
		copyX(closedWound);
	}

	//implicit in y, explicit in x direction
	void halfStepY(bool closedWound)
	{
		setBy();
		thomasY();
		copyX(closedWound);
	}



	void openedWound()
	{
		for(int i=0; i<m; i++)	
		{
			DBL x = i*dx;		

			if(x >= wound.x && x <= wound.y)
				(*X)(i,0) = 1;
		}
	}
};




class Pde5
{
public:
	DBL dx;
	DBL dy;
	DBL dt;
	Mat<DBL> *X;
	Mat<DBL> *Xnew;
	int m;
	int n;
	DBL alpha2;
	DBL k1;
	DBL gamma2;
	DBL C0;
	DBL T0;
	Mat<DBL> B;
	DBL t2, t2x, at2x2, atx2, t2y, at2y2, aty2;
	VecF2 size;
	vector<DBL> posX;
	Mat<VecF2> V;
	Mat<DBL> lastX;
	VecF2 wound;

	Pde5()
	{
		X		= NULL;
		Xnew	= NULL;
		V		= Mat<VecF2>();
		lastX	= Mat<DBL>();
		B		= Mat<DBL>();
	}

	Pde5(Setup &s)
	{
		if(!s.fib.enabled) return;
		size = s.spc.size;
		alpha2 = s.fib.alpha2;
		k1 = s.fib.k1;
		gamma2 = s.fib.gamma2;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dx = s.fib.dx;
		dy = s.fib.dx;
		dt = s.fib.dt;
		wound = VecF2(s.fib.woundBeg, s.fib.woundEnd);
		
		m = size.x / dx + 1;
		n = size.y / dx + 1;
		X = new Mat<DBL>(m, n);
		B = Mat<DBL>(m, n);
		Xnew = new Mat<DBL>(m, n);
		Mat<DBL> B(m, n);
		V = Mat<VecF2>(m,n);
		lastX = Mat<DBL>(m,n);

		V.init(VecF20);
		//initial conditions
		X->init(0);
		Xnew->init(0);
		openedWound();

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha2 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha2 * t2y / dy;
		aty2 = at2y2 * 2;

		lastX.copy(*X);
	}

	void paramHotChange(Setup &s)
	{
		alpha2 = s.fib.alpha2;
		k1 = s.fib.k1;
		gamma2 = s.fib.gamma2;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dt = s.fib.dt;
		wound = VecF2(s.fib.woundBeg, s.fib.woundEnd);

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha2 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha2 * t2y / dy;
		aty2 = at2y2 * 2;
	}

	void nextStep(bool closedWound)
	{
		halfStepX(closedWound);
		halfStepY(closedWound);
	}

	void free()
	{
		if(X != NULL) 
		{
			delete X;
			X = NULL;
		}
		if(Xnew != NULL) 
		{
			delete Xnew;
			Xnew = NULL;
		}
	}

	void saveWave()
	{
		int y = n / 2;
		for(int i=0; i<m; i++)
		{
			if(((*X)(i,y)-0.5)*((*X)(i+1,y)-0.5) < 0)
			{
				DBL aaa = (*X)(i,y)-(*X)(i+1,y);
				posX.push_back(((*X)(i,y)/aaa-0.5/aaa+i)*dx);
				break;
			}
		}
	}


	//returns the concentration on given coordinates p
	inline DBL getU(VecF2 &p)
	{	
		int i = p.x / dx;
		int j = p.y / dy;
		return bilin(p.x, p.y, dx*i, dx*(i+1), dy*j, dy*(j+1), lastX(i, j), lastX(i+1, j), lastX(i, j+1), lastX(i+1, j+1));
		//return (lastX(i, j) + lastX(i+1, j) + lastX(i, j+1) + lastX(i+1, j+1)) / 4;
		
	}

	void update(DBL coeff)
	{
		for(int j=0; j<n; j++)
		{
			DBL y = dy * j;
			for(int i=0; i<m; i++)
			{
				V(i,j) = VecF2(coeff * y*(size.y - y) * 4. / pow(size.y, 2), 0);
				//V(i,j) = VecF2(coeff,0);
			}
		}

		lastX.copy(*X);
	}

	void update(VecF2 v)
	{
		for(int j=0; j<n; j++)
			for(int i=0; i<m; i++)
				V(i,j) = v;
		lastX.copy(*X);
	}

	void update(Mat<Stat> &M, DBL ignoreVal, CProgressCtrl &p, CStatic &desc)
	{
		DBL da = size.x / M.m;
		DBL db = size.y / M.n;
		DBL dah = da / 2;
		DBL dbh = db / 2;
		DBL daL = size.x - da;
		DBL dbL = size.y - db;

		int i1=0, i2=0, j1=0, j2=0;
		DBL dxL=0, dx1=0, dx2=0, dyL=0, dy1=0, dy2=0;

		desc.SetWindowTextA("Adjusting velocity profile for PDE.");
		for(int i=0; i<m; i++)
		{
			bool x0 = false;
			DBL x = dx * i - dah;

			if(x < 0)
			{
				x0 = true;
				i1 = 0;
			}
			else if(x >= daL)
			{
				x0 = true;
				i1 = M.m-1;
			}
			else
			{
				i1 = x / da;
				i2 = i1 + 1;
				dxL = da;
				dx1 = dx * i - (i1 * da + dah);
				dx2 = i2 * da + dah - dx * i;
			} 

			for(int j=0; j<n; j++)
			{
				DBL y = dy * j - dbh;

				if(y < 0)
				{
					dy1 = dy * j;
					if(x0) 
						V(i,j) = M(i1,0).v * dy1 / dbh;
					else 
						V(i,j) = ((M(i1,0).v * dx2 + M(i2,0).v * dx1)) * dy1 / (da * dbh);
				}
				else if(y >= dbL)
				{
					dy2 = size.y - dy * j;
					if(x0)
						V(i,j) = M(i1, M.n-1).v * dy2 / dbh;
					else
						V(i,j) = (M(i1, M.n-1).v * dx2 + M(i2, M.n-1).v * dx1) * dy2 / (da * dbh);
				}
				else
				{
					j1 = y / db;
					j2 = j1 + 1;
					dyL = db;
					dy1 = dy * j - (j1 * db + dbh);
					dy2 = j2 * db + dbh - dy * j;

					if(x0)
						V(i,j) = (M(i1, j1).v * dy2 + M(i1, j1).v * dy1) / db;
					else
						V(i,j) = (M(i1, j1).v * dx2 * dy2 
							+ M(i2, j1).v * dx1 * dy2
							+ M(i1, j2).v * dx2 * dy1
							+ M(i2, j2).v * dx1 * dy1)
							/ (da * db);
				} 
			}
		}

		p.StepIt();

		DBL sigma = 0.1;
		int N = min(V.m, V.n);
		Vec<DBL> f = gaussFilt(N, sigma * N/2);

		desc.SetWindowTextA("Smoothing velocity profile.");
		convAll(V, f, p);

		//ignoring small velocities
		for(int i=0; i<m; i++)
		{
			for(int j=0; j<n; j++)
			{
				//if(fabs(V(i,j).x) < ignoreVal) V(i,j).x = 0;
				if(fabs(V(i,j).y) < ignoreVal) V(i,j).y = 0;
			}
		}

		lastX.copy(*X);
	}

	void drawSmoothed(Mat<DBL>& Vx, Mat<DBL>& Vy)
	{
		VecF2 rx = VecF2(0, size.x * MU::len);
		VecF2 ry = VecF2(0, size.y * MU::len);

		Graph2 g = Graph2();
		g.create(FS::result, "vSmooth", -1, 1, 2);
		g.setLabel("y [m]", "x [m]", "\\v_x [m/s]"); g.setData(Vx); g.setRange(rx, ry); g.plot(VecL3(1,2,0));
		g.setLabel("y [m]", "x [m]", "\\v_y [m/s]"); g.setData(Vy); g.setRange(rx, ry); g.plot(VecL3(1,2,1));
		g.draw();
	}

	void txt(int step, string name)
	{
		MyWriteFile f = MyWriteFile(FS::getFN(FS::result, name, "txt", step), false);
		f.setSep("");
		for(int i=0; i<X->m; i++)
		{
			int k = X->n;
			for(int j=0; j<k; j++)
			{
				f.write((*X)(i,j));
				if(j < k-1) f.write(";");
			}
			f.newLine();
		}
		f.close();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PDE5);
		for(int i=0; i<m*n; i++)
			f.write(X->v[i]);
	}

	void load(MyReadFile &f)
	{		
		for(int i=0; i<m*n; i++) 
			f.read(X->v[i]);
		lastX.copy(*X);
	}

	DBL getTotalCon()
	{
		DBL con = 0;
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				con += (*X)(i,j);
		return con * dx * dx;
	}

private:
	void copyX(bool closedWound)
	{
		Mat<DBL> *p = X;
		X = Xnew;
		Xnew = p;
		Xnew->init(0);
		if(!closedWound) openedWound();
	}

	//implicit in x direction
	void setAx(TDM &Ax, int j)
	{
		//left D bc
		//Ax.a[0] = 0;
		//Ax.b[0] = 1;
		//Ax.c[0] = 0;

		//right D bc
		//Ax.a[m-1] = 0;
		//Ax.b[m-1] = 1;
		//Ax.c[m-1] = 0;

		//left N bc
		int i=0;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. - t2x * V(i,j).x + t2*gamma2;
			Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2*gamma2;
			Ax.c[i] = - at2x2;
		}
		else
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2*gamma2;
			Ax.c[i] = - at2x2;
		}
		//Ax.a[i] = 0;
		//Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ax.c[i] = - at2x2;
		
		//right N bc
		i = m-1;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2*gamma2;
			Ax.c[i] = 0;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
			Ax.b[i] = at2x2 + 1. + t2x * V(i,j).x + t2*gamma2;
			Ax.c[i] = 0;
		}
		else
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2*gamma2;
			Ax.c[i] = 0;
		}
		//Ax.a[i] = - at2x2;
		//Ax.b[i] = at2x2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ax.c[i] = 0;

		for(int i=1; i<m-1; i++)
		{
			if(V(i,j).x < 0)
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. - t2x * V(i,j).x + t2*gamma2;
				Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
			}
			else if(V(i,j).x > 0)
			{
				Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
				Ax.b[i] = atx2 + 1. + t2x * V(i,j).x + t2*gamma2;
				Ax.c[i] = - at2x2;
			}
			else
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. + t2*gamma2;
				Ax.c[i] = - at2x2;
			}
		}
	}

	//implicit in y direction
	void setAy(TDM &Ay, int i)
	{	
		//bottom N bc
		int j=0;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. - t2y * V(i,j).y + t2*gamma2;
			Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2*gamma2;
			Ay.c[j] = - at2y2;
		}
		else
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2*gamma2;
			Ay.c[j] = - at2y2;
		}
		//Ay.a[j] = 0;
		//Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ay.c[j] = - at2y2;

		//top N bc
		j = n-1;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2*gamma2;
			Ay.c[j] = 0;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
			Ay.b[j] = at2y2 + 1. + t2y * V(i,j).y + t2*gamma2;
			Ay.c[j] = 0;
		}
		else
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2*gamma2;
			Ay.c[j] = 0;
		}
		//Ay.a[j] = - at2y2;
		//Ay.b[j] = at2y2 + 1. + t2 * f((*X)(i,j)) + t2*gamma;
		//Ay.c[j] = 0;
			
		for(int j=1; j<n-1; j++)
		{
			if(V(i,j).y < 0)
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. - t2y * V(i,j).y + t2*gamma2;
				Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
			}
			else if(V(i,j).y > 0)
			{
				Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
				Ay.b[j] = aty2 + 1. + t2y * V(i,j).y + t2*gamma2;
				Ay.c[j] = - at2y2;
			}
			else
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. + t2*gamma2;
				Ay.c[j] = - at2y2;
			}
		}
	}

	//explicit in y direction
	void setBx()
	{
		//left D bc
		//B(0,j) = 0;
		//int i0 = 1;

		//right D bc
		//B(m-1,j) = 0;
		//int iLim = m-1;		

		//left N bc
		int i0 = 0;
		//right N bc
		int iLim = m;

		//bottom N bc
		int j = 0;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j) 
					+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j).y) * (*X)(i,j) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) 
					+ at2y2 * (*X)(i,j+1);
			else
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2) * (*X)(i,j) +
					+ at2y2 * (*X)(i,j+1);
		}

		//top N bc
		j = n-1;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)  
					+ (at2y2 - t2y * V(i,j).y) * (*X)(i,j); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j) 
					+ at2y2 * (*X)(i,j);
			else
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2) * (*X)(i,j) + 
					+ at2y2 * (*X)(i,j);
		}
		
		#pragma omp parallel for
		for(int j=1; j<n-1; j++)
		{
			for(int i=i0; i<iLim; i++)
			{
				if(V(i,j).y < 0)
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)  
						+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
				else if(V(i,j).y > 0)
					B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
						+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j)  
						+ at2y2 * (*X)(i,j+1);
				else
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2) * (*X)(i,j)  
						+ at2y2 * (*X)(i,j+1);
			}
		}
		
	}

	//explicit in x direction
	void setBy()
	{
		//bottom D bc
		//B(i,0) = 0;
		//int j0 = 1;
		//top D bc
		//B(i,n-1) = 0;		
		//int jLim = n-1;
			
		//bottom N bc
		int j0 = 0;
		//top N bc
		int jLim = n;
		
		//left N bc
		int i = 0;
		for(int j=j0; j<jLim; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) 
					+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i,j).x) * (*X)(i,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j)
					+ at2x2 * (*X)(i+1,j);
			else
			{
				B(i,j) = at2x2 * (*X)(i,j) 
					+ (1. - atx2) * (*X)(i,j)
					+ at2x2 * (*X)(i+1,j);
			}
		}

		//right N bc
		i = m-1;
		for(int j=j0; j<jLim; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j)
					+ (at2x2 - t2x * V(i,j).x) * (*X)(i,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) 
					+ at2x2 * (*X)(i,j);
			else
			{
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2) * (*X)(i,j)
					+ at2x2 * (*X)(i,j);
			}
		}

		#pragma omp parallel for
		for(int i=1; i<m-1; i++)
		{
			for(int j=j0; j<jLim; j++)
			{
				if(V(i,j).x < 0)
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j) 
						+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
				else if(V(i,j).x > 0)
					B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
						+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j) 
						+ at2x2 * (*X)(i+1,j);
				else
				{
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2) * (*X)(i,j)
						+ at2x2 * (*X)(i+1,j);
				}
			}
				
		}
	}

	//implicit in x, explicit in y direction
	void thomasX()
	{
		//0 <= j < n-1
		#pragma omp parallel for
		for(int j=0; j<n; j++)
		{
			vector<DBL> c = vector<DBL>(m-1, 0);
			vector<DBL> d = vector<DBL>(m, 0);
			
			TDM Ax = TDM(m);
			setAx(Ax, j);

			c[0] = Ax.c[0] / Ax.b[0];
			for(int i=1; i<m-1; i++)
				c[i] = Ax.c[i] / (Ax.b[i] - c[i-1] * Ax.a[i]);

			d[0] = B(0,j) / Ax.b[0];
			for(int i=1; i<m; i++)
			{
				d[i] = (B(i,j) - d[i-1] * Ax.a[i]) / (Ax.b[i] - c[i-1] * Ax.a[i]);
			}

			(*Xnew)(m-1, j) = d[m-1];
			for(int i=m-2; i>=0; i--)
			{
				(*Xnew)(i,j) = d[i] - c[i]*(*Xnew)(i+1,j);
			}

		}	
	}

	//implicit in y, explicit in x direction
	void thomasY()
	{
		//0 < i < m-1
		#pragma omp parallel for
		for(int i=0; i<m; i++)	//Neumann
		//for(int i=1; i<m-1; i++)	//Dirichlet
		{
			vector<DBL> c = vector<DBL>(n-1, 0);
			vector<DBL> d = vector<DBL>(n, 0);

			TDM Ay = TDM(n);
			setAy(Ay, i);
			
			c[0] = Ay.c[0] / Ay.b[0];
			for(int j=1; j<n-1; j++)
				c[j] = Ay.c[j] / (Ay.b[j] - c[j-1] * Ay.a[j]);

			d[0] = B(i,0) / Ay.b[0];
			for(int j=1; j<n; j++)
			{
				d[j] = (B(i,j) - d[j-1] * Ay.a[j]) / (Ay.b[j] - c[j-1] * Ay.a[j]);
			}

			(*Xnew)(i, n-1) = d[n-1];
			for(int j=n-2; j>=0; j--)
			{
				(*Xnew)(i,j) = d[j] - c[j]*(*Xnew)(i,j+1);	
			}
		}	
	}

	//implicit in x, explicit in y direction
	void halfStepX(bool closedWound)
	{
		setBx();
		thomasX();
		copyX(closedWound);
	}

	//implicit in y, explicit in x direction
	void halfStepY(bool closedWound)
	{
		setBy();
		thomasY();
		copyX(closedWound);
	}



	void openedWound()
	{
		for(int i=0; i<m; i++)	
		{
			DBL x = i*dx;		

			if(x >= wound.x && x <= wound.y)
				(*X)(i,0) = 1;
		}
	}
};





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//fibrinogen
// dFg/dt = alpha (d^2Fg/dx^2) - V.(Del Fg) - f(T) Fg
// f(T) = k3 T
class Pde2
{
public:
	DBL dx;
	DBL dy;
	DBL dt;
	Mat<DBL> *X;
	Mat<DBL> *Xnew;
	int m;
	int n;
	DBL alpha2;
	DBL beta3;
	Mat<DBL> B;
	DBL t2, t2x, at2x2, atx2, t2y, at2y2, aty2;
	VecF2 size;
	vector<DBL> posX;
	Mat<DBL> lastX;
	DBL initVal; 

	Pde2()
	{
		X		= NULL;
		Xnew	= NULL;
		lastX	= Mat<DBL>();
		B		= Mat<DBL>();
	}

	Pde2(Setup &s)
	{
		if(!s.fib.enabled) return;

		initVal = s.fib.initFg;

		size = s.spc.size;
		alpha2 = s.fib.alpha2;
		beta3 = s.fib.beta3;
		dx = s.fib.dx;
		dy = s.fib.dx;
		dt = s.fib.dt;
		
		m = size.x / dx + 1;
		n = size.y / dx + 1;
		X = new Mat<DBL>(m, n);
		B = Mat<DBL>(m, n);
		Xnew = new Mat<DBL>(m, n);
		Mat<DBL> B(m, n);
		lastX = Mat<DBL>(m,n);

		//initial conditions
		X->init(initVal);
		Xnew->init(initVal);

		for(int j=0; j<n; j++)
		{
			(*X)(m-1,j) = 0;
			(*Xnew)(m-1,j) = 0;
		}

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha2 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha2 * t2y / dy;
		aty2 = at2y2 * 2;

		lastX.copy(*X);
	}

	void paramHotChange(Setup &s)
	{
		initVal = s.fib.initFg;

		alpha2 = s.fib.alpha2;
		beta3 = s.fib.beta3;
		dt = s.fib.dt;

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha2 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha2 * t2y / dy;
		aty2 = at2y2 * 2;
	}

	void nextStep(Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{
		halfStepX(V, T, Fg, Fp);
		halfStepY(V, T, Fg, Fp);
	}

	void free()
	{
		if(X != NULL) 
		{
			delete X;
			X = NULL;
		}
		if(Xnew != NULL) 
		{
			delete Xnew;
			Xnew = NULL;
		}
	}

	void saveWave()
	{
		int y = n / 2;
		for(int i=0; i<m; i++)
		{
			if(((*X)(i,y)-0.5)*((*X)(i+1,y)-0.5) < 0)
			{
				DBL aaa = (*X)(i,y)-(*X)(i+1,y);
				posX.push_back(((*X)(i,y)/aaa-0.5/aaa+i)*dx);
				break;
			}
		}
	}


	//returns the concentration on given coordinates p
	inline DBL getU(VecF2 &p)
	{	
		int i = p.x / dx;
		int j = p.y / dy;
		return bilin(p.x, p.y, dx*i, dx*(i+1), dy*j, dy*(j+1), lastX(i, j), lastX(i+1, j), lastX(i, j+1), lastX(i+1, j+1));
		//return (lastX(i, j) + lastX(i+1, j) + lastX(i, j+1) + lastX(i+1, j+1)) / 4;
		
	}

	void txt(int step, string name)
	{
		MyWriteFile f = MyWriteFile(FS::getFN(FS::result, name, "txt", step), false);
		f.setSep("");
		for(int i=0; i<X->m; i++)
		{
			int k = X->n;
			for(int j=0; j<k; j++)
			{
				f.write((*X)(i,j));
				if(j < k-1) f.write(";");
			}
			f.newLine();
		}
		f.close();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PDE2);
		for(int i=0; i<m*n; i++)
			f.write(X->v[i]);
	}

	void load(MyReadFile &f)
	{		
		for(int i=0; i<m*n; i++) 
			f.read(X->v[i]);
		lastX.copy(*X);
	}

	DBL getTotalCon()
	{
		DBL con = 0;
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				con += (*X)(i,j);
		return con * dx * dx;
	}

private:
	void copyX()
	{
		Mat<DBL> *p = X;
		X = Xnew;
		Xnew = p;
		Xnew->init(initVal);

		for(int j=0; j<n; j++)
			(*Xnew)(m-1,j) = 0;
	}

	//implicit in x direction
	void setAx(TDM &Ax, int j, Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{
		//left D bc
		Ax.a[0] = 0;
		Ax.b[0] = 1;
		Ax.c[0] = 0;

		//right D bc
		//Ax.a[m-1] = 0;
		//Ax.b[m-1] = 1;
		//Ax.c[m-1] = 0;

		//left N bc
		/*
		int i=0;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2 * f(T(i,j)) - t2x * V(i,j).x;
			Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
		}
		else
		{
			Ax.a[i] = 0;
			Ax.b[i] = at2x2 + 1. + t2 * f(T(i,j));
			Ax.c[i] = - at2x2;
		}
		*/
		
		//right N boundary
		int i = m-1;
		if(V(i,j).x < 0)
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
			Ax.c[i] = 0;
		}
		else if(V(i,j).x > 0)
		{
			Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
			Ax.b[i] = at2x2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) + t2x * V(i,j).x ;
			Ax.c[i] = 0;
		}
		else
		{
			Ax.a[i] = - at2x2;
			Ax.b[i] = at2x2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
			Ax.c[i] = 0;
		}

		for(int i=1; i<m-1; i++)
		{
			if(V(i,j).x < 0)
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) - t2x * V(i,j).x;
				Ax.c[i] = - at2x2 + t2x * V(i+1,j).x;
			}
			else if(V(i,j).x > 0)
			{
				Ax.a[i] = - at2x2 - t2x * V(i-1,j).x;
				Ax.b[i] = atx2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) + t2x * V(i,j).x ;
				Ax.c[i] = - at2x2;
			}
			else
			{
				Ax.a[i] = - at2x2;
				Ax.b[i] = atx2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
				Ax.c[i] = - at2x2;
			}
		}
	}

	//implicit in y direction
	void setAy(TDM &Ay, int i, Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{	
		//bottom D bc
		//int j=0;
		//Ay.a[j] = 0;
		//Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j));
		//Ay.c[j] = - at2y2;

		//top D boundary
		//j = n-1;
		//Ay.a[j] = - at2y2;
		//Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j));
		//Ay.c[j] = 0;
			
		//bottom N boundary
		int j=0;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) - t2y * V(i,j).y;
			Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
			Ay.c[j] = - at2y2;
		}
		else
		{
			Ay.a[j] = 0;
			Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
			Ay.c[j] = - at2y2;
		}

		//top N bc
		j = n-1;
		if(V(i,j).y < 0)
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
			Ay.c[j] = 0;
		}
		else if(V(i,j).y > 0)
		{
			Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
			Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) + t2y * V(i,j).y;
			Ay.c[j] = 0;
		}
		else
		{
			Ay.a[j] = - at2y2;
			Ay.b[j] = at2y2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
			Ay.c[j] = 0;
		}

		for(int j=1; j<n-1; j++)
		{
			if(V(i,j).y < 0)
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) - t2y * V(i,j).y;
				Ay.c[j] = - at2y2 + t2y * V(i,j+1).y;
			}
			else if(V(i,j).y > 0)
			{
				Ay.a[j] = - at2y2 - t2y * V(i,j-1).y;
				Ay.b[j] = aty2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j)) + t2y * V(i,j).y;
				Ay.c[j] = - at2y2;
			}
			else
			{
				Ay.a[j] = - at2y2;
				Ay.b[j] = aty2 + 1. + t2 * f(T(i,j), Fg(i,j), Fp(i,j));
				Ay.c[j] = - at2y2;
			}
		}
	}

	//explicit in y direction
	void setBx(Mat<VecF2> &V)
	{
		//left D bc
		int i0 = 1;
		//right N bc
		int iLim = m;

		for(int j=0; j<n; j++)
		{
			//left D bc
			B(0,j) = initVal;
			//right D bc
			//B(m-1,j) = 0;
		}

		//bottom N boundary
		int j = 0;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)
					+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j).y) * (*X)(i,j) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j)
					+ at2y2 * (*X)(i,j+1);
			else
				B(i,j) = at2y2 * (*X)(i,j) 
					+ (1. - aty2) * (*X)(i,j)
					+ at2y2 * (*X)(i,j+1);
		}

		//top N boundary
		j = n-1;
		for(int i=i0; i<iLim; i++)
		{
			if(V(i,j).y < 0)
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)
					+ (at2y2 - t2y * V(i,j).y) * (*X)(i,j); 
			else if(V(i,j).y > 0)
				B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
					+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j)
					+ at2y2 * (*X)(i,j);
			else
				B(i,j) = at2y2 * (*X)(i,j-1) 
					+ (1. - aty2) * (*X)(i,j)
					+ at2y2 * (*X)(i,j);
		}
		
		#pragma omp parallel for
		for(int j=1; j<n-1; j++)
		{
			for(int i=i0; i<iLim; i++)
			{
				if(V(i,j).y < 0)
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2 + t2y * V(i,j).y) * (*X)(i,j)
						+ (at2y2 - t2y * V(i,j+1).y) * (*X)(i,j+1); 
				else if(V(i,j).y > 0)
					B(i,j) = (at2y2 + t2y * V(i,j-1).y) * (*X)(i,j-1) 
						+ (1. - aty2 - t2y * V(i,j).y) * (*X)(i,j)
						+ at2y2 * (*X)(i,j+1);
				else
					B(i,j) = at2y2 * (*X)(i,j-1) 
						+ (1. - aty2) * (*X)(i,j)
						+ at2y2 * (*X)(i,j+1);
			}
		}		
	}

	//explicit in x direction
	void setBy(Mat<VecF2> &V)
	{
		for(int j=0; j<n; j++)
		{
			//left D bc
			B(0,j) = initVal;
			//right D bc
			//B(m-1,j) = 0;
		}

		int i = m-1;
		for(int j=0; j<n; j++)
		{
			if(V(i,j).x < 0)
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j)
					+ (at2x2 - t2x * V(i,j).x) * (*X)(i,j); 
			else if(V(i,j).x > 0)
				B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
					+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j)
					+ at2x2 * (*X)(i,j);
			else
				B(i,j) = at2x2 * (*X)(i-1,j) 
					+ (1. - atx2) * (*X)(i,j)
					+ at2x2 * (*X)(i,j);
		}

		#pragma omp parallel for
		for(int i=1; i<m-1; i++)
		{
			for(int j=0; j<n; j++)
			{
				if(V(i,j).x < 0)
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2 + t2x * V(i,j).x) * (*X)(i,j)
						+ (at2x2 - t2x * V(i+1,j).x) * (*X)(i+1,j); 
				else if(V(i,j).x > 0)
					B(i,j) = (at2x2 + t2x * V(i-1,j).x) * (*X)(i-1,j) 
						+ (1. - atx2 - t2x * V(i,j).x) * (*X)(i,j)
						+ at2x2 * (*X)(i+1,j);
				else
				{
					B(i,j) = at2x2 * (*X)(i-1,j) 
						+ (1. - atx2) * (*X)(i,j)
						+ at2x2 * (*X)(i+1,j);
				}
			}
		}
	}

	//implicit in x, explicit in y direction
	void thomasX(Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{
		//0 <= j < n-1
		#pragma omp parallel for
		for(int j=0; j<n; j++)
		{
			vector<DBL> c = vector<DBL>(m-1, 0);
			vector<DBL> d = vector<DBL>(m, 0);
			
			TDM Ax = TDM(m);
			setAx(Ax, j, V, T, Fg, Fp);

			c[0] = Ax.c[0] / Ax.b[0];
			for(int i=1; i<m-1; i++)
				c[i] = Ax.c[i] / (Ax.b[i] - c[i-1] * Ax.a[i]);

			d[0] = B(0,j) / Ax.b[0];
			for(int i=1; i<m; i++)
			{
				d[i] = (B(i,j) - d[i-1] * Ax.a[i]) / (Ax.b[i] - c[i-1] * Ax.a[i]);
			}

			(*Xnew)(m-1, j) = d[m-1];
			for(int i=m-2; i>=0; i--)
			{
				(*Xnew)(i,j) = d[i] - c[i]*(*Xnew)(i+1,j);;	
			}

		}	
	}

	//implicit in y, explicit in x direction
	void thomasY(Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{
		//0 < i < m-1
		#pragma omp parallel for
		for(int i=0; i<m; i++)	//Neumann
		//for(int i=1; i<m-1; i++)	//Dirichlet
		{
			vector<DBL> c = vector<DBL>(n-1, 0);
			vector<DBL> d = vector<DBL>(n, 0);

			TDM Ay = TDM(n);
			setAy(Ay, i, V, T, Fg, Fp);
			
			c[0] = Ay.c[0] / Ay.b[0];
			for(int j=1; j<n-1; j++)
				c[j] = Ay.c[j] / (Ay.b[j] - c[j-1] * Ay.a[j]);

			d[0] = B(i,0) / Ay.b[0];
			for(int j=1; j<n; j++)
			{
				d[j] = (B(i,j) - d[j-1] * Ay.a[j]) / (Ay.b[j] - c[j-1] * Ay.a[j]);
			}

			(*Xnew)(i, n-1) = d[n-1];
			for(int j=n-2; j>=0; j--)
			{
				(*Xnew)(i,j) = d[j] - c[j]*(*Xnew)(i,j+1);	
			}
		}	
	}

	//implicit in x, explicit in y direction
	void halfStepX(Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{
		setBx(V);
		thomasX(V, T, Fg, Fp);
		copyX();
	}

	//implicit in y, explicit in x direction
	void halfStepY(Mat<VecF2> &V, Mat<DBL> &T, Mat<DBL> &Fg, Mat<DBL> &Fp)
	{
		setBy(V);
		thomasY(V, T, Fg, Fp);
		copyX();
	}

	inline DBL f(DBL T, DBL Fg, DBL Fp)
	{
		return beta3 * T * Fg * (1. - Fp); 
	}
};






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//fibrin polymer
// dFg/dt = f(T) Fg
// f(T) = k3 T
class Pde3
{
public:
	DBL dx;
	DBL dy;
	DBL dt;
	Mat<DBL> X;
	int m;
	int n;
	DBL beta3;

	Pde3()
	{
		X = Mat<DBL>();
	}

	Pde3(Setup &s)
	{
		size = s.spc.size;
		dx = s.fib.dx;
		dy = s.fib.dx;
		dt = s.fib.dt;
		beta3 = s.fib.beta3;

		m = size.x / dx + 1;
		n = size.y / dx + 1;
		X = Mat<DBL>(m,n);
		X.init(0);
	}

	void paramHotChange(Setup &s)
	{
		dt = s.fib.dt;
		beta3 = s.fib.k3;
	}

	void nextStep(Mat<DBL> &T, Mat<DBL> &Fg)
	{
		int k = m*n;
		#pragma omp parallel
		for(int i=0; i<k; i++)
			X.v[i] = X.v[i] + dt*f(T.v[i], Fg.v[i], X.v[i]);
	}

	void txt(int step, string name)
	{
		MyWriteFile f = MyWriteFile(FS::getFN(FS::result, name, "txt", step), false);
		f.setSep("");
		for(int i=0; i<m; i++)
		{
			for(int j=0; j<n; j++)
			{
				f.write(X(i,j));
				if(j < n-1) f.write(";");
			}
			f.newLine();
		}
		f.close();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PDE3);
		for(int i=0; i<m*n; i++)
			f.write(X.v[i]);
	}

	void load(MyReadFile &f)
	{		
		for(int i=0; i<m*n; i++) 
			f.read(X.v[i]);
	}


	inline DBL getU(VecF2 &p)
	{		
		int i = p.x / dx;
		int j = p.y / dy;
		if(i==m-1) i--;
		if(j==n-1) j--;
		return bilin(p.x, p.y, dx*i, dx*(i+1), dy*j, dy*(j+1), X(i, j), X(i, j+1), X(i+1, j), X(i+1, j+1));
	}

	DBL getTotalCon()
	{
		DBL con = 0;
		for(int i=0; i<m; i++)
			for(int j=0; j<n; j++)
				con += X(i,j);
		return con * dx * dx;
	}

private:
	VecF2 size;

	inline DBL f(DBL T, DBL Fg, DBL Fp)
	{
		return beta3 * T * Fg * (1. - Fp); 
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//thrombin explicit
// dT/dt = alpha (d^2T/dx^2) - V.(Del T) + f(T)(C0-T) - gamma T
// f(T) = k1 (T^2)/(T0+T)
class Pde1E
{
public:
	DBL dx;
	DBL dy;
	DBL dt;
	Mat<DBL> *X;
	Mat<DBL> *Xnew;
	int m;
	int n;
	DBL alpha3;
	DBL k1;
	DBL gamma3;
	DBL C0;
	DBL T0;
	Mat<DBL> B;
	DBL t2, t2x, at2x2, atx2, t2y, at2y2, aty2;
	VecF2 size;
	vector<DBL> posX;
	Mat<VecF2> V;
	Mat<DBL> lastX;

	Pde1E()
	{
		X		= NULL;
		Xnew	= NULL;
		V		= Mat<VecF2>();
		lastX	= Mat<DBL>();
		B		= Mat<DBL>();
	}

	Pde1E(Setup &s)
	{
		if(!s.fib.enabled) return;
		size = s.spc.size;
		alpha3 = s.fib.alpha3;
		k1 = s.fib.k1;
		gamma3 = s.fib.gamma3;
		C0 = s.fib.C0;
		T0 = s.fib.T0;
		dx = s.fib.dx;
		dy = s.fib.dx;
		dt = s.fib.dt;
		
		m = size.x / dx + 1;
		n = size.y / dx + 1;
		X = new Mat<DBL>(m, n);
		B = Mat<DBL>(m, n);
		Xnew = new Mat<DBL>(m, n);
		Mat<DBL> B(m, n);
		V = Mat<VecF2>(m,n);
		lastX = Mat<DBL>(m,n);

		V.init(VecF20);
		//initial conditions
		X->init(0);
		Xnew->init(0);

		for(int i=0; i<m; i++)	
		{
			DBL x = i*dx;		

			if(x >= s.fib.woundBeg && x <= s.fib.woundEnd)
				(*X)(i,0) = s.fib.woundVal;
		}

		t2 = dt / 2;
		t2x = t2 / dx;
		at2x2 = alpha3 * t2x / dx;
		atx2 = at2x2 * 2;
		t2y = t2 / dy;
		at2y2 = alpha3 * t2y / dy;
		aty2 = at2y2 * 2;

		lastX.copy(*X);
	}

	void nextStep()
	{
		DBL atx2 = alpha3 * dt / pow(dx,2);
		DBL aty2 = alpha3 * dt / pow(dy,2);

		for(int i=0; i<m; i++)
		{
			for(int j=0; j<m; j++)
			{
				DBL vx = 0;
				if(V(i,j).x < 0)
				{
					(*Xnew)(i,j);
				}
				else if(V(i,j).x > 0)
				{
				}

				DBL vy = 0;
				if(V(i,j).y < 0)
				{
				}
				else if(V(i,j).y > 0)
				{
				}				

				(*Xnew)(i,j) = (*Xnew)(i,j) 
					+ atx2 * ((*Xnew)(i-1,j) - (*Xnew)(i,j) * 2 + (*Xnew)(i+1,j))
					+ aty2 * ((*Xnew)(i,j-1) - (*Xnew)(i,j) * 2 + (*Xnew)(i,j+1))
					- vx * dt
					- vy * dt
					+ f((*Xnew)(i,j)) * dt
					- gamma3 * (*Xnew)(i,j) * dt;
			}
		}

		copyX();
	}

	void free()
	{
		if(X != NULL) 
		{
			delete X;
			X = NULL;
		}
		if(Xnew != NULL) 
		{
			delete Xnew;
			Xnew = NULL;
		}
	}

	void saveWave()
	{
		int y = n / 2;
		for(int i=0; i<m; i++)
		{
			if(((*X)(i,y)-0.5)*((*X)(i+1,y)-0.5) < 0)
			{
				DBL aaa = (*X)(i,y)-(*X)(i+1,y);
				posX.push_back(((*X)(i,y)/aaa-0.5/aaa+i)*dx);
				break;
			}
		}
	}


	//returns the concentration on given coordinates p
	inline DBL getU(VecF2 &p)
	{	
		int i = p.x / dx;
		int j = p.y / dy;
		return bilin(p.x, p.y, dx*i, dx*(i+1), dy*j, dy*(j+1), lastX(i, j), lastX(i+1, j), lastX(i, j+1), lastX(i+1, j+1));
		//return (lastX(i, j) + lastX(i+1, j) + lastX(i, j+1) + lastX(i+1, j+1)) / 4;
		
	}

	void update(DBL coeff)
	{
		for(int j=0; j<n; j++)
		{
			DBL y = dy * j;
			for(int i=0; i<m; i++)
			{
				V(i,j) = VecF2(coeff * y*(size.y - y) * 4. / pow(size.y, 2), 0);
				//V(i,j) = VecF2(coeff,0);
			}
		}

		lastX.copy(*X);
	}

	void update(Mat<Stat> &M, DBL ignoreVal, CProgressCtrl &p, CStatic &desc)
	{
		DBL da = size.x / M.m;
		DBL db = size.y / M.n;
		DBL dah = da / 2;
		DBL dbh = db / 2;
		DBL daL = size.x - da;
		DBL dbL = size.y - db;

		int i1=0, i2=0, j1=0, j2=0;
		DBL dxL=0, dx1=0, dx2=0, dyL=0, dy1=0, dy2=0;

		desc.SetWindowTextA("Adjusting velocity profile for PDE.");
		for(int i=0; i<m; i++)
		{
			bool x0 = false;
			DBL x = dx * i - dah;

			if(x < 0)
			{
				x0 = true;
				i1 = 0;
			}
			else if(x >= daL)
			{
				x0 = true;
				i1 = M.m-1;
			}
			else
			{
				i1 = x / da;
				i2 = i1 + 1;
				dxL = da;
				dx1 = dx * i - (i1 * da + dah);
				dx2 = i2 * da + dah - dx * i;
			} 

			for(int j=0; j<n; j++)
			{
				DBL y = dy * j - dbh;

				if(y < 0)
				{
					dy1 = dy * j;
					if(x0) 
						V(i,j) = M(i1,0).v * dy1 / dbh;
					else 
						V(i,j) = ((M(i1,0).v * dx2 + M(i2,0).v * dx1)) * dy1 / (da * dbh);
				}
				else if(y >= dbL)
				{
					dy2 = size.y - dy * j;
					if(x0)
						V(i,j) = M(i1, M.n-1).v * dy2 / dbh;
					else
						V(i,j) = (M(i1, M.n-1).v * dx2 + M(i2, M.n-1).v * dx1) * dy2 / (da * dbh);
				}
				else
				{
					j1 = y / db;
					j2 = j1 + 1;
					dyL = db;
					dy1 = dy * j - (j1 * db + dbh);
					dy2 = j2 * db + dbh - dy * j;

					if(x0)
						V(i,j) = (M(i1, j1).v * dy2 + M(i1, j1).v * dy1) / db;
					else
						V(i,j) = (M(i1, j1).v * dx2 * dy2 
							+ M(i2, j1).v * dx1 * dy2
							+ M(i1, j2).v * dx2 * dy1
							+ M(i2, j2).v * dx1 * dy1)
							/ (da * db);
				} 
			}
		}

		p.StepIt();

		DBL sigma = 0.1;
		int N = min(V.m, V.n);
		Vec<DBL> f = gaussFilt(N, sigma * N/2);

		desc.SetWindowTextA("Smoothing velocity profile.");
		convAll(V, f, p);

		//ignoring small velocities
		for(int i=0; i<m; i++)
		{
			for(int j=0; j<n; j++)
			{
				if(fabs(V(i,j).x) < ignoreVal) V(i,j).x = 0;
				if(fabs(V(i,j).y) < ignoreVal) V(i,j).y = 0;
			}
		}

		lastX.copy(*X);
	}

	void drawSmoothed(Mat<DBL>& Vx, Mat<DBL>& Vy)
	{
		VecF2 rx = VecF2(0, size.x * MU::len);
		VecF2 ry = VecF2(0, size.y * MU::len);

		Graph2 g = Graph2();
		g.create(FS::result, "vSmooth", -1, 1, 2);
		g.setLabel("y [m]", "x [m]", "\\v_x [m/s]"); g.setData(Vx); g.setRange(rx, ry); g.plot(VecL3(1,2,0));
		g.setLabel("y [m]", "x [m]", "\\v_y [m/s]"); g.setData(Vy); g.setRange(rx, ry); g.plot(VecL3(1,2,1));
		g.draw();
	}

	void txt(int step, string name)
	{
		MyWriteFile f = MyWriteFile(FS::getFN(FS::result, name, "txt", step), false);
		f.setSep("");
		for(int i=0; i<X->m; i++)
		{
			int k = X->n;
			for(int j=0; j<k; j++)
			{
				f.write((*X)(i,j));
				if(j < k-1) f.write(";");
			}
			f.newLine();
		}
		f.close();
	}

	void save(MyWriteFile &f)
	{
		f.writeTag(Const::PDE1);
		for(int i=0; i<m*n; i++)
			f.write(X->v[i]);
	}

	void load(MyReadFile &f)
	{		
		for(int i=0; i<m*n; i++) 
			f.read(X->v[i]);
		lastX.copy(*X);
	}

private:
	void copyX()
	{
		Mat<DBL> *p = X;
		X = Xnew;
		Xnew = p;
		Xnew->init(0);
	}

	inline DBL f(DBL T)
	{
		return k1 * pow(T, 2) / (T0 + T);
	}
};






#endif
