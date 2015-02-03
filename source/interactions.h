#ifndef INTERACTIONS_H
#define	INTERACTIONS_H



class Interact
{
public:
	static void atomAtom(Setup &s, Atom *a1, Atom *a2, VecF2 &displacement, bool oneWay, PPL &platPairs, PdeSys &ps)
	{
		VecF2 p=a1->pos;
		VecF2 q=a2->pos;
		DBL pp=ps.pdeFp.getU(p);
		DBL qq=ps.pdeFp.getU(q);
		VecF2 r = a1->pos - displacement - a2->pos;
		DBL h = MODUL(r);

		if(h < s.met.rc)
		{
			bool counterFlow1 = !a1->prep && a2->prep;
			bool counterFlow2 = a1->prep && !a2->prep;

			Fn2(s, a1, a2, r, h, !oneWay && !counterFlow2, !counterFlow1, pp, qq);
			
			if(a1->isSPlt() && a2->isSPlt() && s.plt.enabled)
			{
				DBL h0 = s.atom.r * 2;
				if(h < h0)
				{
					VecF2 tmpF = Interact::pltFnMinus(r, h, h0);
						
					if(!oneWay && !counterFlow2) a1->FPL1M += tmpF;
					if(!counterFlow1) a2->FPL1M -= tmpF;
				}

				if(!(a1->prep || a2->prep) && (a1->pos.x > s.spc.sizePrepX && a2->pos.x > s.spc.sizePrepX))
				{
					if(!((a1->oldPlt && a2->oldPlt) || (a1->fix && a2->fix) || (a1->oldPlt && !a1->fix) || (a2->oldPlt && !a2->fix)))
					{
						VecF2 r2 = a1->pos - a2->pos;
						DBL h2 = MODUL(r2);

						h0 = h0 * s.plt.dc;
						if(h2 < h0)
						{
							#pragma omp critical
							{
								platPairs.add(*a1, *a2, h0);
							}
						}
					}
				}
			}
			
			VecF2 p1 = a1->pos;
			VecF2 p2 = a2->pos + displacement;
			VecF2 v1 = a1->v;
			VecF2 v2 = a2->v;

			bool wT = false;
			bool wB = false;

			if(p1.y + p2.y < s.met.rc) wB = true;
			else if(s.spc.size.y * 2 - p1.y - p2.y < s.met.rc) wT = true;

			if(wB || wT)
			{ 
				VecF2 p2m = wB ? VecF2(p2.x, -p2.y) : VecF2(p2.x, s.spc.size.y * 2 - p2.y);
				//mirror(s, a1, a2, p1, p2m, s.met.rc, NONE, NOSLIP, oneWay);
				mirror(s, a1, a2, p1, p2m, s.met.rc, oneWay, ps);
			}
		}
	}

	static VecF2 pltFnPlus(VecF2 &r, DBL h, DBL d)
	{
		if(h > d) return (1. - h / d) * (r/h); else return VecF20;
	}

	static void selfMirror(Setup &s, AtomList &al, PdeSys &ps)
	{
		int n = al.getSize();
		#pragma omp parallel for
		for(int i=0; i<n; i++)
			if(!al.obj[i].free) atomSelfMirror(s, al.obj[i], ps);
	}

private:
//	static const int NONE = 0;		// no influence
//	static const int NOSLIP = 1;	// rebound 
//	static const int FILTER = 2;	// continuous 

	//static void mirror(Setup &s, Atom *a1, Atom *a2, VecF2 p1, VecF2 p2m, DBL h0, int mX, int mY, bool oneWay)
	static void mirror(Setup &s, Atom *a1, Atom *a2, VecF2 p1, VecF2 p2m, DBL h0, bool oneWay, PdeSys &ps)
	{

		VecF2 p=a1->pos;
		VecF2 q=a2->pos;
		DBL pp=ps.pdeFp.getU(p);
		DBL qq=ps.pdeFp.getU(q);
		VecF2 rm = p1 - p2m;
		DBL hm = MODUL(rm);
		if(hm < h0)
		{
			//bool vxChange = (mX == FILTER && mY == NONE) || (mX == NOSLIP && mY == FILTER) || (mX == NOSLIP && mY == NOSLIP);
			//bool vyChange = (mX == NONE && mY == FILTER) || (mX == FILTER && mY == NOSLIP) || (mX == NOSLIP && mY == NOSLIP);
			//bool pxChange = mX == NONE;
			//bool pyChange = mY == NONE;
			
			bool vxChange = false;
			bool vyChange = false;
			bool pxChange = true;
			bool pyChange = false;

			VecF2 vm = VecF20;
			vm.x = vxChange ? a2->v.x : -a2->v.x;
			vm.y = vyChange ? a2->v.y : -a2->v.y;

			if(!oneWay) Fn1(s, a1, rm, hm, vm, pp);

			if(pxChange) rm.x = -rm.x;
			if(pyChange) rm.y = -rm.y;
			vm.x = vxChange ? a1->v.x : -a1->v.x;
			vm.y = vyChange ? a1->v.y : -a1->v.y;

			Fn1(s, a2, rm, hm, vm, qq);
		}
	}

	static void atomSelfMirror(Setup &s, Atom &a, PdeSys &ps)
	{
		bool wallTB = false;
		bool wallLR = false;
		DBL h0 = s.met.rc;
		VecF2 p = a.pos;
		VecF2 v = a.v;
		DBL q=ps.pdeFp.getU(p);
		if(p.y < s.met.rc)
		{
			DBL t = s.met.rc - p.y;
			if(p.y < t)
			{
				VecF2 pm = VecF2(p.x, -p.y);
				VecF2 rm = p - pm;
				DBL hm = MODUL(rm);
				if(hm < h0)
				{
					Fn1(s, &a, rm, hm, -a.v, q);

					wallTB = true;
				}
			}
		}
			
		if(p.y > s.spc.size.y - s.met.rc)
		{
			DBL	t = s.met.rc - (s.spc.size.y - p.y);
			if(p.y > s.spc.size.y - t)
			{
				VecF2 pm = VecF2(p.x, s.spc.size.y * 2 - p.y);
				VecF2 rm = p - pm;
				DBL hm = MODUL(rm);
				if(hm < h0)
				{
					Fn1(s, &a, rm, hm, -a.v, q);

					wallTB = true;
				}
			}
		}

		if(s.spc.prepArea)
		{
			if(p.x > s.spc.size.x - s.met.rc)
			{
				DBL	t = s.met.rc - (s.spc.size.x - p.y);
				if(p.y > s.spc.size.x - t)
				{
					VecF2 pm = VecF2(s.spc.size.x * 2 - p.x, p.y);
					VecF2 rm = p - pm;
					DBL hm = MODUL(rm);
					if(hm < h0)
					{
						Fn1(s, &a, rm, hm, VecF2(a.v.x, -a.v.y), q);

						wallLR = true;
					}
				}
			}
		}

		if(wallTB && wallLR)
		{
			VecF2 pm = VecF20;
			if(p.y < s.met.rc)
				pm = VecF2(s.spc.size.x * 2 - p.x, -p.y);
			else if(p.y > s.spc.size.y - s.met.rc)
				pm = VecF2(s.spc.size.x * 2 - p.x, s.spc.size.y * 2 - p.y);
				
			VecF2 rm = p - pm;
			DBL hm = MODUL(rm);
			if(hm < h0)
			{
				Fn1(s, &a, rm, hm, VecF2(-a.v.x, a.v.y), q);
			}
		}
	} 

	static VecF2 pltFnMinus(VecF2 &r, DBL h, DBL d)
	{
		if(h < d) return (1. - h / d) * (r/h); else return VecF20;
	}


	static void Fn1(Setup &s, Atom *a1, VecF2 &r, DBL h, VecF2 v, DBL p)
	{
		DBL tmp = 1. - h / s.met.rc;
		DBL oR = s.met.p == 1 ? tmp : pow(tmp, s.met.p);
		DBL oD = pow(oR, 2);
		VecF2 r1 = r / h;
		VecF4 r1r1(r1.x*r1.x, r1.x*r1.y, r1.y*r1.x, r1.y*r1.y); 
		VecF2 pp=VecF2(p,p);
		VecF2 FC = tmp * r1;
		VecF2 FR = oR * Random(-1, 1) * r1;
		VecF2 FD1 = oD * (v * r1) * r1;
		VecF4 FD2 = oD * r1r1;

		a1->FC += FC;
		a1->FR += FR;
		a1->FD1 += FD1;
		a1->FD2 += FD2;
	}

	static void Fn2(Setup &s, Atom *a1, Atom *a2, VecF2 &r, DBL h, bool f1, bool f2, DBL p, DBL q)
	{
		DBL tmp = 1. - h / s.met.rc;
		DBL oR = s.met.p == 1 ? tmp : pow(tmp, s.met.p);
		DBL oD = pow(oR, 2);
		VecF2 pp=VecF2(p,p);
		VecF2 qq=VecF2(q,q);
		VecF2 r1 = r / h;
		VecF4 r1r1(r1.x*r1.x, r1.x*r1.y, r1.y*r1.x, r1.y*r1.y); 

		VecF2 FC = tmp * r1;
		VecF2 FR = oR * Random(-1, 1) * r1;
		VecF2 FD1 = oD * (a2->v * r1) * r1;
		VecF4 FD2 = oD * r1r1;

		if(f1)
		{

			a1->FC += FC;
			a1->FR += FR;
			a1->FD1 += FD1;
			a1->FD2 += FD2;
		}

		if(f2)
		{

			r1 = -r1;
			FD1 = oD * (a1->v * r1) * r1;
			a2->FC -= FC;
			a2->FR -= FR;
			a2->FD1 += FD1;
			a2->FD2 += FD2;
		}
	}

};

#endif