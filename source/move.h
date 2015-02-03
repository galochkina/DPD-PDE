#ifndef MOVE_H
#define	MOVE_H

class Move
{
public: 

	static void move(Setup &s, Atom &a, AtomList &al, bool plt, DBL dt)
	{
		if(plt && !a.isSPlt()) return;
		if(!plt && a.isSPlt()) return;
		if(a.fix) return;
		VecF2 posNew = a.pos;
		VecF2 Fprep = VecF2(120*s.atom.v_maxPF/pow(50,2)*300,0);
		//VecF2 Fprep = VecF2(80000,0);
		VecF2 F1;

		VecF2 G = (!s.spc.prepArea || a.prep) ? s.met.G : VecF20;
		if (a.pos.x <= s.spc.sizePrepX){
		F1 = ((a.FC * s.met.A + a.FR * s.met.sigmaSqrtdt + a.FD1 * s.met.gamma + (a.FPL1M + a.FPL1) * s.plt.F2weak + a.FPL2 * s.plt.F2strong + a.FPL3 * s.plt.F2fib + Fprep) / s.atom.m + G) * dt + a.v;
		}
		else {
        F1 = ((a.FC * s.met.A + a.FR * s.met.sigmaSqrtdt + a.FD1 * s.met.gamma + (a.FPL1M + a.FPL1) * s.plt.F2weak + a.FPL2 * s.plt.F2strong + a.FPL3 * s.plt.F2fib) / s.atom.m + G ) * dt + a.v;
		}
		VecF4 FD2 = a.FD2 * s.met.gamma * dt / s.atom.m;
		FD2.x += 1;
		FD2.w += 1;
		DBL kk = 1. / (FD2.x * FD2.w - FD2.y * FD2.z);
		VecF4 Inv(kk * FD2.w, -kk * FD2.y, -kk * FD2.z, kk * FD2.x);
		a.v.x = F1.x * Inv.x + F1.y * Inv.z;
		a.v.y = F1.x * Inv.y + F1.y * Inv.w;

		posNew = a.pos + dt * a.v;

		if(a.isSPlt()) a.resetFplt();

		if(s.met.maxV != 0)
		{
			DBL modul = MODUL(a.v);
			if(modul > s.met.maxV) 
				a.v = a.v / modul * s.met.maxV;
		}

		if(posNew.y < 0 || posNew.y >= s.spc.size.y) 
		{
			DBL dist = min(a.pos.y, s.spc.size.y - a.pos.y);
			DBL time = fabs(dist / a.v.y);
			posNew.x = a.pos.x + a.v.x * time;
			posNew.y = a.pos.y + a.v.y * (time*2 - dt);
			a.v.x = 0;
			a.v.y = - a.v.y;
			a.acc.x = 0;
			a.acc.y = 0;
		}

		    DBL more=s.spc.size.x - 1;
			DBL prepHalf = s.spc.sizePrepX/2;
			//if(a.pos.x <= prepHalf && posNew.x > s.spc.sizePrepX) a.copied = false;
			if(posNew.x <= s.spc.size.x) a.copied = false;
			if(posNew.x < 0) 
			{
				posNew.x = mod(posNew.x, s.spc.size.x);
			}
			if(a.pos.x<= s.spc.sizePrepX && posNew.x > s.spc.sizePrepX) 
			{
				DBL r=rand()%100;
				if(r<=1.05){
					a.type=1;}
			}
			else if(posNew.x > s.spc.size.x)
			{

			
			bool plBool = false;
				//delete particle
				//Atom aNew(a);
				//aNew.prep = false;
				//aNew.copied = false;
				//aNew.pos = posNew;

				

						//a.copied = true;
						//new copy of particle

					//}
				//if(a.copied=true) al.remove(a.id);
					//periodic movement
				posNew.x = mod(posNew.x, s.spc.size.x);
				
				//al.addLater(aNew);
				
				//al.remove(aNew.id);
				

			}
			//else if(posNew.x > more)
			//{  
				//if(a.prep)
				//{
					//if(!a.copied)
					//{
						//Atom aNew(a);
						//aNew.prep = false;
						//aNew.copied = false;
						//aNew.pos = posNew;

						//al.addLater(aNew);

						//a.copied = true;
						//new copy of particle

					//}
					//if(a.copied=true) al.remove(a.id);
					//periodic movement
					//posNew.x = mod(posNew.x, s.spc.size.x);
				//} 
			//}	
		
		

		//if(s.plt.constInflow && a.isSPlt() && a.prep) posNew.y = a.pos.y;

		a.pos = posNew;

	}
};

#endif