#ifndef DATA_H
#define	DATA_H

namespace MU
{
	static DBL len	= 0;
	static DBL area = 0;
	static DBL vol	= 0;
	static DBL time	= 0;
	static DBL mass	= 0;
	static DBL dens	= 0;
	static DBL visc	= 0;
	static DBL vel	= 0;
	static DBL acc	= 0;


	static void init(int kL, int kM, int kT)
	{
		len		= pow(10., (DBL)(kL));	
		area	= pow(10., (DBL)(2*kL));
		vol		= pow(10., (DBL)(3*kL));
		time	= pow(10., (DBL)kT);
		mass	= pow(10., (DBL)kM);
		dens	= pow(10., (DBL)(kM - 3*kL));
		vel		= pow(10., (DBL)(kL - kT));
		acc		= pow(10., (DBL)(kL - 2*kT));
		visc	= pow(10., (DBL)(kM - kL - kT));
	}
}

#endif