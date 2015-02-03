#ifndef CODES_H
#define	CODES_H

//for method
static const int M_MD			= 201;	// Molecular Dynamics
static const int M_DPD			= 202;	// Dissipative Particle Dynamics
static const int M_MIX			= 203;	// DPD with MD conservative force

//sim option
static const int SIM_DPD		= 301;
static const int SIM_TRANS		= 302;
static const int SIM_PDE		= 303;

//platelet force type
static const int PL_LIN			= 0;
static const int PL_STEP		= 1;
static const int PL_HYB			= 2;

//for PDE model type
static const int HYB1			= 0;	//fibrin with CRDE
static const int HYB2			= 1;	//thrombin with CRDE and fibrin

//for file draw
static const int DRAW_PLOT = 1;
static const int DRAW_SURF = 2;
static const int DRAW_DENS = 3;
static const int DRAW_VECT = 4;

#endif