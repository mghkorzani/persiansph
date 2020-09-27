/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "./Source/Domain.h"

using std::cout;
using std::endl;
using std::ifstream;

int check = 0;
double DampF,DampS,Cs,u;
double t = 0.5;
double tim = 0.0;
void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<t)
    #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3) domi.Particles[i]->a -= DampS * domi.Particles[i]->v;

	}
	tim = domi.Time;

	if (domi.Time>0.0 && check==0)
	{
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
			if (domi.Particles[i]->ID == 6 && domi.Particles[i]->Material == 3)
			{
	    		domi.Particles[i]->c 	= 0.0;
	    		domi.Particles[i]->phi	= 15.0/180.0*M_PI;
				domi.Particles[i]->k 	= 10.0;
			}
//			if (domi.Particles[i]->ID == 3 && domi.Particles[i]->Material == 3)
//	    		domi.Particles[i]->c	= 5.0e3;
		}

		check = 1;
	}
}

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (tim<1.0)
	{
		if (position(1)<(0.5*4.5))
			Vel = pow((position(1)/(0.32*4.5)),(1.0/7.0))*u*tim,0.0,0.0;
		else
			Vel = 1.07*u*tim,0.0,0.0;
	}
	else
	{
		if (position(1)<(0.5*4.5))
			Vel = pow((position(1)/(0.32*4.5)),(1.0/7.0))*u,0.0,0.0;
		else
			Vel = 1.07*u,0.0,0.0;
	}

//	Den = 998.21*((1+9.81*(4.5-position(1))/(Cs*Cs)));
	Den = 998.21*pow((1+7.0*9.81*(4.5-position(1))/(Cs*Cs)),(1.0/7.0));
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
		if (position(1)<(0.5*4.5))
			Vel = pow((position(1)/(0.32*4.5)),(1.0/7.0))*u,0.0,0.0;
		else
			Vel = 1.07*u,0.0,0.0;
//	Den = 998.21*((1+9.81*(4.5-position(1))/(Cs*Cs)));
	Den = 998.21*pow((1+7.0*9.81*(4.5-position(1))/(Cs*Cs)),(1.0/7.0));
}

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.KernelType	= 0;
	dom.Nproc		= 24;

//	dom.XSPH		= 0.5;
	dom.Scheme		= 0;
	dom.Gravity		= 0.0,-9.81,0.0;
	dom.GeneralAfter= & UserDamping;

	u					= 0.18;
	dom.BC.InOutFlow	= 3;
	dom.BC.inv			= u,0.0,0.0;
//	dom.BC.allv			= u,0.0,0.0;
//	dom.BC.allDensity	= 998.21;
	dom.BC.inDensity	= 998.21;
	dom.InCon			= & UserInFlowCon;
//	dom.AllCon			= & UserAllFlowCon;

	double xb,yb,h,dx,T;

	dx		= 0.15;
	h		= dx*1.3;
	dom.InitialDist	= dx;

	double rhoF,Mu,CsF,Tf;
	rhoF	= 998.21;
	Mu		= 1.002e-3;
	CsF		= 10.0*sqrt(2.0*9.81*4.5);
	CsF		= 10.0*5.0;
	Tf		= (0.2*h/CsF);
	DampF 	= 0.04*CsF/h;
	dom.VisEq= 1;

	Cs = CsF;


	dom.AddBoxLength(1 ,Vec3_t ( -16.0 - 3.0*dx , -3.0*dx , 0.0 ), 32.0 + 6.0*dx + dx/10.0 , 4.5 + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoF, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->PresEq	= 1;
		dom.Particles[a]->Alpha		= 0.03;
//		dom.Particles[a]->Beta		= 0.05;
		dom.Particles[a]->Mu		= Mu;
		dom.Particles[a]->MuRef		= Mu;
		dom.Particles[a]->Material	= 1;
		dom.Particles[a]->Cs		= CsF;
//		dom.Particles[a]->Shepard	= true;
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (yb<0.0)
		{
			dom.Particles[a]->ID	= 2;
			dom.Particles[a]->NoSlip= true;
			dom.Particles[a]->IsFree= false;
		}
//		if (yb>(-0.5*(xb+2.0)+4.5) && dom.Particles[a]->ID == 1)
//			dom.Particles[a]->ID	= 5;
//		if (yb>(-0.5*(xb-11.0)) && dom.Particles[a]->ID == 1)
//			dom.Particles[a]->ID	= 5;
		if (yb>(-0.346*(xb-11.0)) && dom.Particles[a]->ID == 1)
			dom.Particles[a]->ID	= 5;

		if (dom.Particles[a]->ID==1) dom.Particles[a]->Density  = rhoF*pow((1+7.0*9.81*(4.5-dom.Particles[a]->x(1))/(CsF*CsF)),(1.0/7.0));
//		if (dom.Particles[a]->ID==1) dom.Particles[a]->Density  = rhoF*((1+9.81*(4.5-dom.Particles[a]->x(1))/(CsF*CsF)));
	}

	double K,G,Nu,E,rhoS,CsS,Phi,c,Psi,k,Ts;

	rhoS	= 2038.7;
	Phi		= 25.0;
	Psi		= 0.0;
	k		= 15.83;
	c		= 5.0e3;
	E		= 25.0e6;
	Nu		= 0.3;
	K		= E/(3.0*(1.0-2.0*Nu));
	G		= E/(2.0*(1.0+Nu));
    CsS		= sqrt(E/rhoS);
    CsS		= 200.0;
    Ts		= (0.2*h/CsS);
   	DampS	= 0.05*sqrt(E/(rhoS*h*h));

	dom.AddBoxLength(3 ,Vec3_t ( -16.0 , -3.0*dx , 0.0 ), 32.0 + dx/10.0 , 5.0 + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoS, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==3)
		{
			dom.Particles[a]->n			= 1.0;
			dom.Particles[a]->k			= k;
			dom.Particles[a]->RhoF		= rhoF;
			dom.Particles[a]->Cs		= CsS;
			dom.Particles[a]->G			= G;
			dom.Particles[a]->K			= K;
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Fail		= 3;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
			dom.Particles[a]->TI		= 0.5;
			dom.Particles[a]->TIn		= 2.55;
    		dom.Particles[a]->c			= c;
    		dom.Particles[a]->phi		= Phi/180.0*M_PI;
    		dom.Particles[a]->psi		= Psi/180.0*M_PI;
			xb=dom.Particles[a]->x(0);
			yb=dom.Particles[a]->x(1);
			if (yb<0.0)
			{
				dom.Particles[a]->ID	= 4;
				dom.Particles[a]->IsFree= false;
				dom.Particles[a]->NoSlip= true;
			}
			if (yb>(-0.5*(xb-11.0)) && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 5;
			if (yb>( 0.5*(xb+11.0)) && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 5;
			if (yb>=1.0 && yb<=1.5 && xb>0.0 && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 6;

		}
	}
	dom.DelParticles(5);

    T		= std::min(Tf,Ts);
    std::cout<<"Tf = "<<Tf<<std::endl;
    std::cout<<"Ts = "<<Ts<<std::endl;
    std::cout<<"T  = "<<T<<std::endl;

	dom.Solve(/*tf*/50000.0,/*dt*/T,/*dtOut*/0.2,"test06",2000);
	return 0;
}
MECHSYS_CATCH
