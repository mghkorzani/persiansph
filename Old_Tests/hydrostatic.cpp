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
#include "./Source/Interaction.h"


using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.PresEq		= 0;
	dom.KernelType	= 0;
	dom.Nproc		= 8;
	dom.Alpha		= 1.0;
	dom.Beta		= 1.0;
	dom.NoSlip		= true;
	dom.Shepard		= true;
	dom.TI			= 0.3;
	dom.XSPH		= 0.5;

	double xb,yb,h,rho;
	double dx,H,L;

	H				= 1.2;
	L				= 0.8;
	rho				= 998.21;
	dx				= 0.02;
	h				= dx*1.5;

	dom.Gravity		= 0.0,-9.81,0.0;
	dom.Cs			= 10.0*sqrt(2.0*9.81*1.5);
	dom.InitialDist	= dx;
	double maz;
	maz=(0.1*h/dom.Cs);

	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , -3.0*dx , 0.0 ), L , H + 3.0*dx + dx/10.0  ,  0 , dx/2.0 ,rho, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->Material = 1;
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (xb<(-0.4+3.0*dx))
		{
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
		}
		if (xb>(0.4-3.0*dx))
		{
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
		}
		if (yb<0.0)
		{
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
		}
//		if (dom.Particles[a]->ID==1)
//		{
//			dom.Particles[a]->Density  = rho*pow((1+rho*9.81*(H-dom.Particles[a]->x(1))/(rho*dom.Cs*dom.Cs/7.0)),(1.0/7.0));
//			dom.Particles[a]->Densityb = rho*pow((1+rho*9.81*(H-dom.Particles[a]->x(1))/(rho*dom.Cs*dom.Cs/7.0)),(1.0/7.0));
			dom.Particles[a]->Density  = rho*((1+9.81*(H-dom.Particles[a]->x(1))/(dom.Cs*dom.Cs)));
			dom.Particles[a]->Densityb = rho*((1+9.81*(H-dom.Particles[a]->x(1))/(dom.Cs*dom.Cs)));

//		}
	}

	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(400.0*maz),"test06",10000);
	return 0;
}
MECHSYS_CATCH
