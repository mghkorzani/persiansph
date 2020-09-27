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

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension		= 2;
	dom.BC.Periodic[1]	= true;
//	dom.BC.Periodic[0]	= true;
	dom.VisEq		= 3;
	dom.KernelType		= 4;
	dom.Nproc		= 24;
	dom.Scheme		= 0;

	double xb,yb,h,rho,mass,U,Cs,P0,Mu,dx,R,Rc,Re,T;
	size_t no;

	Mu		= 1.002e-3;
	rho		= 998.21;
	dx		= 0.003;
	h		= dx*1.1;
	Rc		= 0.05;
	mass		= dx*dx*rho;
	Re		= 100.0;
	U		= Re*Mu/(rho*2.0*Rc);
	Cs		= U*10.0;
	P0		= Cs*Cs*rho*0.08;
	T		= (0.25*h/(Cs+U));
	dom.DtAcc	= 0.1;

	dom.BC.InOutFlow	= 3;
	dom.BC.allv		= U,0.0,0.0;
	dom.BC.inDensity	= rho;
	dom.BC.inv		= U,0.0,0.0;
//	dom.BC.outv		= U,0.0,0.0;
//	dom.BC.outDensity	= rho;
	dom.InitialDist 	= dx;
//	dom.BC.MassConservation = true;

	std::cout<<"Re = "<<Re<<std::endl;
	std::cout<<"V  = "<<U<<std::endl;
	std::cout<<"Cs = "<<Cs<<std::endl;
	std::cout<<"P0 = "<<P0<<std::endl;
	std::cout<<"Time Step = "<<T<<std::endl;
	std::cout<<"Resolution = "<<(2.0*Rc/dx)<<std::endl;

	dom.AddBoxLength(1 ,Vec3_t ( -8.0*Rc , -8.0*Rc , 0.0 ), 16.0*Rc , 16.0*Rc  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if ((xb*xb+yb*yb)<((Rc+h/2.0)*(Rc+h/2.0)))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);


	for (size_t j=0;j<6;j++)
	{
		R = Rc-dx*j;
		no = ceil(2*M_PI*R/dx);
		for (size_t i=0; i<no; i++)
		{
			xb = R*cos(2*M_PI/no*i);
			yb = R*sin(2*M_PI/no*i);
			dom.AddSingleParticle(2,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
			dom.Particles[dom.Particles.Size()-1]->NoSlip = true;
		}
	}

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->Shepard	= true;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->P0		= P0;
		dom.Particles[a]->Mu		= Mu;
		dom.Particles[a]->MuRef		= Mu;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Material	= 1;
	}


	dom.Solve(/*tf*/20000.0,/*dt*/T,/*dtOut*/(200.0*T),"test06",100);
	return 0;
}
MECHSYS_CATCH
