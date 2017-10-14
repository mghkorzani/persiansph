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

void UserAcc(SPH::Domain & domi)
{
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
		if (domi.Particles[i]->IsFree)
			domi.Particles[i]->a += Vec3_t(0.002,0.0,0.0);
}

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
//    size_t Nproc;
//    if (argc<=1)
//    	Nproc = 12;
//    else
//    	Nproc = atoi(argv[1]);
    double m0;
    if (argc<=1)
    	m0 = 300.0;
    else
    	m0 = atof(argv[1]);

	SPH::Domain		dom;

	dom.Dimension	= 2;
	dom.Nproc		= 8;
	dom.VisEq		= 1;
	dom.KernelType	= 4;
	dom.Scheme		= 0;
	dom.GeneralBefore	= & UserAcc;
	dom.BC.Periodic[0]	= true;

	double yb,h,rho,Cs,dx,t,Mu,T0,m;

	rho	= 998.21;
	dx 	= 2.5e-5;
	h 	= dx*1.1;
	Cs	= 0.07;
	Mu 	= 1.002e-3;
	m	= m0;
	T0 	= 4.0e-4;
	t	= 0.2*0.125*h*h*rho/(Mu+m*T0);
	dom.InitialDist 	= dx;

	std::cout<<"t = "<<t<<std::endl;
	std::cout<<"m = "<<m<<std::endl;
	std::cout<<"T0 = "<<T0<<std::endl;
	std::cout<<"Mu = "<<Mu<<std::endl;

	dom.AddBoxLength(3 ,Vec3_t ( 0.0 , -0.0006 , 0.0 ), 0.0005 , 0.00121 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->Mu 		= Mu;
		dom.Particles[a]->MuRef 	= Mu;
		dom.Particles[a]->T0 		= T0;
   		dom.Particles[a]->m 		= m;
		dom.Particles[a]->P0		= Cs*Cs*rho*0.01;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Material	= 1;

		yb=dom.Particles[a]->x(1);
		if (yb>=0.0005 || yb<=-0.0005)
		{
			dom.Particles[a]->ID		= 4;
			dom.Particles[a]->IsFree	= false;
			dom.Particles[a]->NoSlip	= true;
		}
	}


//    	dom.WriteXDMF("maz");
	dom.Solve(/*tf*/100.0,/*dt*/t,/*dtOut*/0.005,"test06",1000);
	return 0;
}
MECHSYS_CATCH
