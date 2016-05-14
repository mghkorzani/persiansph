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
        SPH::Domain		dom;

        dom.Dimension		= 2;
        dom.BC.Periodic[0]	= true;
        dom.Nproc		= 12;
    	dom.VisEq		= 3;
    	dom.KernelType		= 4;
    	dom.Scheme		= 0;

        double yb,h,rho;
    	double dx;

    	rho = 998.21;
    	dx = 2.5e-5;
    	h = dx*1.1;
    	dom.InitialDist 	= dx;

        double maz;
        maz=(0.2*h/(0.07));
        std::cout<<maz<<std::endl;
    	dom.GeneralBefore	= & UserAcc;

     	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -0.0006 , 0.0 ), 0.0005*50.0 , 0.00121 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
     		dom.Particles[a]->Cs		= 0.07;
     		dom.Particles[a]->PresEq	= 0;
     		dom.Particles[a]->Mu		= 1.002e-3;
    		dom.Particles[a]->MuRef		= 1.002e-3;
    		dom.Particles[a]->Material	= 1;

    		yb=dom.Particles[a]->x(1);
    		if (yb>=0.0005 || yb<=-0.0005)
    		{
    			dom.Particles[a]->ID		= 2;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    	}

//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/100.0,/*dt*/maz,/*dtOut*/0.05,"test06",1000);
        return 0;
}
MECHSYS_CATCH
