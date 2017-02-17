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

void NewUserOutput(SPH::Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = Particles->ShearRate;
}


int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Dimension		= 2;
        dom.BC.Periodic[0]	= true;
        dom.Nproc		= 8;
    	dom.VisEq		= 3;
    	dom.KernelType		= 4;
    	dom.Scheme		= 0;

        double yb,h,rho;
    	double dx;

    	rho	= 998.21;
    	dx	= 2.5e-5;
    	h	= dx*1.1;
    	dom.InitialDist 	= dx;

        double maz;
        maz	=(0.2*h/(0.08));
        std::cout<<maz<<std::endl;

    	dom.AddBoxLength(3 ,Vec3_t ( 0.0 , -0.0000875001 , 0.0 ), 0.0005 , 0.00118  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
		dom.Particles[a]->LES		= true; //Just used to activate ShearRate calculation
		dom.Particles[a]->CSmag		= 0.0;  //To deactive LES for the above purpose (No LES used for this simulation)

    		dom.Particles[a]->Cs		= 0.08;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Mu		= 1.002e-3;
    		dom.Particles[a]->MuRef		= 1.002e-3;
    		dom.Particles[a]->Material	= 1;

    		yb=dom.Particles[a]->x(1);
    		if (yb>=0.0009901)
    		{
    			dom.Particles[a]->ID		= 4;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->v			= 2.5e-5,0.0,0.0;
    			dom.Particles[a]->vb		= 2.5e-5,0.0,0.0;
    			dom.Particles[a]->NoSlip	= true;
   		}
    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID		= 5;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->NoSlip	= true;
    		}
    	}

	dom.OutputName[0]	= "ShearRate";
        dom.UserOutput		= & NewUserOutput;


//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(1000.0*maz),"test06",250);
        return 0;
}
MECHSYS_CATCH
