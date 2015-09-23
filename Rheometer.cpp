/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo Torres  *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#include "./Source/Domain.h"


using std::cout;
using std::endl;

void UserMazi(SPH::Domain & domi)
{
	for (size_t a=0; a<domi.Particles.Size(); a++)
	{
		domi.Particles[a]->x(1);
		if (domi.Particles[a]->x(1)<0.0005)
		{
			domi.Particles[a]->IsFree=false;
		}
	}
}


int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Dimension	= 2;
        dom.Cs			= 0.0025;
        dom.BC.Periodic[0]	= true;
//        dom.MU			= 1.002e-3;
        dom.Nproc		= 8;
    	dom.PresEq		= 0;
    	dom.VisEq		= 3;
    	dom.KernelType	= 4;
    	dom.NoSlip		= true;

    	dom.RigidBody	= true;
    	dom.RBTag		= 4;
    	dom.vel			= 2.5e-5,0.0,0.0;
    	dom.vellim		= 2.5e-5,0.0,0.0;

        double xb,yb,h,rho;
    	double dx;

    	rho = 998.21;
    	dx = 2.5e-5;
    	h = dx*1.1;
    	dom.InitialDist 	= dx;

        double maz;
        maz=(0.005*h/(dom.Cs+0.00025));
        std::cout<<maz<<std::endl;
//    	dom.GeneralBefore	= & UserMazi;
//    	dom.GeneralAfter	= & UserMazi;

    	dom.AddBoxLength(3 ,Vec3_t ( 0.0 , -0.0000875001 , 0.0 ), 0.0005 , 0.00118  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->Mu = 1.002e-3;
    		dom.Particles[a]->MuRef = 1.002e-3;

    		yb=dom.Particles[a]->x(1);
    		if (yb>=0.0009901)
    		{
    			dom.Particles[a]->ID=4;
    			dom.Particles[a]->IsFree=false;
    		}
    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID=5;
    			dom.Particles[a]->IsFree=false;
    		}
    	}


    	for (size_t i=0; i<100; i++)
    	{
    		xb = -0.0001+0.0007/100.0*i;
    		dom.AddNSSingleParticle(5,Vec3_t ( xb ,  0.0   , 0.0 ),true);
    		dom.AddNSSingleParticle(4,Vec3_t ( xb ,  0.001 , 0.0 ),true);
    	}


//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(1000.0*maz),"test06",1500);
        return 0;
}
MECHSYS_CATCH
