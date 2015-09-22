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

int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Gravity		= 0.002,0.0,0.0;
        dom.Dimension	= 2;
        dom.Cs			= 0.0020;
        dom.P0			= dom.Cs*dom.Cs*998.21*0.01;
        dom.BC.Periodic[0]	= true;
//        dom.MU			= 1.002e-3;
        dom.Nproc		= 8;
    	dom.PresEq		= 0;
    	dom.VisEq		= 1;
    	dom.KernelType	= 4;
    	dom.NoSlip		= true;
    	dom.Shepard		= false;

        double xb,yb,h,rho;
    	double dx;

    	rho = 998.21;
    	dx = 2.5e-5;
    	h = dx*1.1;
    	dom.InitialDist 	= dx;

        double maz;
        maz=(0.00002*h/(dom.Cs));
        std::cout<<maz<<std::endl;

     	dom.AddBoxLength(3 ,Vec3_t ( 0.0 , -0.0006 , 0.0 ), 0.0005 , 0.00121 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->Mu = 1.002e-3;
    		dom.Particles[a]->MuRef = 1.002e-3;
    		dom.Particles[a]->T0 = 4.0e-4;

    		yb=dom.Particles[a]->x(1);
    		if (yb>=0.0005 || yb<=-0.0005)
    		{
    			dom.Particles[a]->ID=4;
    			dom.Particles[a]->IsFree=false;
    		}
    	}


    	for (size_t i=0; i<100; i++)
    	{
    		xb = -dx+(0.0005+dx)/(((0.0005/dx)+2.0)*2.0)*i;
    		dom.AddNSSingleParticle(4,Vec3_t ( xb ,  -0.0005   , 0.0 ),true);
    		dom.AddNSSingleParticle(4,Vec3_t ( xb ,   0.0005   , 0.0 ),true);
    	}


//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/100.0,/*dt*/maz,/*dtOut*/0.005,"test06",1000);
        return 0;
}
MECHSYS_CATCH
