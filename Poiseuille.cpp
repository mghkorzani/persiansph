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

#include <Domain.h>


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Gravity		= 0.002,0.0,0.0;
        dom.Dimension	= 2;
        dom.Cs			= 0.0025;
        dom.PeriodicX	= true;
        dom.MU			= 1.002e-3;
        dom.Nproc		= 24;
    	dom.PresEq		= 0;
    	dom.VisEq		= 3;
    	dom.KernelType	= 4;
    	dom.NoSlip		= true;

        double xb,yb,h,rho;
    	double dx;

    	rho = 998.21;
    	dx = 2.5e-5;
    	h = dx*1.1;
    	dom.InitialDist 	= dx;

        double maz;
        maz=(0.005*h/(dom.Cs+0.00025));
        std::cout<<maz<<std::endl;

    	dom.AddRandomBox(3 ,Vec3_t ( 0.0 , -0.0000875001 , 0.0 ), 0.0005 , 0.00118  ,  0 , dx/2.0 ,rho, h, 1 , 0 );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		yb=dom.Particles[a]->x(1);
    		if (yb>=0.0009901 || yb<0.0)
    		{
    			dom.Particles[a]->ID=4;
    			dom.Particles[a]->IsFree=false;
    		}
    	}


    	for (size_t i=0; i<50; i++)
    	{
    		xb = -0.0001+0.0007/50.0*i;
    		dom.AddSingleParticle(Vec3_t ( xb ,  0.0 , 0.0 ));
    		dom.AddSingleParticle(Vec3_t ( xb ,  0.001 , 0.0 ));
    	}


//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(100.0*maz),"test06",1500);
        return 0;
}
MECHSYS_CATCH
