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

        dom.Gravity		= 0.0,0.0,0.002;
        dom.Dimension	= 3;
        dom.Cs			= 0.0025;
        dom.PeriodicZ	= true;
        dom.MU			= 1.002e-3;
        dom.Nproc		= 21;
    	dom.PresEq		= 0;
    	dom.VisEq		= 1;
    	dom.KernelType	= 0;
//    	dom.NoSlip		= true;

        double xb,yb,h,rho;
    	double dx;

    	rho = 998.21;
    	dx = 2.5e-5;
    	h = dx*1.1;
    	dom.InitialDist 	= dx;

        double maz;
        maz=(0.001*h/(dom.Cs+0.00025));

    	dom.AddBoxLength(3 ,Vec3_t ( -0.0006 , -0.0006 , 0.0 ), 0.0012 , 0.0012  ,  0.0005 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		yb=dom.Particles[a]->x(1);
    		xb=dom.Particles[a]->x(0);
    		if ((yb*yb+xb*xb)>=(0.0005*0.0005))
    		{
    			dom.Particles[a]->ID=4;
    			dom.Particles[a]->IsFree=false;
    		}

    		if ((yb*yb+xb*xb)>(0.0006*0.0006))
    		{
    			dom.Particles[a]->ID=5;
    			dom.Particles[a]->IsFree=false;
    		}
    	}
    	dom.DelParticles(5);


//    	for (size_t i=0; i<50; i++)
//    	{
//    		xb = -0.0001+0.0007/50.0*i;
//        	for (size_t j=0; j<50; j++)
//        	{
//				yb = 0.0005*cos(j*2.0*M_PI/50.0);
//				zb = 0.0005*sin(j*2.0*M_PI/50.0);
//				dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , zb ),true);
//				dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , zb ),true);
//        	}
//    	}


    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(100.0*maz),"test06",1500);
        return 0;
}
MECHSYS_CATCH
