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

        dom.Gravity		= 0.0,-9.81,0.0;
        dom.Dimension	= 2;
        dom.MU			= 1.002e-3;
        dom.Nproc		= 24;
    	dom.PresEq		= 0;
    	dom.VisEq		= 3;
    	dom.KernelType	= 4;
//    	dom.NoSlip		= true;
    	dom.Shepard		= false;


    	dom.BC.InFlow[0]=1;
    	dom.BC.allv= 1.0,0.0,0.0;
    	dom.BC.inv[0] = 1.0,0.0,0.0;
//    	dom.BC.ina[0] = 1.0,0.0,0.0;
    	dom.BC.inDensity[0] = 1000.0;

        double xb,yb,h,rho,Re,H,U;
    	double dx;

    	rho = 998.21;
    	Re = 10.0;
    	H = pow((Re*2.0*dom.MU*dom.MU*3/(2*rho*rho*9.81*.001)),(1.0/3.0));
    	U = rho*9.81*.001*H*H/(2.0*dom.MU);
//    	dom.Cs = 10.0*U;
    	dom.Cs = rho*9.81*1.455e-3;

    	dx = H*2.0/125.0;
    	h = dx*1.1;
    	dom.InitialDist	= dx;
//    	dom.TI = 0.05;

        double maz;
        maz=(0.1*h/(dom.Cs+U));

    	dom.AddBoxLength(3 ,Vec3_t ( 0.0 , -5.0*dx , 0.0 ), 2.0*H , H+5.0*dx  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		yb=dom.Particles[a]->x(1);
    		xb=dom.Particles[a]->x(0);
    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID=4;
    			dom.Particles[a]->IsFree=false;
    		}
    	}

    	xb = -3.0*dx;
    	while (xb<(2.0*H+3.0*dx))
    	{
    		dom.AddNSSingleParticle(4 , Vec3_t ( xb , 0.0 , 0.0 ) , true);
    		xb += dx/10.0;
    	}

//    	dom.WriteXDMF("maz");

    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(100.0*maz),"test06");
        return 0;
}
MECHSYS_CATCH
