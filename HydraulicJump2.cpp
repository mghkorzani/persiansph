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

double Cs,Mu,Density,h1;

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(Density*9.81*0.0547/(2.0*Mu)*(2*0.001*position(1)-position(1)*position(1)),0.0,0.0);
	Den = Density*(1+9.81*(0.001-position(1))/(Cs*Cs));
}
void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(Density*9.81*0.006/(2.0*Mu)*(2*0.00209*position(1)-position(1)*position(1)),0.0,0.0);
	Den = Density*(1+9.81*(0.00209-position(1))/(Cs*Cs));
}
void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(0)<(5.0*h1))
	{
		Vel = Vec3_t(Density*9.81*0.0547/(2.0*Mu)*(2*0.001*position(1)-position(1)*position(1)),0.0,0.0);
		Den = Density*(1+9.81*(0.001-position(1))/(Cs*Cs));
	}
	else
	{
		Vel = Vec3_t(Density*9.81*0.006/(2.0*Mu)*(2*0.00209*position(1)-position(1)*position(1)),0.0,0.0);
		Den = Density*(1+9.81*(0.00209-position(1))/(Cs*Cs));
	}
}

int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Gravity		= 0.0,-9.81,0.0;
        dom.Dimension	= 2;
        dom.MU			= 1.002e-3;
        dom.Nproc		= 12;
    	dom.PresEq		= 0;
    	dom.VisEq		= 1;
    	dom.KernelType	= 2;
    	dom.Shepard		= true;


    	dom.BC.InOutFlow =1;
    	dom.BC.allv = 0.148,0.0,0.0;
    	dom.BC.inv = 0.1783,0.0,0.0;
    	dom.BC.inDensity = 998.21;
//    	dom.BC.outDensity = 998.21;
    	dom.BC.allDensity = 998.21;
//    	dom.BC.outv = 0.0852,0.0,0.0;

        double xb,yb,h,rho,Re,H,U;
    	double dx;

    	rho = 998.21;
    	H = 0.001;
    	U = 0.2673;
//    	dom.Cs = rho*9.81*H*2.0;
    	dom.Cs = 2.7;
//    	dom.TI = 0.05;

    	dx = H/70.0;
    	h = dx*1.1;
    	dom.InitialDist	= dx;

    	Cs = dom.Cs;
    	Mu = dom.MU;
    	Density = rho;
    	h1 = H;
    	dom.InCon = & UserInFlowCon;
//    	dom.OutCon = & UserOutFlowCon;
    	dom.AllCon = & UserAllFlowCon;

        double maz;
        maz=(0.10*h/(dom.Cs+U));

    	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -5.0*dx , 0.0 ), 15.0*H , 2.09*H+5.0*dx  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
    	dom.DomMax(1) = 2.5*H;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		xb=dom.Particles[a]->x(0);
    		yb=dom.Particles[a]->x(1);
    		if (yb<0.0)
    		{
    			dom.Particles[a]->ID=2;
    			dom.Particles[a]->IsFree=false;
    		}
    		if (yb<(1.5*H) && xb >(15.0*H-5.0*dx))
    		{
    			dom.Particles[a]->ID=2;
    			dom.Particles[a]->IsFree=false;
    		}
    		if (yb>0.001 && xb<(5.0*H))
    		{
    			dom.Particles[a]->ID=3;
    		}
   	}
    	dom.DelParticles(3);

    	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(0.0005),"test06");
        return 0;
}
MECHSYS_CATCH
