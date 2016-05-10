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

#include <Domain.h>


using std::cout;
using std::endl;
using std::ifstream;

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 3;

	dom.MU			= 1.002e-3;
	dom.PresEq		= 0;
	dom.VisEq		= 0;
	dom.KernelType	= 0;
	dom.Nproc		= 20;

	double xb,yb,zb,h,rho,mass;
	double dx;

	rho = 998.21;
	dx = 0.016;
	h = dx*1.1;
	mass = dx*dx*dx*rho;

	dom.Cs				= 10.0 * sqrt(2.0*9.81*0.5);
	dom.InitialDist 	= dx;
	double maz = (0.2*h/(dom.Cs + dom.Cs/10.0));

	dom.AddBoxLength(10,Vec3_t ( -7.7-2.0*dx , -1.8-2.0*dx , -2.0*dx ), 18+4.1*dx , 3.6+4.1*dx , 2.1*dx     , dx/2.0 , rho , h , 1 , 0 , false , true );
	dom.AddBoxLength(8 ,Vec3_t ( -7.7-2.0*dx , -1.8-2.0*dx , 0.0     ), 2.1*dx    , 3.6+4.1*dx , 0.5+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , true );
	dom.AddBoxLength(7 ,Vec3_t ( 10.3        , -1.8-2.0*dx , 0.0     ), 2.1*dx    , 3.6+4.1*dx , 0.5+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , true );
	dom.AddBoxLength(6 ,Vec3_t ( -7.7        , -1.8-2.0*dx , 0.0     ), 18+0.1*dx , 2.1*dx     , 0.5+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , true );
	dom.AddBoxLength(5 ,Vec3_t ( -7.7        ,  1.8        , 0.0     ), 18+0.1*dx , 2.1*dx     , 0.5+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , true );
//  If dx is changed, -0.804 should be revised.
	dom.AddBoxLength(2 ,Vec3_t ( -0.804      , -1.8        , 0.0     ), 0.8+0.1*dx, 3.6+0.1*dx , 0.5+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , true );
	dom.AddBoxLength(1 ,Vec3_t ( -7.7        , -1.8        , 0.0     ), 6.9+0.1*dx, 3.6+0.1*dx , 0.4+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , false);
//  If dx is changed, -0.004 should be revised.
	dom.AddBoxLength(1 ,Vec3_t ( -0.004      , -1.8        , 0.0     ), 10.3+0.004+0.1*dx, 3.6+0.1*dx , 0.02+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , false);
	dom.AddBoxLength(3 ,Vec3_t (  3.8        , -0.23       , 0.0     ), 0.8+0.1*dx, 0.4+0.1*dx , 0.4+0.1*dx , dx/2.0 , rho , h , 1 , 0 , false , true);

	double teta,inc,m;
	size_t no;
	m    = 0.155/0.34;
	teta = atan(m);
	inc  = dx * cos (teta);
	no   = ceil(0.34/inc);
	inc  = 0.34/no;

	for (size_t j=0; j<no-1; j++)
		for (size_t i=0; i<(18/dx-dx); i++)
		{
			dom.AddSingleParticle(4,Vec3_t ( -7.7+dx/2.0+i*dx ,  1.46+inc*(j+2)    , (m*((1.46+inc*(j+1))-1.46)) )  , mass , rho , h , true);
			dom.AddSingleParticle(4,Vec3_t ( -7.7+dx/2.0+i*dx ,  -(1.46+inc*(j+2)) , (-m*(-(1.46+inc*(j+1))+1.46)) ), mass , rho , h , true);
		}
	for (size_t j=0; j<no-3; j++)
		for (size_t i=0; i<(18/dx-dx); i++)
		{
			dom.AddSingleParticle(4,Vec3_t ( -7.7+dx/2.0+i*dx ,  1.46+inc*(j+4)    , (m*((1.46+inc*(j+1))-1.46)) )  , mass , rho , h , true);
			dom.AddSingleParticle(4,Vec3_t ( -7.7+dx/2.0+i*dx ,  -(1.46+inc*(j+4)) , (-m*(-(1.46+inc*(j+1))+1.46)) ), mass , rho , h , true);
		}

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		zb=dom.Particles[a]->x(2);

		if (xb<0 && xb>-0.8 && (yb>-0.5 && yb<0.5) && dom.Particles[a]->ID==2)
		{
			dom.Particles[a]->ID=1;
			dom.Particles[a]->IsFree=true;
		}

		if (zb>0.4 && dom.Particles[a]->ID==1)
		{
			dom.Particles[a]->ID=11;
		}

		if (yb<-1.46  && zb<(-m*(yb+1.46)) && dom.Particles[a]->ID==1 )
		{
			dom.Particles[a]->ID=12;
		}

		if (yb>1.46  && zb<(m*(yb-1.46)) && dom.Particles[a]->ID==1 )
		{
			dom.Particles[a]->ID=12;
		}

		if (dom.Particles[a]->ID==3)
		{
			dom.Particles[a]->x(0)=3.8+(xb-3.8)*cos(64.0/180.0*M_PI)-(yb+0.23)*sin(64.0/180.0*M_PI);
			dom.Particles[a]->x(1)=-0.23+(xb-3.8)*sin(64.0/180.0*M_PI)+(yb+0.23)*cos(64.0/180.0*M_PI);
		}

		if (xb>((yb-0.67)/2.057+3.79) && xb<((yb-0.49)/2.057+4.15) && yb>((xb-3.44)*-0.5-0.05)  && yb<((xb-4.15)*-0.5+0.49) && dom.Particles[a]->ID==1 )
		{
			dom.Particles[a]->ID=13;
		}

		if (xb<0 && xb>-0.8 && dom.Particles[a]->ID==4)
		{
			dom.Particles[a]->ID=14;
			dom.Particles[a]->IsFree=true;
		}
	}

	dom.DelParticles(11);
	dom.DelParticles(12);
	dom.DelParticles(13);
	dom.DelParticles(14);

	dom.WriteXDMF("maz");
	std::cout<<"No of Particles =  "<<dom.Particles.Size()<<std::endl;
//	dom.Solve(/*tf*/2500000.0,/*dt*/maz,/*dtOut*/(500.0*maz),"test06",80);
	return 0;
}
MECHSYS_CATCH
