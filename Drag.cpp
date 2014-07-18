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
	dom.Dimension	= 2;

	dom.Periodic	= true;
	dom.ConstVelPeriodic= 10;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.Cs			= 500;
	dom.Alpha		= 0.03;
	dom.P0			= 5.0;
	dom.PresEq		= 0;


	dom.TI			= 0.05;
	dom.InitialDist = 0.01;

	dom.AddBoxLength(1 ,Vec3_t ( 1.0-113*0.01 , 1.0-36.5*0.01 , 0.0 ), 226*0.01 , 0 , 0 , 226 , 1 , 1 , 1.0e-7 , 1.0e-3 , 0.011 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.0-113*0.01 , 1.0-35.5*0.01 , 0.0 ), 226*0.01 , 0 , 0 , 226 , 1 , 1 , 1.0e-7 , 1.0e-3 , 0.011 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.0-113*0.01 , 1.0+35*0.01   , 0.0 ), 226*0.01 , 0 , 0 , 226 , 1 , 1 , 1.0e-7 , 1.0e-3 , 0.011 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.0-113*0.01 , 1.0+36*0.01   , 0.0 ), 226*0.01 , 0 , 0 , 226 , 1 , 1 , 1.0e-7 , 1.0e-3 , 0.011 , true);

	double xa,ya,xb,yb,yc;

	dom.AddRandomBox(3 ,Vec3_t ( 1.0-113.5*0.01 , 1.0-35*0.01 , 0.0 ), 226*0.01 ,70*0.01  ,  0 , 0.005 ,1.0e-3, 0.011);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (((xb-1.0)*(xb-1.0)+(yb-1.0)*(yb-1.0)<=0.014))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

	xa=0.8875;
	ya=sqrt (0.01265625-(xa-1.0)*(xa-1.0))+1.0;

	while (xa<=1.1125)
	{
		xb= xa+0.0001;
		yb=sqrt (0.01265625-(xb-1.0)*(xb-1.0))+1.0;
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.0001;
			yb=sqrt (0.01265625-(xb-1.0)*(xb-1.0))+1.0;
		}
		xa=xb;
		ya=yb;
		yc=(-sqrt (0.01265625-(xa-1.0)*(xa-1.0))+1.0);

		if (yc>0.0)
		{
			dom.AddSingleParticle(4,Vec3_t ( xa ,         yc , 0.0 ), 1.0e-7 , 1.0e-3 , 0.011 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa , 1.0+1.0-yc , 0.0 ), 1.0e-7 , 1.0e-3 , 0.011 , true);
		}
	}

	xa=0.8975;
	ya=sqrt (0.01050625-(xa-1.0)*(xa-1.0))+1.0;

	while (xa<=1.1025)
	{
		xb= xa+0.0001;
		yb=sqrt (0.01050625-(xb-1.0)*(xb-1.0))+1.0;
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.0001;
			yb=sqrt (0.01050625-(xb-1.0)*(xb-1.0))+1.0;
		}
		xa=xb;
		ya=yb;

		yc=(-sqrt (0.01050625-(xa-1.0)*(xa-1.0))+1.0);

		if (yc>0.0)
		{
			dom.AddSingleParticle(4,Vec3_t ( xa ,         yc , 0.0 ), 1.0e-7 , 1.0e-3 , 0.011 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa , 1.0+1.0-yc , 0.0 ), 1.0e-7 , 1.0e-3 , 0.011 , true);
		}
	}
	dom.AddSingleParticle(4,Vec3_t ( 1.1025 , 1.0 , 0.0 ), 1.0e-7 , 1.0e-3 , 0.011 , true);

//	dom.WriteXDMF("maz");

	dom.Solve(/*tf*/0.3,/*dt*/0.000001,/*dtOut*/0.00001,"test06");
	return 0;
}
MECHSYS_CATCH
