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

	dom.PeriodicX	= true;
	dom.PeriodicY	= true;
	dom.ConstVelPeriodic= 85.0;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.Cs			= 1500.0;
//	dom.Alpha		= 0.05;
	dom.MU			= 1.002e-3;
	dom.P0			= 100.0;
	dom.PresEq		= 0;
//	dom.Shepard		= false;


	dom.TI			= 0.05;
	dom.InitialDist = 0.01;

	double xa,ya,xb,yb,yc;

	dom.AddRandomBox(3 ,Vec3_t ( -80*0.01 , -41.2157*0.01 , 0.0 ), 225*0.01 ,83*0.01  ,  0 , 0.005 ,9.9821e-4, 0.011);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if ((xb*xb+yb*yb<=0.018))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

	xa=-0.1275;
	ya=sqrt (0.01625625-xa*xa);

	while (xa<=0.1275)
	{
		dom.AddSingleParticle(3,Vec3_t ( xa ,  ya , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , false);
		xb= xa+0.00001;
		yb=sqrt (0.01625625-xb*xb);
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.00001;
			yb=sqrt (0.01625625-xb*xb);
		}
		xa=xb;
		ya=yb;
		yc=sqrt (0.01625625-xa*xa);

		if (yc>0.0)
		{
			dom.AddSingleParticle(3,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , false);
		}
	}

	xa=-0.12;
	ya=sqrt (0.0144-xa*xa);

	while (xa<=0.12)
	{
		dom.AddSingleParticle(3,Vec3_t ( xa ,  ya , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , false);
		xb= xa+0.00001;
		yb=sqrt (0.0144-xb*xb);
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.00001;
			yb=sqrt (0.0144-xb*xb);
		}
		xa=xb;
		ya=yb;
		yc=sqrt (0.0144-xa*xa);

		if (yc>0.0)
		{
			dom.AddSingleParticle(3,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , false);
		}
	}
	dom.AddSingleParticle(3,Vec3_t ( 0.12 , 0.0 , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , false);

	xa=-0.1125;
	ya=sqrt (0.01265625-xa*xa);

	while (xa<=0.1125)
	{
		dom.AddSingleParticle(4,Vec3_t ( xa ,  ya , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , true);
		xb= xa+0.00001;
		yb=sqrt (0.01265625-xb*xb);
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.00001;
			yb=sqrt (0.01265625-xb*xb);
		}
		xa=xb;
		ya=yb;
		yc=sqrt (0.01265625-xa*xa);

		if (yc>0.0)
		{
			dom.AddSingleParticle(4,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , true);
		}
	}

	xa=-0.1050;
	ya=sqrt (0.015025-xa*xa);

	while (xa<=0.1050)
	{
		xb= xa+0.0001;
		yb=sqrt (0.011025-xb*xb);
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.0001;
			yb=sqrt (0.011025-xb*xb);
		}
		xa=xb;
		ya=yb;

		yc=sqrt (0.011025-xa*xa);

		if (yc>0.0)
		{
			dom.AddSingleParticle(4,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa ,  yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.011 , true);
		}
	}

//	dom.WriteXDMF("maz");

	dom.Solve(/*tf*/0.2,/*dt*/1.5e-6,/*dtOut*/5.0e-5,"test06");
	return 0;
}
MECHSYS_CATCH
