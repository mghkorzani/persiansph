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
	dom.ConstVelPeriodic= 44.6;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.Cs			= 800.0;
	dom.MU			= 1.002e-3;
	dom.P0			= 50.0;
	dom.PresEq		= 0;
	dom.KernelType	= 2;
	dom.Nproc		= 2;

	dom.TI			= 0.05;
	dom.InitialDist = 0.01;

	double xa,ya,xb,yb,yc;

	dom.AddRandomBox(3 ,Vec3_t ( -56.25*0.01 , -56.1*0.01 , 0.0 ), 135.0*0.01 ,112.5*0.01  ,  0 , 0.005 ,9.9821e-4, 0.012);

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
		dom.AddSingleParticle(3,Vec3_t ( xa ,  ya , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
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
			dom.AddSingleParticle(3,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
		}
	}

	xa=-0.12;
	ya=sqrt (0.0144-xa*xa);

	while (xa<=0.12)
	{
		dom.AddSingleParticle(3,Vec3_t ( xa ,  ya , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
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
			dom.AddSingleParticle(3,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
		}
	}
	dom.AddSingleParticle(3,Vec3_t ( 0.12 , 0.0 , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);

	xa=-0.1125;
	ya=sqrt (0.01265625-xa*xa);

	while (xa<=0.1125)
	{
		dom.AddSingleParticle(4,Vec3_t ( xa ,  ya , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
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
			dom.AddSingleParticle(4,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
		}
	}

	xa=-0.1051;
	ya=sqrt (0.011025-xa*xa);

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
			dom.AddSingleParticle(4,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa ,  yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
		}
	}
	dom.AddSingleParticle(4,Vec3_t ( 0.1050 , 0.0 , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);

	xa=-0.0975;
	ya=sqrt (0.00950625-xa*xa);

	while (xa<=0.0975)
	{
		xb= xa+0.0001;
		yb=sqrt (0.00950625-xb*xb);
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.0001;
			yb=sqrt (0.00950625-xb*xb);
		}
		xa=xb;
		ya=yb;

		yc=sqrt (0.00950625-xa*xa);

		if (yc>0.0)
		{
			dom.AddSingleParticle(4,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa ,  yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
		}
	}
	dom.AddSingleParticle(4,Vec3_t ( -0.0975 , 0.0 , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);

	xa=-0.0925;
	ya=sqrt (0.00855625-xa*xa);
	while (xa<=0.0925)
	{
		xb= xa+0.0001;
		yb=sqrt (0.00855625-xb*xb);
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.009)
		{
			xb=xb+0.0001;
			yb=sqrt (0.00855625-xb*xb);
		}
		xa=xb;
		ya=yb;

		yc=sqrt (0.00855625-xa*xa);

		if (yc>0.0)
		{
			dom.AddSingleParticle(4,Vec3_t ( xa , -yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa ,  yc , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
		}
	}
	dom.AddSingleParticle(4,Vec3_t ( -0.0925 , 0.0 , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);


//	dom.WriteXDMF("maz");

	dom.Solve(/*tf*/0.2,/*dt*/3.5e-6,/*dtOut*/5.0e-5,"test06");
	return 0;
}
MECHSYS_CATCH
