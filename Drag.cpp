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
	dom.ConstVelPeriodic= 0.2;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.Cs			= 4;
	dom.MU			= 0.001;
	dom.P0			= 2000.0;
	dom.PresEq		= 0;


	dom.AddBoxLength(1 ,Vec3_t ( 1.00001125 , 1.00158625+0.0000225 , 0.0 ), 0.0056475 , 0 , 0 , 251 , 1 , 1 , 5.0625e-7 , 1000 , 2.475e-5 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.00001125 , 1.00160875+0.0000225 , 0.0 ), 0.0056475 , 0 , 0 , 251 , 1 , 1 , 5.0625e-7 , 1000 , 2.475e-5 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.00001125 , 1.00021375+0.0000225 , 0.0 ), 0.0056475 , 0 , 0 , 251 , 1 , 1 , 5.0625e-7 , 1000 , 2.475e-5 , true);
	dom.AddBoxLength(1 ,Vec3_t ( 1.00001125 , 1.00019125+0.0000225 , 0.0 ), 0.0056475 , 0 , 0 , 251 , 1 , 1 , 5.0625e-7 , 1000 , 2.475e-5 , true);

	double xa,ya,xb,yb,yc;

	dom.AddRandomBox(3 ,Vec3_t ( 1.000 , 1.000225+0.0000225 , 0.0 ), 0.005625 , 0.00135 ,  0 , 250 , 60 , 1 , 5.0625e-7 , 1000, 2.475e-5);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (((xb-1.001125)*(xb-1.001125)+(yb-1.0009)*(yb-1.0009)<=0.000000018))
			{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
			}
	}
	dom.DelParticles(4);

	xa=1.0010125;
	ya=sqrt (0.00000001265625-(xa-1.001125)*(xa-1.001125))+1.0009;

	while (xa<=1.0012375)
	{
		dom.AddSingleParticle(4,Vec3_t ( xa , ya , 0.0 ), 5.0625e-7 , 1000 , 2.475e-5 , true);

		xb= xa+0.000000000001;
		yb=sqrt (0.00000001265625-(xb-1.001125)*(xb-1.001125))+1.0009;
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.000019)
			{
				xb=xb+0.000000000001;
				yb=sqrt (0.00000001265625-(xb-1.001125)*(xb-1.001125))+1.0009;
			}
		xa=xb;
		ya=yb;
		yc=(-sqrt (0.00000001265625-(xa-1.001125)*(xa-1.001125))+1.0009);

		if (yc>0.0) dom.AddSingleParticle(4,Vec3_t ( xa , yc , 0.0 ), 5.0625e-7 , 1000 , 2.475e-5 , true);
	}

	xa=1.00103505;
	ya=sqrt (8.1e-9-(xa-1.001125)*(xa-1.001125))+1.0009;

	while (xa<=1.001215)
	{
		xb= xa+0.000000000001;
		yb=sqrt (8.1e-9-(xb-1.001125)*(xb-1.001125))+1.0009;
		while(sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))<=0.000019)
			{
				xb=xb+0.000000000001;
				yb=sqrt (8.1e-9-(xb-1.001125)*(xb-1.001125))+1.0009;
			}
		xa=xb;
		ya=yb;

		yc=(-sqrt (8.1e-9-(xa-1.001125)*(xa-1.001125))+1.0009);

		if (yc>0.0)
			{
			dom.AddSingleParticle(4,Vec3_t ( xa , yc , 0.0 ), 5.0625e-7 , 1000 , 2.475e-5 , true);
			dom.AddSingleParticle(4,Vec3_t ( xa , 1.0009+1.0009-yc , 0.0 ), 5.0625e-7 , 1000 , 2.475e-5 , true);
			}
	}
	dom.AddSingleParticle(4,Vec3_t ( 1.001035 , 1.0009 , 0.0 ), 5.0625e-7 , 1000 , 2.475e-5 , true);


//	dom.WriteXDMF("maz");

	dom.Solve(/*tf*/0.5,/*dt*/0.000001,/*dtOut*/0.0002,"test06");
	return 0;
}
MECHSYS_CATCH
