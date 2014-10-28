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
using std::ifstream;

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.Alpha		= 0.05;
	dom.PresEq		= 1;
	dom.KernelType	= 4;
	dom.Nproc		= 64;

	dom.TI			= 0.1;

	double xb,yb,h,rho;
	double dx;

	rho = 998.21;
	dx = 0.002;
	h = dx*1.1;

	dom.Gravity			= 0.0,-9.81,0.0;
	dom.Cs				= 10.0*sqrt(2.0*9.81*1.25);
	dom.InitialDist 	= dx;
	double maz;
	maz=(0.2*h/dom.Cs);


	dom.AddRandomBox(3 ,Vec3_t ( 0 , 0.0 , 0.0 ), 1.0 , 1.25  ,  0 , dx/2.0 ,rho, h);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->Density  = rho*pow((1+rho*9.81*(1.25-dom.Particles[a]->x(1))/(rho*dom.Cs*dom.Cs/7.0)),(1.0/7.0));
		dom.Particles[a]->Densityb = rho*pow((1+rho*9.81*(1.25-dom.Particles[a]->x(1))/(rho*dom.Cs*dom.Cs/7.0)),(1.0/7.0));
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (xb<0.007 || xb>0.993 || yb<0.007)
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}

	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/(100.0*maz),"test06",1500);
	return 0;
}
MECHSYS_CATCH
