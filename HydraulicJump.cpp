
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

double Cs,h1,h2,u1,u2,Fr1;

void UserInFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(u1,0.0,0.0);
	Den = 998.21*pow((1+7.0*9.81*(h1-position(1))/(Cs*Cs)),1.0/7.0);
}

void UserAllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	if (position(0)<30.0*h1)
	{
		Vel = Vec3_t(u1,0.0,0.0);
		Den = 998.21*pow((1+7.0*9.81*(h1-position(1))/(Cs*Cs)),1.0/7.0);
	}
	else
	{
		Vel = Vec3_t(u2,0.0,0.0);
		Den = 998.21*pow((1+7.0*9.81*(h2-position(1))/(Cs*Cs)),1.0/7.0);
	}
}
void UserOutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, SPH::Boundary & bdry)
{
	Vel = Vec3_t(u2,0.0,0.0);
	Den = 998.21*pow((1+7.0*9.81*(h2-position(1))/(Cs*Cs)),1.0/7.0);
	if (position(1)>h2) Den=998.21;
}

int main(int argc, char **argv) try
{
	SPH::Domain		dom;

	dom.Gravity		= 0.0,-9.81,0.0;
	dom.Dimension	= 2;
//	dom.MU			= 1.002e-3;
	dom.Nproc		= 24;
	dom.PresEq		= 1;
//	dom.VisEq		= 1;
	dom.Alpha		= 0.05;
	dom.KernelType	= 0;
//	dom.Shepard		= true;
//	dom.TI			= 0.025;
//	dom.XSPH		= 0.5;
	dom.BCDensityUpdate = false;
	dom.BC.MassConservation = true;

	h1		= 1.00;
	Fr1		= 5.0;
	h2		= h1/2.0*(-1+sqrt(1+8*Fr1*Fr1));
	u1		= Fr1*sqrt(9.81*h1);
	u2		= u1*h1/h2;
	std::cout<<"h1 = "<<h1<<std::endl;
	std::cout<<"Fr1 = "<<Fr1<<std::endl;
	std::cout<<"h2 = "<<h2<<std::endl;
	std::cout<<"u1 = "<<u1<<std::endl;
	std::cout<<"u2 = "<<u2<<std::endl;

	dom.BC.InOutFlow	= 3;
	dom.BC.allv			= u2,0.0,0.0;
	dom.BC.inv			= u1,0.0,0.0;
	dom.BC.outv			= u2,0.0,0.0;
	dom.BC.inDensity	= 998.21;
	dom.BC.outDensity	= 998.21;
	dom.BC.allDensity	= 998.21;

	dom.InCon	= & UserInFlowCon;
	dom.AllCon	= & UserAllFlowCon;
	dom.OutCon	= & UserOutFlowCon;

	double dx, xb, yb, h, rho;

	rho		= 998.21;
	dom.Cs	= 10.0 * u1;
	dx		= h1 / 20.0;
	h		= dx * 1.1;
	Cs		= dom.Cs;

	dom.InitialDist	= dx;
	dom.DomMax(1)	= 2 * h2;
	double maz		= (0.1*h/(dom.Cs+u1));

	dom.AddBoxLength(1 ,Vec3_t ( 0.0 , -3.0*dx , 0.0 ), 70.0*h1 , h2+3.0*dx  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb = dom.Particles[a]->x(0);
		yb = dom.Particles[a]->x(1);
		if (yb<0.0)
		{
			dom.Particles[a]->ID		= 2;
			dom.Particles[a]->Density	= 998.21*pow((1+7.0*9.81*(h1-yb)/(Cs*Cs)),1.0/7.0)*1.005;
			dom.Particles[a]->Densityb	= 998.21*pow((1+7.0*9.81*(h1-yb)/(Cs*Cs)),1.0/7.0)*1.005;
			dom.Particles[a]->IsFree	= false;
		}
		if (yb>h1 && xb<30.0*h1 && dom.Particles[a]->ID == 1)
		{
			dom.Particles[a]->ID		= 3;
		}
	}
	dom.DelParticles(3);

	dom.Solve(/*tf*/50000.0,/*dt*/maz,/*dtOut*/0.1,"test06",2000);
	return 0;
}
MECHSYS_CATCH

