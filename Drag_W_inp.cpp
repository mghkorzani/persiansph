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
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    ifstream infile(filename.CStr());
    double velocity;
    double Cs;
    double Mu;
    double P;
    double PressEq;
    double Kernel;
    double h;
    double tim;
    double timestep;
    infile >> velocity;             infile.ignore(200,'\n');
    infile >> Cs;             infile.ignore(200,'\n');
    infile >> Mu;             infile.ignore(200,'\n');
    infile >> P;             infile.ignore(200,'\n');
    infile >> PressEq;             infile.ignore(200,'\n');
    infile >> Kernel;             infile.ignore(200,'\n');
    infile >> h;             infile.ignore(200,'\n');
    infile >> tim;             infile.ignore(200,'\n');
    infile >> timestep;             infile.ignore(200,'\n');

    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.NoSlip		= true;
	dom.PeriodicX	= true;
	dom.PeriodicY	= true;
	dom.ConstVelPeriodic= velocity*4.461319;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.Cs			= 4.461319*velocity*20.0;
	dom.MU			= Mu;

	if (velocity==100.0) dom.P0	= 3900;
	if (velocity==80.0) dom.P0	= 2500;
	dom.PresEq		= abs(PressEq);
	dom.VisEq		= 0;
	dom.KernelType	= abs(Kernel);
	dom.Nproc		= 12;

	dom.TI			= 0.05;
	dom.InitialDist = 0.01;

	double xb,yb;

	dom.AddRandomBox(3 ,Vec3_t ( -125*0.01 , -125*0.01 , 0.0 ), 250.0*0.01 ,250*0.01  ,  0 , 0.005 ,9.9821e-4, 0.012);

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

	double dx,R;
	size_t no;
	dx = 0.01;

	R = 0.12982;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
	}

	R = 0.12116;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , false);
	}

	R = 0.1125;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = 0.10384;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i+M_PI/no);
		yb = R*sin(2*M_PI/no*i+M_PI/no);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = 0.09518;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = 0.08652;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i+M_PI/no);
		yb = R*sin(2*M_PI/no*i+M_PI/no);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	R = 0.07786;
//	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), 4.32238e-8 , 9.9821e-4 , 0.012 , true);
	}

	dom.Solve(/*tf*/tim,/*dt*/(0.2*h/(dom.Cs+dom.ConstVelPeriodic)),/*dtOut*/5.0e-5,"test06");
	return 0;
}
MECHSYS_CATCH
