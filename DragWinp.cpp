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
    double Ref;
    double Csf;
    double P0f;
    infile >> Ref;		infile.ignore(200,'\n');
    infile >> Csf;		infile.ignore(200,'\n');
    infile >> P0f;		infile.ignore(200,'\n');

    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.NoSlip		= true;
	dom.PeriodicX	= true;
	dom.PeriodicY	= true;

	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.MU			= 1.002e-3;
	dom.PresEq		= 0;
	dom.VisEq		= 3;
	dom.KernelType	= 4;
	dom.Nproc		= 24;

//	dom.TI			= 0.05;

	double xb,yb,h,rho,mass;
	double dx,R,Rc,Re;
	size_t no;

	rho = 998.21;
	dx = 0.002;
	h = dx*1.1;
	Rc = 0.05;
	mass = (sqrt(3.0)*dx*dx/4.0)*rho;
	Re = Ref;

//	dom.ConstVelPeriodic= Re*dom.MU/(rho*2.0*Rc);
	dom.vel				= Re*dom.MU/(rho*2.0*Rc),0.0,0.0;
	double temp 		= norm(dom.vel);
	dom.Acc				= 6e-6,0.0,0.0;
	dom.Cs				= temp*Csf;
	dom.P0				= dom.Cs*dom.Cs*rho*P0f;
	dom.InitialDist 	= dx;
	double maz;
	maz=(0.2*h/(dom.Cs+temp));

	std::cout<<"Re = "<<Re<<std::endl;
	std::cout<<"V  = "<<temp<<std::endl;
	std::cout<<"Cs = "<<dom.Cs<<std::endl;
	std::cout<<"P0 = "<<dom.P0<<std::endl;
	std::cout<<"Time Step = "<<maz<<std::endl;
	std::cout<<"Resolution = "<<(2.0*Rc/dx)<<std::endl;

	dom.AddRandomBox(3 ,Vec3_t ( -5.0*Rc , -10.0*Rc , 0.0 ), 30.0*Rc , 20.0*Rc  ,  0 , dx/2.0 ,rho, h);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if ((xb*xb+yb*yb)<((Rc+3.0/2.0*sqrt(3.0)*dx)*(Rc+3.0/2.0*sqrt(3.0)*dx)))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

	R = Rc+sqrt(3.0)*dx;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , false);
	}

	R = Rc+sqrt(3.0)/2.0*dx;
	no = ceil(2*M_PI*R/dx);
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(3,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , false);
	}

	//No-Slip BC
	R = Rc;
	no = ceil(2*M_PI*R/(dx/5.0));
	for (size_t i=0; i<no; i++)
	{
		xb = R*cos(2*M_PI/no*i);
		yb = R*sin(2*M_PI/no*i);
		dom.AddSingleParticle(Vec3_t ( xb ,  yb , 0.0 ));
	}

	no = ceil(2*M_PI*Rc/dx);

	for (size_t j=0;j<6;j++)
	{
		R = Rc-sqrt(3.0)/2.0*dx*j;
//		no = ceil(2*M_PI*R/dx);
		for (size_t i=0; i<no; i++)
		{
			xb = R*cos(2*M_PI/no*i+(j%2)*M_PI/no);
			yb = R*sin(2*M_PI/no*i+(j%2)*M_PI/no);
			dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
		}
	}


	dom.Solve(/*tf*/20000.0,/*dt*/maz,/*dtOut*/(100.0*maz),"test06",1500);
	return 0;
}
MECHSYS_CATCH
