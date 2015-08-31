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
    double N;
    double Tf;
    double Rcf;
    double Kf;
    double Visf;
    double DomSizex;
    double DomSizey;
    double ax;
    double hf;
    double dxf;
    double Nop;
    infile >> Ref;		infile.ignore(200,'\n');
    infile >> Csf;		infile.ignore(200,'\n');
    infile >> P0f;		infile.ignore(200,'\n');
    infile >> N;		infile.ignore(200,'\n');
    infile >> Tf;		infile.ignore(200,'\n');
    infile >> Rcf;		infile.ignore(200,'\n');
    infile >> Visf;		infile.ignore(200,'\n');
    infile >> Kf;		infile.ignore(200,'\n');
    infile >> DomSizex;	infile.ignore(200,'\n');
    infile >> DomSizey;	infile.ignore(200,'\n');
    infile >> ax;		infile.ignore(200,'\n');
    infile >> hf;		infile.ignore(200,'\n');
    infile >> dxf;		infile.ignore(200,'\n');
    infile >> Nop;		infile.ignore(200,'\n');

    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.BC.Periodic[1] = true;
	dom.RigidBody	= true;
	dom.RBTag		= 4;
//	dom.NoSlip		= true;

	dom.MU			= 1.002e-3;
//	dom.MU			= 2.0e-5;
	dom.PresEq		= 0;
	dom.VisEq		= abs(Visf);
	dom.KernelType	= abs(Kf);
	dom.Nproc		= abs(N);

//	dom.TI			= 0.05;

	double xb,yb,h,rho,mass,U;
	double dx,R,Rc,Re;
	size_t no;

	rho = 998.21;
	dx = dxf;
	h = dx*hf;
	Rc = Rcf;
	mass = dx*dx*rho;
	Re = Ref;
	U = Re*dom.MU/(rho*2.0*(Rc+dx/2.0));

	dom.BC.InOutFlow =3;
	dom.BC.allv = U,0.0,0.0;
	dom.BC.inDensity = rho;
	dom.BC.inv = U,0.0,0.0;
	dom.BC.outv = U,0.0,0.0;
	dom.BC.outDensity = rho;
	dom.BC.MassConservation = true;

	dom.Cs				= U*Csf;
	dom.P0				= dom.Cs*dom.Cs*rho*P0f;
	dom.InitialDist 	= dx;
	double maz = (Tf*h/(dom.Cs+U));

	std::cout<<"Re = "<<Re<<std::endl;
	std::cout<<"V  = "<<U<<std::endl;
	std::cout<<"Cs = "<<dom.Cs<<std::endl;
	std::cout<<"P0 = "<<dom.P0<<std::endl;
	std::cout<<"Time Step = "<<maz<<std::endl;
	std::cout<<"Resolution = "<<(2.0*(Rc+dx/2.0)/dx)<<std::endl;

	dom.AddBoxLength(3 ,Vec3_t ( -ax*Rc , -DomSizey/2.0*Rc , 0.0 ), DomSizex*Rc , DomSizey*Rc  ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if ((xb*xb+yb*yb)<((Rc+h/2.0)*(Rc+h/2.0)))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

//	no = 3*ceil(2*M_PI*Rc/dx);
//	for (size_t i=0; i<no; i++)
//	{
//		xb = R*cos(2*M_PI/no*i);
//		yb = R*sin(2*M_PI/no*i);
//		dom.AddNSSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), true);
//	}


	for (size_t j=0;j<7;j++)
	{
		R = Rc-dx*j;
		no = ceil(2*M_PI*R/dx);
		for (size_t i=0; i<no; i++)
		{
			xb = R*cos(2*M_PI/no*i);
			yb = R*sin(2*M_PI/no*i);
			dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , 0.0 ), mass , rho , h , true);
		}
	}

	dom.Solve(/*tf*/30000000.0,/*dt*/maz,/*dtOut*/(500.0*maz),"test06",abs(Nop));
	return 0;
}
MECHSYS_CATCH
