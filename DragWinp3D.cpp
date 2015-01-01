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
    infile >> Ref;		infile.ignore(200,'\n');
    infile >> Csf;		infile.ignore(200,'\n');
    infile >> P0f;		infile.ignore(200,'\n');
    infile >> N;		infile.ignore(200,'\n');
    infile >> Tf;		infile.ignore(200,'\n');
    infile >> Rcf;		infile.ignore(200,'\n');
    infile >> Visf;		infile.ignore(200,'\n');
    infile >> Kf;		infile.ignore(200,'\n');

    SPH::Domain		dom;
	dom.Dimension	= 3;

	dom.BC.Periodic[1] = true;
	dom.RigidBody	= true;
	dom.RBTag		= 4;

	dom.MU			= 1.002e-3;
	dom.PresEq		= 0;
	dom.VisEq		= abs(Visf);
	dom.KernelType	= abs(Kf);
	dom.Nproc		= abs(N);

//	dom.TI			= 0.05;

	double xb,yb,zb,h,rho,mass,U;
	double dx,R,Rc,Re;
	size_t no,no1;

	rho = 998.21;
	dx = 0.002;
	h = dx*1.1;
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

	dom.AddBoxLength(3 ,Vec3_t ( -5.0*Rc , -5.0*Rc , -5.0*Rc ), 15.0*Rc , 10.0*Rc  ,  10.0*Rc , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		zb=dom.Particles[a]->x(2);
		if ((xb*xb+yb*yb+zb*zb)<((Rc+h/2.0)*(Rc+h/2.0)))
		{
			dom.Particles[a]->ID=4;
			dom.Particles[a]->IsFree=false;
		}
	}
	dom.DelParticles(4);

	for (size_t j=0;j<6;j++)
	{
		R = Rc-dx*j;
		no = ceil(2*M_PI*R/dx);
		no1 = ceil(M_PI*R/dx);
		for (size_t k=0; k<=no1; k++)
		for (size_t i=0; i<no; i++)
		{
			xb = R*cos(2*M_PI/no*i)*sin(M_PI/no1*k);
			yb = R*sin(2*M_PI/no*i)*sin(M_PI/no1*k);
			zb = R*cos(M_PI/no1*k);
			dom.AddSingleParticle(4,Vec3_t ( xb ,  yb , zb ), mass , rho , h , true);
		}
	}
//	dom.WriteXDMF("maz");
	dom.Solve(/*tf*/20000.0,/*dt*/maz,/*dtOut*/(500.0*maz),"test06");
	return 0;
}
MECHSYS_CATCH
