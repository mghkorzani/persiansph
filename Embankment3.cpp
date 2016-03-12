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

#include "./Source/Domain.h"
#include "./Source/Interaction.h"
#include <fstream>


using std::cout;
using std::endl;

int check = 0;
double DampF,DampS,Cs,u;
double t = 2.0;
double tim = 0.0;
Array< size_t > Forced;
double tout = t,dtout=0.5;
Vec3_t force,loc;
size_t nop;

void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<t)
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.Particles.Size(); i++)
		{
			if (domi.Particles[i]->IsFree && domi.Particles[i]->Material == 3) domi.Particles[i]->a -= DampS * domi.Particles[i]->v;

		}
	if (domi.Time>t)
	{
        if (domi.Time>=tout)
        {
        	if (nop != domi.Particles.Size())
        	{
         		Forced.Clear();
        		for (size_t a=0; a<domi.Particles.Size(); a++)
         			if (domi.Particles[a]->ID == 7)
        				Forced.Push(a);
        		nop = domi.Particles.Size();
        	}

        	force = 0.0;
        	loc = 0.0;
    		for (size_t i=0; i<Forced.Size(); i++)
    		{
    			force += domi.Particles[Forced[i]]->a*domi.Particles[Forced[i]]->Mass;
    			loc += domi.Particles[Forced[i]]->x;
    		}
    		std::fstream Force ("Force.txt", std::ios::out | std::ios::app );
    		std::fstream X ("Loc.txt", std::ios::out | std::ios::app );
    		Force << force(1) <<std::endl;
    		X << (loc(1)/Forced.Size()) <<std::endl;
    		Force.close();
    		X.close();
            tout += dtout;
        }
		for (size_t i=0; i<Forced.Size(); i++)
			domi.Particles[Forced[i]]->a = Vec3_t(0.0 , -0.0001 , 0.0);
	}
}

int main(int argc, char **argv) try
{
    SPH::Domain		dom;
	dom.Dimension	= 2;

	dom.KernelType	= 0;
	dom.SeepageType = 2;
	dom.VisEq		= 0;
	dom.Nproc		= 8;
	dom.Excemption	= 4;
//	dom.TimestepConstrain1 = false;

//	dom.XSPH		= 0.5;
	dom.Scheme		= 0;
	dom.Gravity		= 0.0,-9.81,0.0;
	dom.GeneralAfter= & UserDamping;

	double xb,yb,h,dx;

	dx		= 0.125;
	h		= dx*1.3;
	dom.InitialDist	= dx;

	double Nu,CsS,Ts;
	double E1,K1,G1,CsS1,rhoS1,c1,Phi1,Psi1;
	double E2,K2,G2,CsS2,rhoS2,c2,Phi2,Psi2;

	Nu		= 0.3;

	E1		= 25.0e6;
	K1		= E1/(3.0*(1.0-2.0*Nu));
	G1		= E1/(2.0*(1.0+Nu));
	rhoS1	= 18.0e3/9.81;
    CsS1	= sqrt(K1/rhoS1);
    c1		= 6.0e3;
    Phi1	= 30.0;
    Psi1	= 0.0;

	E2		= 50.0e6;
	K2		= E2/(3.0*(1.0-2.0*Nu));
	G2		= E2/(2.0*(1.0+Nu));
	rhoS2	= 20.0e3/9.81;
    CsS2	= sqrt(K2/rhoS2);
    c2		= 10.0e3;
    Phi2	= 20.0;
    Psi2	= 0.0;

    CsS		= std::max(CsS1,CsS2);
    Ts		= (0.25*h/CsS);
  	DampS	= 0.02*sqrt(E1/(rhoS1*h*h));

	dom.AddBoxLength(3 ,Vec3_t (-12.5 - 3.0*dx ,-3.0 - 3.0*dx , 0.0 ), 25.0 + 7.0*dx + dx/10.0 , 7.0 + 7.0*dx + dx/10.0  ,  0 , dx/2.0 ,rhoS1, h,1 , 0 , false,false);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		if (dom.Particles[a]->ID==3)
		{
//			dom.Particles[a]->Shepard	= true;
			dom.Particles[a]->Material	= 3;
			dom.Particles[a]->Alpha		= 0.1;
			dom.Particles[a]->Beta		= 0.1;
			dom.Particles[a]->TI		= 0.5;
			dom.Particles[a]->TIn		= 2.55;

			xb=dom.Particles[a]->x(0);
			yb=dom.Particles[a]->x(1);
			if (yb<-3.0)
			{
				dom.Particles[a]->ID	= 5;
				dom.Particles[a]->IsFree= false;
				dom.Particles[a]->NoSlip= true;
			}
			if (yb>(-(1.0/1.5)*(xb-7.5)) && yb>0.0 && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 8;
			if (yb>( (1.0/1.5)*(xb+7.5)) && yb>0.0 && dom.Particles[a]->ID == 3)
				dom.Particles[a]->ID	= 8;
			if ((xb>12.5 || xb<-12.5) && dom.Particles[a]->ID == 3)
			{
				dom.Particles[a]->ID	= 6;
				dom.Particles[a]->IsFree= false;
				dom.Particles[a]->NoSlip= false;
			}

			if (yb<=0.0)
			{
				dom.Particles[a]->Density	= rhoS2;
				dom.Particles[a]->Densitya	= rhoS2;
				dom.Particles[a]->Densityb	= rhoS2;
				dom.Particles[a]->RefDensity= rhoS2;
				dom.Particles[a]->Cs		= CsS2;
				dom.Particles[a]->G			= G2;
				dom.Particles[a]->K			= K2;
				dom.Particles[a]->Fail		= 3;
				dom.Particles[a]->c			= c2;
				dom.Particles[a]->phi		= Phi2/180.0*M_PI;
				dom.Particles[a]->psi		= Psi2/180.0*M_PI;
			}
			else
			{
				if (dom.Particles[a]->ID == 3)
					dom.Particles[a]->ID	= 4;
				dom.Particles[a]->Cs		= CsS1;
				dom.Particles[a]->G			= G1;
				dom.Particles[a]->K			= K1;
				dom.Particles[a]->Fail		= 3;
				dom.Particles[a]->c			= c1;
				dom.Particles[a]->phi		= Phi1/180.0*M_PI;
				dom.Particles[a]->psi		= Psi1/180.0*M_PI;
			}
			if (dom.Particles[a]->ID == 4 && (xb>1.0 || xb<-1.0) && yb>=4.0)
				dom.Particles[a]->ID	= 8;
		}
	}
	dom.DelParticles(8);

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		xb=dom.Particles[a]->x(0);
		yb=dom.Particles[a]->x(1);
		if (dom.Particles[a]->ID == 4 && xb<1.0 && xb>-1.0 && yb>=4.0)
		{
			Forced.Push(a);
			dom.Particles[a]->ID	= 7;
		}
	}
	nop = dom.Particles.Size();
	std::cout<<Forced<<std::endl;


	dom.Solve(/*tf*/50000.0,/*dt*/Ts,/*dtOut*/0.5,"test06",2000);
	return 0;
}
MECHSYS_CATCH
