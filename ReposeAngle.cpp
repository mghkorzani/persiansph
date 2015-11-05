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

using std::cout;
using std::endl;

// psi*sqrt(E/(rho*h*h))
double check = 0.0;
double Damp;
void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<0.025)
    #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
		if (domi.Particles[i]->IsFree) domi.Particles[i]->a -= Damp * domi.Particles[i]->v;
	if (domi.Time>0.025 && domi.Time<(0.025+1.0*domi.deltat) && check == 0.0)
	{
		domi.TI	= 0.5;
		domi.TIn= 2.55;
//		std::cout<<"yes"<<std::endl;
//		domi.DelParticles(3);
		check = 1.0;
	}
}

void UserBC(SPH::Domain & domi)
{
	if (domi.Time>(0.025+1.0*domi.deltat))
	{
		#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
		for (size_t i=0; i<domi.FixedParticles.Size(); i++)
		{
			if (domi.Particles[domi.FixedParticles[i]]->x(0)>=0.2)
			{
				domi.Particles[domi.FixedParticles[i]]->Sigma(0,0) = 0.0;
				domi.Particles[domi.FixedParticles[i]]->Sigma(1,1) = 0.0;
				domi.Particles[domi.FixedParticles[i]]->Sigma(2,2) = 0.0;
			}
		}
	}


}

int main(int argc, char **argv) try
{
        SPH::Domain		dom;

        dom.Dimension	= 2;
        dom.Nproc		= 24;
    	dom.PresEq		= 0;
    	dom.KernelType	= 0;
//    	dom.Shepard		= true;

//    	dom.TI			= 0.3;
//    	dom.TIn			= 4.0;
    	dom.Alpha		= 0.1;
//    	dom.Beta		= 0.1;
//    	dom.XSPH		= 0.5;
    	dom.Gravity		= 0.0, -9.81, 0.0;

        double dx,h,rho,K,G;
    	double E,Nu,H,n,L;

    	H	= 0.1;
    	n	= 25.0;
    	L	= 0.2;

    	rho	= 1850.0;
    	E	= 20.0e6;
    	Nu	= 0.2;
    	K	= E/(3.0*(1.0-2.0*Nu));
    	G	= E/(2.0*(1.0+Nu));
    	dx	= H / n;
    	h	= dx*1.5;
        dom.Cs			= sqrt(K/rho);
//        dom.Cs			= sqrt(E/rho);
    	dom.InitialDist	= dx;

    	Damp = 0.02*sqrt(E/(rho*h*h));

        double timestep;
        timestep = (0.05*h/(dom.Cs));

//        dom.GeneralBefore = & UserBC;
        dom.GeneralAfter = & UserDamping;

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = "<<dom.Cs<<endl;
        cout<<"G = "<<G<<endl;
        cout<<"K = "<<K<<endl;
        cout<<"h = "<<h<<endl;

     	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx ,  0.0    , 0.0 ), L + 6.0*dx + dx/10.0 , H + dx/10.0      ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx , -3.0*dx , 0.0 ), 3.5*L + dx/10.0      , 3.0*dx + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

		double K0 = 1-sin(30.0/180.0*M_PI);

     	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G			= G;
    		dom.Particles[a]->K			= K;
    		dom.Particles[a]->Material	= 3;
    		dom.Particles[a]->Fail		= 3;
    		dom.Particles[a]->c			= 0.0;
    		dom.Particles[a]->phi		= 30.0/180.0*M_PI;
    		dom.Particles[a]->psi		= 1.0 /180.0*M_PI;
    		if (dom.Particles[a]->x(1)<0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 4;
    		}
    		if (dom.Particles[a]->x(0)<0.0 && dom.Particles[a]->x(1)>0.0)
    		{
    			dom.Particles[a]->NoSlip	= false;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 2;
    		}
    		if (dom.Particles[a]->x(0)>L && dom.Particles[a]->x(1)>0.0)
    		{
    			dom.Particles[a]->NoSlip	= false;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 3;
    		}

    		dom.Particles[a]->Sigma(1,1)	= rho * -9.81 *(H-dom.Particles[a]->x(1));
    		dom.Particles[a]->Sigmab(1,1)	= rho * -9.81 *(H-dom.Particles[a]->x(1));
    		dom.Particles[a]->Sigma(0,0)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
    		dom.Particles[a]->Sigmab(0,0)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
    		dom.Particles[a]->Sigma(2,2)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
    		dom.Particles[a]->Sigmab(2,2)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
    	}
     			dom.DelParticles(3);

    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
        return 0;
}
MECHSYS_CATCH
