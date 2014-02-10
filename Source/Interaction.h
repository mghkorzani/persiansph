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

#ifndef MECHSYS_SPH_INTERACTION_H
#define MECHSYS_SPH_INTERACTION_H


#include <Particle.h>

namespace SPH {

class Interaction
{
public:
    // Constructor
    Interaction (Particle * Pt1, Particle * Pt2,size_t dim, double VisAlpha, double VisBeta, double Vel, double DV);

    // Methods
    void CalcForce      (double dt = 0.0);                ///< Calculates the contact force between particles
    double Kernel       (double r,double h);	          ///< Kernel function
    double GradKernel   (double r,double h);              ///< Gradient of kernel function
    double Pressure     (double Density, double Density0);///< Equation of state for weakly compressible fluid
    double SoundSpeed   (double Density, double Density0);///< Speed of sound in the fluid (dP/drho)
    void ShearStrainCal	(Mat3_t GradVel,Mat3_t SStrain);  ///< Calculate shear strain from velocity gradient

    // Data
    Particle	*P1;			///< Pointer to first particle
    Particle	*P2;			///< Pointer to second particle
    double		alpha;			///< Coefficient of bulk viscosity
    double		beta;			///< Coefficient of Neumann-Richtmyer viscosity
    double		h;				///< Smoothing length
    size_t		Dim;			///< Dimension of the problem
    double		V2;				///< Squared maximum velocity of the fluid for pressure and sound speed
    double		MU;				///< Dynamic Viscosity
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Interaction::Interaction (Particle * Pt1, Particle * Pt2,size_t dim, double VisAlpha, double VisBeta, double Vel, double DV)
{
    P1		=	Pt1;
    P2		=	Pt2;
    h		=	P1->h;                                  ///< It should be revised
    Dim		=	dim;
    alpha	=	VisAlpha;
    beta	=	VisBeta;
    V2		=	Vel*Vel;
    MU		=	DV;
}

inline void Interaction::CalcForce(double dt)
{
	double di = P1->Density;
    double dj = P2->Density;
    double mi = P1->Mass;
    double mj = P2->Mass;
    double Pi,Pj;

    if (P1->IsFree)
    {
    	if (P2->IsFree)
    	{
    	Pi = P1->Pressure = Pressure(di,P1->RefDensity);
    	Pj = P2->Pressure = Pressure(dj,P2->RefDensity);
    	}
    	else
    		{
    		Pi = P1->Pressure = Pressure(di,P1->RefDensity);
    		Pj = P2->Pressure;
    		}
    }
    else
    	{
    	Pj = P2->Pressure = Pressure(dj,P2->RefDensity);
    	Pi = P1->Pressure;
    	}

    Vec3_t vij = P1->v - P2->v;
    Vec3_t rij = P1->x - P2->x;

    if (MU!=0.0)
    {
    Mat3_t Kdelta;
    Identity(Kdelta);

    Mat3_t Sigmai,Sigmaj;
    Mat3_t GradientVelocity,ShearStrain;
    Vec3_t temp;


    Dyad((mj/di),-1*vij,GradKernel(norm(rij),h)*(rij/norm(rij)),GradientVelocity);
    ShearStrainCal(GradientVelocity,ShearStrain);
    Sigmaj			 = Pj*Kdelta-MU*ShearStrain;

    Dyad((mi/dj),vij,-1*GradKernel(norm(rij),h)*(rij/norm(rij)),GradientVelocity);
    ShearStrainCal(GradientVelocity,ShearStrain);
    Sigmaj			 = Pj*Kdelta-MU*ShearStrain;

    Mult (((1/(di*di))*Sigmai+(1/(dj*dj))*Sigmaj),GradKernel(norm(rij),h)*(rij/norm(rij)),temp);
    P1->a			+= -1*mj*temp;                     ///<(2.73) Li, Liu Book
    P1->dDensity	+= (di*mj/dj)*dot(vij,(rij/norm(rij)))*GradKernel(norm(rij),h);                                  ///<(2.58) Li, Liu Book

    P2->a			-= -1*mi*temp;
    P2->dDensity	+= (dj*mi/di)*dot(vij,(rij/norm(rij)))*GradKernel(norm(rij),h);
    }
    else
    {
    	double MUij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    	double Cij	= 0.5*(SoundSpeed(di,P1->RefDensity)+SoundSpeed(dj,P1->RefDensity));
    	double PIij;
    	if (dot(vij,rij)<0) PIij = (-alpha*Cij*MUij+beta*MUij*MUij)/(0.5*(di+dj));                          ///<(2.74) Li, Liu Book
    	else                PIij = 0.0;

    	P1->a			+= -1*mj*(Pi/(di*di)+Pj/(dj*dj)+PIij)*GradKernel(norm(rij),h)*(rij/norm(rij));                     ///<(2.73) Li, Liu Book
        P1->dDensity	+= (di*mj/dj)*dot(vij,(rij/norm(rij)))*GradKernel(norm(rij),h);                                  ///<(2.58) Li, Liu Book

        P2->a			-= -1*mi*(Pi/(di*di)+Pj/(dj*dj)+PIij)*GradKernel(norm(rij),h)*(rij/norm(rij));
        P2->dDensity	+= (dj*mi/di)*dot(vij,(rij/norm(rij)))*GradKernel(norm(rij),h);
    }
}

inline double Interaction::Kernel(double r,double h)
{
	double C;
	switch (Dim)
    {case 1:
       C = 2.0/(3.0*h);
       break;
    case 2:
       C = 10.0/(7.0*h*h*M_PI);
       break;
    case 3:
       C = 1.0/(h*h*h*M_PI);
       break;
    default:
       std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
       C = 1.0/(h*h*h*M_PI);
       break;
    }

    double q = r/h;
    if ((q>=0.0)&&(q<1)) return C*(1-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
    else if (q<=2)       return C*((1.0/4.0)*(2-q)*(2-q)*(2-q));
    else                 return 0.0;
}

inline double Interaction::GradKernel(double r, double h)
{
	double C;
	switch (Dim)
    {case 1:
       C = 2.0/(3.0*h*h);
       break;
    case 2:
       C = 10.0/(7.0*h*h*h*M_PI);
       break;
    case 3:
       C = 1.0/(h*h*h*h*M_PI);
       break;
    default:
       std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
       C = 1.0/(h*h*h*h*M_PI);
       break;
    }
    double q = r/h;
    if ((q>=0.0)&&(q<1)) return C*(-3.0*q+(9.0/4.0)*q*q);
    else if (q<=2)       return C*(-1*(3.0/4.0)*(2-q)*(2-q));
    else                 return 0.0;
}

inline double Interaction::Pressure(double Density, double Density0)
{
	return (100*Density0*V2/7)*(pow(Density/Density0,7)-1);
}

inline double Interaction::SoundSpeed(double Density, double Density0)
{
	return sqrt(7*(100*Density0*V2/7)*(pow(Density/Density0,6)/Density0));
}

inline void Interaction::ShearStrainCal (Mat3_t GradVel,Mat3_t SStrain)
{
    SStrain(0,0)=2*GradVel(0,0);  				SStrain(0,1)=GradVel(1,0)+GradVel(0,1);  	SStrain(0,2)=GradVel(2,0)+GradVel(0,2);
    SStrain(1,0)=GradVel(1,0)+GradVel(0,1);  	SStrain(1,1)=2*GradVel(1,1);				SStrain(1,2)=GradVel(2,1)+GradVel(1,2);
    SStrain(2,0)=GradVel(2,0)+GradVel(0,2);		SStrain(2,1)=GradVel(2,1)+GradVel(1,2);		SStrain(2,2)=2*GradVel(2,2);
}


}; // namespace SPH

#endif // MECHSYS_SPH_INTERACTION_H
