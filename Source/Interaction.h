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
    Interaction (Particle * Pt1, Particle * Pt2,size_t dim, double VisAlpha, double VisBeta, double Vel, double DV, double XSPHfac);

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
    double 		X;				///< Factor of XSPH
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Interaction::Interaction (Particle * Pt1, Particle * Pt2,size_t dim, double VisAlpha, double VisBeta, double Vel, double DV, double XSPHfac)
{
    P1		=	Pt1;
    P2		=	Pt2;
    h		=	P1->h;                                  ///< It should be revised
    Dim		=	dim;
    alpha	=	VisAlpha;
    beta	=	VisBeta;
    V2		=	Vel*Vel;
    MU		=	DV;
    X		=	XSPHfac;

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

    Mat3_t Kdelta,PI,RF;
    Mat3_t Sigmai,Sigmaj;
    Mat3_t GradientVelocity,ShearStrain;
    GradientVelocity=0.0;
    ShearStrain=0.0;

    Vec3_t temp;
	Mat3_t A,B;

    Kdelta = 	1.0, 0.0, 0.0,
    			0.0, 1.0, 0.0,
    			0.0, 0.0, 1.0;

    double TI;
//    TI= ((Pi/(di*di)+Pj/(dj*dj))*0.01)*(Kernel(norm(rij),h)/Kernel(0.002,h))*(Kernel(norm(rij),h)/Kernel(0.002,h))*(Kernel(norm(rij),h)/Kernel(0.002,h))*(Kernel(norm(rij),h)/Kernel(0.002,h));
    TI=0.0;


    double MUij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);                                                ///<(2.75) Li, Liu Book
	double Cij	= 0.5*(SoundSpeed(di,P1->RefDensity)+SoundSpeed(dj,P1->RefDensity));
	double PIij;
	if (dot(vij,rij)<0) PIij = (-alpha*Cij*MUij+beta*MUij*MUij)/(0.5*(di+dj));                          ///<(2.74) Li, Liu Book
	else                PIij = 0.0;

    Dyad((mj/di),-vij,GradKernel(norm(rij),h)*(rij/norm(rij)),GradientVelocity);
    ShearStrainCal(GradientVelocity,ShearStrain);
	Mult((Pi/(di*di)),Kdelta,A);
	Mult((-MU/(di*di)),ShearStrain,B);
	Add (A,B,Sigmai);

    GradientVelocity=0.0;
    ShearStrain=0.0;

    Dyad((mi/dj),vij,-GradKernel(norm(rij),h)*(rij/norm(rij)),GradientVelocity);
    ShearStrainCal(GradientVelocity,ShearStrain);
	Mult((Pj/(dj*dj)),Kdelta,A);
	Mult((-MU/(dj*dj)),ShearStrain,B);
	Add (A,B,Sigmaj);

	Mult(PIij,Kdelta,PI);
	Mult(TI,Kdelta,RF);

	Add (Sigmai,Sigmaj,A);
	Add (A,PI,B);
	Add (B,RF,A);


    Mult (A,(GradKernel(norm(rij),h)*(rij/norm(rij))),temp);

    omp_set_lock(&P1->my_lock);
	P1->a			+= -mj*temp;                     ///<(2.73) Li, Liu Book
    P1->dDensity	+= (mj)*dot(vij,(rij/norm(rij)))*GradKernel(norm(rij),h);                                  ///<(2.58) Li, Liu Book
    P1->VXSPH		+= X*mj/(0.5*(di+dj))*Kernel(norm(rij),h)*-vij;
    omp_unset_lock(&P1->my_lock);

    omp_set_lock(&P2->my_lock);
    P2->a			-= -mi*temp;
    P2->dDensity	+= (mi)*dot(vij,(rij/norm(rij)))*GradKernel(norm(rij),h);
    P2->VXSPH		+= X*mi/(0.5*(di+dj))*Kernel(norm(rij),h)*vij;
    omp_unset_lock(&P2->my_lock);
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
