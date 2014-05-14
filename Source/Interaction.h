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

    	Pi = P1->Pressure = Pressure(di,P1->RefDensity);
    	Pj = P2->Pressure = Pressure(dj,P2->RefDensity);

    Vec3_t vij = P1->v - P2->v;
    Vec3_t rij = P1->x - P2->x;

    double MUij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    double Cij = 0.5*(SoundSpeed(di,P1->RefDensity)+SoundSpeed(dj,P1->RefDensity));
    double PIij;
    if (dot(vij,rij)<0) PIij = (-alpha*Cij*MUij+beta*MUij*MUij)/(0.5*(di+dj));                          ///<(2.74) Li, Liu Book
    else                PIij = 0.0;

    omp_set_lock(&P2->my_lock);
    P2->VXSPH		+= X*mi/(0.5*(di+dj))*Kernel(norm(rij),h)*vij;
    omp_unset_lock(&P2->my_lock);


//    double pa,pb;
//    if (Pi<0) pa = abs(Pi);
//    		else pa=0;
//    if (Pj<0) pa = abs(Pj);
//    		else pb=0;
//    double TI = 0.3*(pa/(di*di)+pb/(dj*dj))*pow((Kernel(norm(rij),h)/Kernel(0.0000225,h)),4);
    double TI=0.0;

    omp_set_lock(&P1->my_lock);
    P1->VXSPH		+= X*mj/(0.5*(di+dj))*Kernel(norm(rij),h)*-vij;
	P1->a			+= -mj*(Pi/(di*di)+Pj/(dj*dj)+PIij+TI)*GradKernel(norm(rij),h)*(rij/norm(rij))+ mj*8*MU/((di+dj)*(di+dj)*dot(rij,rij))*dot(rij,GradKernel(norm(rij),h)*(rij/norm(rij)))*vij;  ///<(2.73) Li, Liu Book
    P1->dDensity	+= (mj)*dot((vij+P1->VXSPH-P2->VXSPH),(rij/norm(rij)))*GradKernel(norm(rij),h);                                  ///<(2.58) Li, Liu Book
    omp_unset_lock(&P1->my_lock);


    omp_set_lock(&P2->my_lock);
    P2->a			-= -mi*(Pi/(di*di)+Pj/(dj*dj)+PIij+TI)*GradKernel(norm(rij),h)*(rij/norm(rij))+ mi*8*MU/((di+dj)*(di+dj)*dot(rij,rij))*dot(rij,GradKernel(norm(rij),h)*(rij/norm(rij)))*vij;
    P2->dDensity	+= (mi)*dot((-vij+P2->VXSPH-P1->VXSPH),(-rij/norm(rij)))*GradKernel(norm(rij),h);
    omp_unset_lock(&P2->my_lock);
}

//inline double Interaction::Kernel(double r,double h)
//{
//	double C;
//	switch (Dim)
//    {case 1:
//       C = 1.0/(7.0*h);
//       break;
//    case 2:
//       C = 1.0/(3.0*h*h*M_PI);
//       break;
//    case 3:
//       C = 15.0/(62.0*h*h*h*M_PI);
//       break;
//    default:
//       std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
//       C = 1.0/(h*h*h*M_PI);
//       break;
//    }
//
//    double q = r/h;
//    if ((q>=0.0)&&(q<1)) return C*(6.0-6.0*q+(q*q*q));
//    else if (q<=2)       return C*((2-q)*(2-q)*(2-q));
//    else                 return 0.0;
//}
//
//inline double Interaction::GradKernel(double r, double h)
//{
//	double C;
//	switch (Dim)
//    {case 1:
//       C = 1.0/(7.0*h*h);
//       break;
//    case 2:
//       C = 1.0/(3.0*h*h*h*M_PI);
//       break;
//    case 3:
//       C = 15.0/(62.0*h*h*h*h*M_PI);
//       break;
//    default:
//       std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
//       C = 1.0/(h*h*h*h*M_PI);
//       break;
//    }
//    double q = r/h;
//    if ((q>=0.0)&&(q<1)) return C*(-6.0+3.0*q*q);
//    else if (q<=2)       return C*(-3.0*(2-q)*(2-q));
//    else                 return 0.0;
//}


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
//		return (100*Density0*V2/7)*(pow(Density/Density0,7)-1);
		return 3000+10*10*(Density-Density0);
}

inline double Interaction::SoundSpeed(double Density, double Density0)
{
		return sqrt(7*(100*Density0*V2/7)*(pow(Density/Density0,6)/Density0));
}

}; // namespace SPH

#endif // MECHSYS_SPH_INTERACTION_H
