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
    		Interaction	(Particle * Pt1, Particle * Pt2,size_t Dim0, double Alpha0, double Beta0, double Cs0, double MU0, double XSPHC0, double TIC0, double InitialDist0, double P00, int PresEq0);
    // Methods
    void	CalcForce	(double dt);						///< Calculates the contact force between particles
    double	Kernel		(double r,double h);				///< Kernel function
    double	GradKernel	(double r,double h);				///< Gradient of the kernel function
    double	Pressure	(double Density, double Density0);	///< Equation of state for a weakly compressible fluid
    double	SoundSpeed	(double Density, double Density0);	///< Speed of sound in a fluid (dP/drho)

    // Data
    Particle	*P1;			///< Pointer to the first particle
    Particle	*P2;			///< Pointer to the second particle

    double		alpha;			///< Coefficient of the bulk viscosity
    double		beta;			///< Coefficient of Neumann-Richtmyer viscosity
    double		MU;				///< Dynamic Viscosity coefficient

    double		h;				///< Smoothing length
    size_t		Dim;			///< Dimension of problem
    double		Cs;				///< Speed of Sound

    double 		XSPHC;			///< XSPH coefficient
    double		TIC;			///< Monagham Tensile Instability coefficient
    double		InitialDist;	///< Initial distance of particles for Tensile Instability calculation

    double		P0;				///< Background pressure
    int			PresEq;			///< Selection function for the various equation of state
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Interaction::Interaction (Particle * Pt1, Particle * Pt2,size_t Dim0, double Alpha0, double Beta0, double Cs0, double MU0, double XSPHC0, double TIC0, double InitialDist0, double P00, int PresEq0)
{
    P1		= Pt1;
    P2		= Pt2;

    h		= P1->h;                                  ///< It should be revised
    Dim		= Dim0;
    Cs		= Cs0;

    alpha	= Alpha0;
    beta	= Beta0;
    MU		= MU0;

    XSPHC	= XSPHC0;
    TIC		= TIC0;
    InitialDist = InitialDist0;

    P0		= P00;
    PresEq	= PresEq0;
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

    //Artificial Viscosity
    double PIij = 0.0;
    if (alpha!=0.0 || beta!=0.0)
    {
    	double MUij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    	double Cij = 0.5*(SoundSpeed(di,P1->RefDensity)+SoundSpeed(dj,P1->RefDensity));
    	if (dot(vij,rij)<0) PIij = (-alpha*Cij*MUij+beta*MUij*MUij)/(0.5*(di+dj));                          ///<(2.74) Li, Liu Book
    	else                PIij = 0.0;
    }

    //Tensile Instability
    double TI = 0.0;
    if ((TIC > 0.0) && (Pi < 0.0) && (Pj < 0.0))
    {
        double pa,pb;
    	pa = abs(Pi);
    	pb = abs(Pj);
        TI = TIC*(pa/(di*di)+pb/(dj*dj))*pow((Kernel(norm(rij),h)/Kernel(InitialDist,h)),4);
    }
    else TI = 0.0;

    //Real Viscosity
    Vec3_t VI = 0.0;
    if (MU!=0.0) VI = 8*MU/((di+dj)*(di+dj)*dot(rij,rij))*dot(rij,GradKernel(norm(rij),h)*(rij/norm(rij)))*vij;

    if (XSPHC != 0.0)
    {
        omp_set_lock(&P1->my_lock);
        P1->VXSPH		+= XSPHC*mj/(0.5*(di+dj))*Kernel(norm(rij),h)*-vij;
        omp_unset_lock(&P1->my_lock);

        omp_set_lock(&P2->my_lock);
		P2->VXSPH		+= XSPHC*mi/(0.5*(di+dj))*Kernel(norm(rij),h)*vij;
		omp_unset_lock(&P2->my_lock);
    }

    omp_set_lock(&P1->my_lock);
	P1->a			+= -mj*(Pi/(di*di)+Pj/(dj*dj)+PIij+TI)*GradKernel(norm(rij),h)*(rij/norm(rij))+ mj*VI;
    P1->dDensity	+= (di*mj/dj)*dot((vij+P1->VXSPH-P2->VXSPH),(rij/norm(rij)))*GradKernel(norm(rij),h);
    omp_unset_lock(&P1->my_lock);


    omp_set_lock(&P2->my_lock);
    P2->a			-= -mi*(Pi/(di*di)+Pj/(dj*dj)+PIij+TI)*GradKernel(norm(rij),h)*(rij/norm(rij))+ mi*VI;
    P2->dDensity	+= (dj*mi/di)*dot((-vij+P2->VXSPH-P1->VXSPH),(-rij/norm(rij)))*GradKernel(norm(rij),h);
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
	switch (PresEq)
    {
	case 0:
		return P0+(Cs*Cs)*(Density-Density0);
		break;
	case 1:
		return P0+(Density0*Cs*Cs/7)*(pow(Density/Density0,7)-1);
		break;
	case 2:
		return (Cs*Cs)*Density;
		break;
	default:
		std::cout << "Please correct Pressure Equation No, otherwise Eq. (0) is used" << std::endl;
		return P0+(Cs*Cs)*(Density-Density0);
		break;
    }
}

inline double Interaction::SoundSpeed(double Density, double Density0)
{
	switch (PresEq)
    {
	case 0:
		return Cs;
		break;
	case 1:
		return sqrt((Cs*Cs)*pow(Density/Density0,6));
		break;
	case 2:
		return Cs;
		break;
	default:
		std::cout << "Please correct Pressure Equation No, otherwise Eq. (0) is used" << std::endl;
		return Cs;
		break;
    }
}

}; // namespace SPH

#endif // MECHSYS_SPH_INTERACTION_H
