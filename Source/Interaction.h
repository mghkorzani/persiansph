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
    		Interaction	(Particle * Pt1, Particle * Pt2,size_t Dim0, double Alpha0, double Beta0, double Cs0, double MU0,
    				double XSPHC0, double TIC0, double InitialDist0, double P00, size_t PresEq0, size_t KernelType0, Vec3_t domsize);
    // Methods
    void	CalcForce	(double dt, double CF, bool ShepardFilter);	///< Calculates the contact force between particles
    double	Kernel		(double r,double h);				///< Kernel function
    double	GradKernel	(double r,double h);				///< Gradient of the kernel function
    double LaplaceKernel(double r, double h);				///< Laplacian of the kernel function
    double SecDerivativeKernel(double r, double h);			///< Second derivative of the kernel function
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
    size_t		PresEq;			///< Selection function for the various equation of state
    size_t		KernelType;		///< Type of the kernel which is going to be used
    Vec3_t		DomainSize;		///< Domain size which also show the periodic Boundary as well
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Interaction::Interaction (Particle * Pt1, Particle * Pt2,size_t Dim0, double Alpha0, double Beta0, double Cs0, double MU0,
		double XSPHC0, double TIC0, double InitialDist0, double P00, size_t PresEq0, size_t KernelType0, Vec3_t domsize)
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
    KernelType = KernelType0;
    DomainSize = domsize;
}

inline void Interaction::CalcForce(double dt, double CF, bool ShepardFilter)
{
	double di = P1->Density;
    double dj = P2->Density;
    double mi = P1->Mass;
    double mj = P2->Mass;
    double Pi;
    double Pj;

    Pi = P1->Pressure = Pressure(P1->Density,P1->RefDensity);
    Pj = P2->Pressure = Pressure(P2->Density,P2->RefDensity);

    Vec3_t vij = P1->v - P2->v;
    Vec3_t rij = P1->x - P2->x;

    //Correction of rij for Periodic BC
    if (DomainSize(0)>0.0) {if (rij(0)>2*CF*h || rij(0)<-2*CF*h) {(P1->CC[0]>P2->CC[0]) ? rij(0) -= DomainSize(0) : rij(0) += DomainSize(0);}}
	if (DomainSize(1)>0.0) {if (rij(1)>2*CF*h || rij(1)<-2*CF*h) {(P1->CC[1]>P2->CC[1]) ? rij(1) -= DomainSize(1) : rij(1) += DomainSize(1);}}
	if (DomainSize(2)>0.0) {if (rij(2)>2*CF*h || rij(2)<-2*CF*h) {(P1->CC[2]>P2->CC[2]) ? rij(2) -= DomainSize(2) : rij(2) += DomainSize(2);}}

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
    double TIij = 0.0, TIji = 0.0;
    if ((TIC > 0.0) && (Pi < 0.0) && (Pj < 0.0))
    {
        TIij = TIji= TIC*(-Pi/(di*di)-Pj/(dj*dj))*pow((Kernel(norm(rij),h)/Kernel(InitialDist,h)),4);
        if (!P1->IsFree) TIij = 0.0;
        if (!P2->IsFree) TIji = 0.0;
    }

    //Real Viscosity
    Vec3_t VI = 0.0;
//    if (MU!=0.0) VI = 2.0*MU/(di*dj*norm(rij))*GradKernel(norm(rij),h)*vij;  //Morris et al 1997
//    if (MU!=0.0) VI = 8.0*MU/((di+dj)*(di+dj)*(dot(rij,rij)+0.01*h*h))*dot(rij,GradKernel(norm(rij),h)*(rij/norm(rij)))*vij; //Shao et al 2003
//    if (MU!=0.0) VI = MU/(di*dj)*((Dim+1.0/3.0)*GradKernel(norm(rij),h)/norm(rij)*vij+(dot(vij,rij)/3.0*rij+dot(rij,rij)*vij)/norm(rij)*(-1.0/dot(rij,rij)*GradKernel(norm(rij),h)+1.0/norm(rij)*LaplaceKernel(norm(rij),h)));
//    if (MU!=0.0) VI = MU/(di*dj)*LaplaceKernel(norm(rij),h)*vij;  //Real Viscosity

    // XSPH Monaghan
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
    P1->a += -mj*(Pi/(di*di)+Pj/(dj*dj)+PIij+TIij)*GradKernel(norm(rij),h)*(rij/norm(rij))+ mj*VI;
    P1->Vis += mj*VI;
    if (P1->ct==30 && ShepardFilter)
    {
    	P1->SumDen += mj*Kernel(norm(rij),h);
    	P1->ZWab +=mj/dj*Kernel(norm(rij),h);
    }
    P1->dDensity += (di*mj/dj)*dot((vij+P1->VXSPH-P2->VXSPH),(rij/norm(rij)))*GradKernel(norm(rij),h);
    omp_unset_lock(&P1->my_lock);


    omp_set_lock(&P2->my_lock);
    P2->a -= -mi*(Pi/(di*di)+Pj/(dj*dj)+PIij+TIji)*GradKernel(norm(rij),h)*(rij/norm(rij))+ mi*VI;
    P2->Vis -= mi*VI;
    if (P2->ct==30 && ShepardFilter)
    {
    	P2->SumDen += mi*Kernel(norm(rij),h);
    	P2->ZWab +=mi/di*Kernel(norm(rij),h);
    }
    P2->dDensity += (dj*mi/di)*dot((-vij+P2->VXSPH-P1->VXSPH),(-rij/norm(rij)))*GradKernel(norm(rij),h);
    omp_unset_lock(&P2->my_lock);
}

inline double Interaction::Kernel(double r,double h)
{
	double C;
	double q = r/h;
    if (Dim<=1 || Dim>3)
    {
    	std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
    	Dim = 3;
    }
	switch (KernelType)
    {case 0:
    	// Qubic Spline
    	Dim ==2 ? C = 10.0/(7.0*h*h*M_PI) : C = 1.0/(h*h*h*M_PI);

		if ((q>=0.0)&&(q<1)) return C*(1.0-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
		else if (q<=2)       return C*((1.0/4.0)*(2.0-q)*(2.0-q)*(2.0-q));
		else                 return 0.0;

		break;
    case 1:
    	// Quadratic
    	Dim ==2 ? C = 2.0/(h*h*M_PI) : C = 5.0/(4.0*h*h*h*M_PI);

    	if (q<=2) 			 return C*(3.0/4.0-(3.0/4.0)*q+(3.0/16.0)*q*q);
		else                 return 0.0;

    	break;
    case 2:
    	// Quintic
    	Dim ==2 ? C = 7.0/(4.0*h*h*M_PI) : C = 7.0/(8.0*h*h*h*M_PI);

    	if (q<=2) 			 return C*pow((1.0-q/2.0),4)*(2.0*q+1.0);
		else                 return 0.0;

    	break;
    case 3:
    	// Gaussian with compact support
    	Dim ==2 ? C = 1.0/(h*h*M_PI) : C = 1.0/(h*h*h*pow(M_PI,(3.0/2.0)));

    	if (q<=2) 			 return C*exp(-q*q);
		else                 return 0.0;

    	break;
   default:
		// Qubic Spline
	    std::cout << "Kernel No is out of range so Cubic Spline is used" << std::endl;
		Dim ==2 ? C = 10.0/(7.0*h*h*M_PI) : C = 1.0/(h*h*h*M_PI);

		if ((q>=0.0)&&(q<1)) return C*(1.0-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
		else if (q<=2)       return C*((1.0/4.0)*(2.0-q)*(2.0-q)*(2.0-q));
		else                 return 0.0;

       break;
    }
}

inline double Interaction::GradKernel(double r, double h)
{
	double C;
	double q = r/h;
    if (Dim<=1 || Dim>3)
    {
    	std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
    	Dim = 3;
    }
	switch (KernelType)
    {case 0:
    	// Qubic Spline
    	Dim ==2 ? C = 10.0/(7.0*h*h*h*M_PI) : C = 1.0/(h*h*h*h*M_PI);

        if ((q>=0.0)&&(q<1)) return C*(-3.0*q+(9.0/4.0)*q*q);
        else if (q<=2)       return C*((-3.0/4.0)*(2.0-q)*(2.0-q));
        else                 return 0.0;

		break;
    case 1:
    	// Quadratic
    	Dim ==2 ? C = 2.0/(h*h*h*M_PI) : C = 5.0/(4.0*h*h*h*h*M_PI);

    	if (q<=2) 			 return C*(-3.0/4.0+(3.0/8.0)*q);
		else                 return 0.0;

    	break;
    case 2:
    	// Quintic
    	Dim ==2 ? C = 7.0/(4.0*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*M_PI);

    	if (q<=2) 			 return C*-5.0*q*pow((1.0-q/2.0),3);
		else                 return 0.0;

    	break;
    case 3:
    	// Gaussian with compact support
    	Dim ==2 ? C = 1.0/(h*h*h*M_PI) : C = 1.0/(h*h*h*h*pow(M_PI,(3.0/2.0)));

    	if (q<=2) 			 return -2.0*q*C*exp(-q*q);
		else                 return 0.0;

    	break;
   default:
		// Qubic Spline
	    std::cout << "Kernel No is out of range so Cubic Spline is used" << std::endl;
    	Dim ==2 ? C = 10.0/(7.0*h*h*h*M_PI) : C = 1.0/(h*h*h*h*M_PI);

        if ((q>=0.0)&&(q<1)) return C*(-3.0*q+(9.0/4.0)*q*q);
        else if (q<=2)       return C*(-1.0*(3.0/4.0)*(2.0-q)*(2.0-q));
        else                 return 0.0;

       break;
    }
}

inline double Interaction::LaplaceKernel(double r, double h)
{
	double C;
	double q = r/h;
    if (Dim<=1 || Dim>3)
    {
    	std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
    	Dim = 3;
    }
	switch (KernelType)
    {case 0:
    	// Qubic Spline
    	Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

        if ((q>=0.0)&&(q<1)) return C*(-3.0+(9.0/2.0)*q)+(Dim-1)/q*C*(-3.0*q+(9.0/4.0)*q*q);
        else if (q<=2)       return C*((3.0/2.0)*(2.0-q))+(Dim-1)/q*C*((-3.0/4.0)*(2.0-q)*(2.0-q));
        else                 return 0.0;

		break;
    case 1:
    	// Quadratic
    	Dim ==2 ? C = 2.0/(h*h*h*h*M_PI) : C = 5.0/(4.0*h*h*h*h*h*M_PI);

    	if (q<=2) 			 return C*(-3.0/4.0)+(Dim-1)/q*C*(-3.0/4.0+(3.0/8.0)*q);
		else                 return 0.0;

    	break;
    case 2:
    	// Quintic
    	Dim ==2 ? C = 7.0/(4.0*h*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*h*M_PI);

    	if (q<=2) 			 return -3.0*C*pow((1.0-q/2.0),2)*(1.0+3.0*q-3.0*q*q)+(Dim-1)/q*C*-5.0*q*pow((1.0-q/2.0),3);
		else                 return 0.0;

    	break;
    case 3:
    	// Gaussian with compact support
    	Dim ==2 ? C = 1.0/(h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*pow(M_PI,(3.0/2.0)));

    	if (q<=2) 			 return C*2.0*(2.0*q*q-1.0)*exp(-q*q)+(Dim-1)/q*-2.0*q*C*exp(-q*q);
		else                 return 0.0;

    	break;
   default:
		// Qubic Spline
	    std::cout << "Kernel No is out of range so Cubic Spline is used" << std::endl;
    	Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

        if ((q>=0.0)&&(q<1)) return C*(-3.0+(9.0/2.0)*q)+(Dim-1)/q*C*(-3.0*q+(9.0/4.0)*q*q);
        else if (q<=2)       return C*((3.0/2.0)*(2.0-q))+(Dim-1)/q*C*((-3.0/4.0)*(2.0-q)*(2.0-q));
        else                 return 0.0;

       break;
    }
}

inline double Interaction::SecDerivativeKernel(double r, double h)
{
	double C;
	double q = r/h;
    if (Dim<=1 || Dim>3)
    {
    	std::cout << "Please correct dimension for kernel and run again, otherwise 3D is used" << std::endl;
    	Dim = 3;
    }
	switch (KernelType)
    {case 0:
    	// Qubic Spline
    	Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

        if ((q>=0.0)&&(q<1)) return C*(-3.0+(9.0/2.0)*q);
        else if (q<=2)       return C*((3.0/2.0)*(2.0-q));
        else                 return 0.0;

		break;
    case 1:
    	// Quadratic
    	Dim ==2 ? C = 2.0/(h*h*h*h*M_PI) : C = 5.0/(4.0*h*h*h*h*h*M_PI);

    	if (q<=2) 			 return C*(-3.0/4.0);
		else                 return 0.0;

    	break;
    case 2:
    	// Quintic
    	Dim ==2 ? C = 7.0/(4.0*h*h*h*h*M_PI) : C = 7.0/(8.0*h*h*h*h*h*M_PI);

    	if (q<=2) 			 return -3.0*C*pow((1.0-q/2.0),2)*(1.0+3.0*q-3.0*q*q);
		else                 return 0.0;

    	break;
    case 3:
    	// Gaussian with compact support
    	Dim ==2 ? C = 1.0/(h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*pow(M_PI,(3.0/2.0)));

    	if (q<=2) 			 return C*2.0*(2.0*q*q-1.0)*exp(-q*q);
		else                 return 0.0;

    	break;
   default:
		// Qubic Spline
	    std::cout << "Kernel No is out of range so Cubic Spline is used" << std::endl;
    	Dim ==2 ? C = 10.0/(7.0*h*h*h*h*M_PI) : C = 1.0/(h*h*h*h*h*M_PI);

        if ((q>=0.0)&&(q<1)) return C*(-3.0+(9.0/2.0)*q);
        else if (q<=2)       return C*((3.0/2.0)*(2.0-q));
        else                 return 0.0;

       break;
    }
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
