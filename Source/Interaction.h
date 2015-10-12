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

// Std Lib
#include <stdio.h>			/// for NULL
#include <algorithm>		/// for min,max

#include "Particle.h"
#include "Functions.h"
#include "Domain.h"

namespace SPH {

inline void Domain::CalcForceFF(Particle * P1, Particle * P2)
{
	double h = (P1->h+P2->h)/2;
	double di,dj;
    double mi = P1->Mass;
    double mj = P2->Mass;
    double Mui = P1->Mu;
    double Muj = P2->Mu;

    if (!P1->IsFree) di = DensitySolid(PresEq, Cs, P0,P1->Pressure, P2->RefDensity); else di = P1->Density;
    if (!P2->IsFree) dj = DensitySolid(PresEq, Cs, P0,P2->Pressure, P1->RefDensity); else dj = P2->Density;

    Vec3_t vij = P1->v - P2->v;
    Vec3_t xij = P1->x - P2->x;

    // Correction of xij for Periodic BC
    if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

    double rij	= norm(xij);
    double GK	= GradKernel(Dimension, KernelType, rij, h) / rij;
    double K	= Kernel(Dimension, KernelType, rij, h);

    //Artificial Viscosity
    double PIij = 0.0;
    if (Alpha!=0.0 || Beta!=0.0)
    {
    	double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    	double Cij = 0.5*(SoundSpeed(PresEq, Cs, di, P1->RefDensity)+SoundSpeed(PresEq, Cs, dj, P2->RefDensity));
    	if (dot(vij,xij)<0) PIij = (-Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj));                          ///<(2.74) Li, Liu Book
    }

    //Tensile Instability
    double TIij = 0.0;
    if (TI > 0.0)
    {
    	double Ri,Rj;
    	Ri = 0.0;
    	Rj = 0.0;
    	if (P1->Pressure < 0.0) Ri = -P1->Pressure/(di*di);
    	if (P2->Pressure < 0.0) Rj = -P2->Pressure/(dj*dj);
        TIij = TI*(Ri + Rj)*pow((K/Kernel(Dimension, KernelType, InitialDist, h)),4);
    }

	//Real Viscosity
    Mat3_t StrainRate;
    Vec3_t VI = 0.0;
    if (NoSlip || P1->IsFree*P2->IsFree)
    {
		Vec3_t vab;
		if (P1->IsFree*P2->IsFree)
		{
			vab = vij;
		}
		else
		{
			// No-Slip velocity correction
			if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->vb); else vab = (2.0*P1->v-P1->vb) - P2->v;
		}

		StrainRate = 2.0*vab(0)*xij(0)           , vab(0)*xij(1)+vab(1)*xij(0) , vab(0)*xij(2)+vab(2)*xij(0) ,
					 vab(0)*xij(1)+vab(1)*xij(0) , 2.0*vab(1)*xij(1)           , vab(1)*xij(2)+vab(2)*xij(1) ,
					 vab(0)*xij(2)+vab(2)*xij(0) , vab(1)*xij(2)+vab(2)*xij(1) , 2.0*vab(2)*xij(2)           ;
		StrainRate = -GK * StrainRate;

		if (!P1->IsFree) Mui = Muj;
		if (!P2->IsFree) Muj = Mui;

		if (VisEq==0) VI = (Mui+Muj)/(di*dj)*GK*vab;																		//Morris et al 1997
		if (VisEq==1) VI = 4.0*(Mui+Muj)/((di+dj)*(di+dj))*GK*vab;																//Shao et al 2003
		if (VisEq==2) VI = -(Mui+Muj)/2.0/(di*dj)*LaplaceKernel(Dimension, KernelType, rij, h)*vab;								//Real Viscosity (considering incompressible fluid)
		if (VisEq==3) VI = -(Mui+Muj)/2.0/(di*dj)*( LaplaceKernel(Dimension, KernelType, rij, h)*vab +
				1.0/3.0*(GK*vij + dot(vij,xij) * xij / (rij*rij) * (-GK+SecDerivativeKernel(Dimension, KernelType, rij, h) ) ) );	//Takeda et al 1994 (Real viscosity considering 1/3Mu for compressibility as per Navier Stokes but ignore volumetric viscosity)
		if ((VisEq<0 || VisEq>3))
		{
			std::cout << "Viscosity Equation No is out of range. Please correct it and run again" << std::endl;
			std::cout << "0 => Morris et al 1997" << std::endl;
			std::cout << "1 => Shao et al 2003" << std::endl;
			std::cout << "2 => Real viscosity for incompressible fluids" << std::endl;
			std::cout << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)" << std::endl;
			abort();
		}
    }

    // XSPH Monaghan
    if (XSPH != 0.0)
    {
        omp_set_lock(&P1->my_lock);
        P1->VXSPH		+= XSPH*mj/(0.5*(di+dj))*K*-vij;
        omp_unset_lock(&P1->my_lock);

        omp_set_lock(&P2->my_lock);
		P2->VXSPH		+= XSPH*mi/(0.5*(di+dj))*K*vij;
		omp_unset_lock(&P2->my_lock);
    }


    omp_set_lock(&P1->my_lock);
    P1->a   += -mj * ( P1->Pressure/(di*di) + P2->Pressure/(dj*dj) + PIij + TIij ) * GK*xij + mj*VI;
//    P1->Vis +=  mj * VI;
    if (P1->IsFree) P1->StrainRate = P1->StrainRate + mj/dj*StrainRate; else P1->StrainRate = 0.0;
    if (P1->ct==30 && Shepard)
    {
    	P1->SumDen += mj*    K;
    	P1->ZWab   += mj/dj* K;
    }
    P1->dDensity += di * (mj/dj) * dot( vij , GK*xij );
    omp_unset_lock(&P1->my_lock);


    omp_set_lock(&P2->my_lock);
    P2->a   -= -mi * ( P1->Pressure/(di*di) + P2->Pressure/(dj*dj) + PIij + TIij ) * GK*xij + mi*VI;
//    P2->Vis -=  mi * VI;
    if (P2->IsFree) P2->StrainRate = P2->StrainRate + mi/di*StrainRate; else P2->StrainRate = 0.0;
    if (P2->ct==30 && Shepard)
    {
    	P2->SumDen += mi*    K;
    	P2->ZWab   += mi/di* K;
    }
    P2->dDensity += dj * (mi/di) * dot( -vij , -GK*xij );
    omp_unset_lock(&P2->my_lock);
}

inline void Domain::CalcForceSS(Particle * P1, Particle * P2)
{
	double h = (P1->h+P2->h)/2;
	double di,dj;
    double mi = P1->Mass;
    double mj = P2->Mass;

    if (!P1->IsFree) di = DensitySolid(PresEq, Cs, P0,P1->Pressure, P2->RefDensity); else di = P1->Density;
    if (!P2->IsFree) dj = DensitySolid(PresEq, Cs, P0,P2->Pressure, P1->RefDensity); else dj = P2->Density;

    Vec3_t vij = P1->v - P2->v;
    Vec3_t xij = P1->x - P2->x;

    // Correction of xij for Periodic BC
    if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

    double rij	= norm(xij);
    double GK	= GradKernel(Dimension, KernelType, rij, h) / rij;
    double K	= Kernel(Dimension, KernelType, rij, h);

    //Artificial Viscosity
    Mat3_t PIij;
    set_to_zero(PIij);
    if (Alpha!=0.0 || Beta!=0.0)
    {
    	double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);                                                ///<(2.75) Li, Liu Book
    	double Cij = 0.5*(SoundSpeed(PresEq, Cs, di, P1->RefDensity)+SoundSpeed(PresEq, Cs, dj, P2->RefDensity));
    	if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;                          ///<(2.74) Li, Liu Book
    }

    //Tensile Instability
    Mat3_t TIij;
    set_to_zero(TIij);

    if (TI > 0.0)
    {
        Mat3_t Ri, Rj;
        set_to_zero(Ri);
        set_to_zero(Rj);

        // XY plane must be used, It is very slow in 3D
        if (Dimension == 2)
        {
			double teta, Sigmaxx, Sigmayy, C, S;

			if ((P1->Sigma(0,0)-P1->Sigma(1,1))!=0.0) teta = 0.5*atan(2.0*P1->Sigma(0,1)/(P1->Sigma(0,0)-P1->Sigma(1,1))); else teta = M_PI/4.0;
			C = cos(teta);
			S = sin(teta);
			Sigmaxx = C*C*P1->Sigma(0,0) + 2.0*C*S*P1->Sigma(0,1) + S*S*P1->Sigma(1,1);
			Sigmayy = S*S*P1->Sigma(0,0) - 2.0*C*S*P1->Sigma(0,1) + C*C*P1->Sigma(1,1);
			if (Sigmaxx>0) Sigmaxx = -TI * Sigmaxx/(di*di); else Sigmaxx = 0.0;
			if (Sigmayy>0) Sigmayy = -TI * Sigmayy/(di*di); else Sigmayy = 0.0;
			Ri(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
			Ri(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
			Ri(0,1) = Ri(1,0) = S*C*(Sigmaxx-Sigmayy);

			if ((P2->Sigma(0,0)-P2->Sigma(1,1))!=0.0) teta = 0.5*atan(2.0*P2->Sigma(0,1)/(P2->Sigma(0,0)-P2->Sigma(1,1))); else teta = M_PI/4.0;
			C = cos(teta);
			S = sin(teta);
			Sigmaxx = C*C*P2->Sigma(0,0) + 2.0*C*S*P2->Sigma(0,1) + S*S*P2->Sigma(1,1);
			Sigmayy = S*S*P2->Sigma(0,0) - 2.0*C*S*P2->Sigma(0,1) + C*C*P2->Sigma(1,1);
			if (Sigmaxx>0) Sigmaxx = -TI * Sigmaxx/(dj*dj); else Sigmaxx = 0.0;
			if (Sigmayy>0) Sigmayy = -TI * Sigmayy/(dj*dj); else Sigmayy = 0.0;
			Rj(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
			Rj(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
			Rj(0,1) = Rj(1,0) = S*C*(Sigmaxx-Sigmayy);
        }
        else
        {
        	Mat3_t Vec,Val,VecT,temp;

        	Rotation(P1->Sigma,Vec,VecT,Val);
			if (Val(0,0)>0) Val(0,0) = -TI * Val(0,0)/(di*di); else Val(0,0) = 0.0;
			if (Val(1,1)>0) Val(1,1) = -TI * Val(1,1)/(di*di); else Val(1,1) = 0.0;
			if (Val(2,2)>0) Val(2,2) = -TI * Val(2,2)/(di*di); else Val(2,2) = 0.0;
			Mult(Vec,Val,temp);
			Mult(temp,VecT,Ri);

			Rotation(P2->Sigma,Vec,VecT,Val);
			if (Val(0,0)>0) Val(0,0) = -TI * Val(0,0)/(di*di); else Val(0,0) = 0.0;
			if (Val(1,1)>0) Val(1,1) = -TI * Val(1,1)/(di*di); else Val(1,1) = 0.0;
			if (Val(2,2)>0) Val(2,2) = -TI * Val(2,2)/(di*di); else Val(2,2) = 0.0;
			Mult(Vec,Val,temp);
			Mult(temp,VecT,Rj);

        }

        TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),4)*(Ri+Rj);
    }

    Mat3_t StrainRate,Rotation;
    set_to_zero(StrainRate);
    set_to_zero(Rotation);

    Vec3_t vab;
	if (P1->IsFree*P2->IsFree)
	{
		vab = vij;
	}
	else
	{
		if (NoSlip)
		{
			// No-Slip velocity correction
			if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->vb); else vab = (2.0*P1->v-P1->vb) - P2->v;
		}
		else vab = vij;
	}

    // Strain Rate Pattern
//	StrainRate = 2.0*vab(0)*xij(0)           , vab(0)*xij(1)+vab(1)*xij(0) , vab(0)*xij(2)+vab(2)*xij(0) ,
//				 vab(0)*xij(1)+vab(1)*xij(0) , 2.0*vab(1)*xij(1)           , vab(1)*xij(2)+vab(2)*xij(1) ,
//				 vab(0)*xij(2)+vab(2)*xij(0) , vab(1)*xij(2)+vab(2)*xij(1) , 2.0*vab(2)*xij(2)           ;
	StrainRate(0,0) = 2.0*vab(0)*xij(0);
	StrainRate(0,1) = vab(0)*xij(1)+vab(1)*xij(0);
	StrainRate(0,2) = vab(0)*xij(2)+vab(2)*xij(0);
	StrainRate(1,0) = StrainRate(0,1);
	StrainRate(1,1) = 2.0*vab(1)*xij(1);
	StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
	StrainRate(2,0) = StrainRate(0,2);
	StrainRate(2,1) = StrainRate(1,2);
	StrainRate(2,2) = 2.0*vab(2)*xij(2);
	StrainRate = -0.5 * GK * StrainRate;

	// Rotation Pattern
//	Rotation   =   0.0                        , vab(0)*xij(1)-vab(1)*xij(0) , vab(0)*xij(2)-vab(2)*xij(0) ,
//				  vab(1)*xij(0)-vab(0)*xij(1) , 0.0                         , vab(1)*xij(2)-vab(2)*xij(1) ,
//				  vab(2)*xij(0)-vab(0)*xij(2) , vab(2)*xij(1)-vab(1)*xij(2) , 0.0                         ;
	Rotation(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
	Rotation(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
	Rotation(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
	Rotation(1,0) = -Rotation(0,1);
	Rotation(2,0) = -Rotation(0,2);
	Rotation(2,1) = -Rotation(1,2);
	Rotation = -0.5 * GK * Rotation;

    // XSPH Monaghan
    if (XSPH != 0.0)
    {
        omp_set_lock(&P1->my_lock);
        P1->VXSPH		+= XSPH*mj/(0.5*(di+dj))*K*-vij;
        omp_unset_lock(&P1->my_lock);

        omp_set_lock(&P2->my_lock);
		P2->VXSPH		+= XSPH*mi/(0.5*(di+dj))*K*vij;
		omp_unset_lock(&P2->my_lock);
    }

    Vec3_t temp;
    Mult( GK*xij , mj * ( 1.0/(di*di)*P1->Sigma + 1.0/(dj*dj)*P2->Sigma + PIij + TIij ) , temp);
    omp_set_lock(&P1->my_lock);
    P1->a   += temp;
    if (P1->IsFree) P1->StrainRate = P1->StrainRate + mj/dj*StrainRate;
    if (P1->IsFree) P1->Rotation = P1->Rotation + mj/dj*Rotation;
    if (P1->ct==30 && Shepard)
    {
    	P1->SumDen += mj*    K;
    	P1->ZWab   += mj/dj* K;
    }
    P1->dDensity += di * (mj/dj) * dot( vij , GK*xij );
    omp_unset_lock(&P1->my_lock);

    Mult( GK*xij , mi * ( 1.0/(di*di)*P1->Sigma + 1.0/(dj*dj)*P2->Sigma + PIij + TIij ) ,temp);
    omp_set_lock(&P2->my_lock);
    P2->a   -= temp;
    if (P2->IsFree) P2->StrainRate = P2->StrainRate + mi/di*StrainRate;
    if (P2->IsFree) P2->Rotation = P2->Rotation + mi/di*Rotation;
    if (P2->ct==30 && Shepard)
    {
    	P2->SumDen += mi*    K;
    	P2->ZWab   += mi/di* K;
    }
    P2->dDensity += dj * (mi/di) * dot( -vij , -GK*xij );
    omp_unset_lock(&P2->my_lock);
}

}; // namespace SPH
