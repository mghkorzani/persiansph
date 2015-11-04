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

inline void Domain::CalcForce11(Particle * P1, Particle * P2)
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
        TIij = TI*(Ri + Rj)*pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn);
    }

	//Real Viscosity
    Mat3_t StrainRate;
    Vec3_t VI = 0.0;
    if ((P1->NoSlip || P2->NoSlip) || P1->IsFree*P2->IsFree)
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

inline void Domain::CalcForce22(Particle * P1, Particle * P2)
{
	double h = (P1->h+P2->h)/2;
	double di,dj;
    double mi = P1->Mass;
    double mj = P2->Mass;
    Mat3_t Sigmaj,Sigmai;

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

	Sigmai = P1->Sigma;
	Sigmaj = P2->Sigma;

    //Tensile Instability
    Mat3_t TIij;
    set_to_zero(TIij);

    if (TI > 0.0)
    {
    	if (P1->IsFree*P2->IsFree)
             TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn)*(P1->TIR+P2->TIR);
     	else
    	{
    		if (P1->IsFree) TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn)*(2.0*P1->TIR);
    		if (P2->IsFree) TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn)*(2.0*P2->TIR);
    	}
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
		if (P1->NoSlip || P2->NoSlip)
		{
			// No-Slip velocity correction
			if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->vb); else vab = (2.0*P1->v-P1->vb) - P2->v;
		}
		if (!(P1->NoSlip || P2->NoSlip))
		{
			if (P1->IsFree)	vab = P1->v - P2->vb; else vab = P1->vb - P2->v;

		}

	}

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
    Mult( GK*xij , mj * ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
    omp_set_lock(&P1->my_lock);
		P1->a		+= temp;
		P1->dDensity+= di * (mj/dj) * dot( vij , GK*xij );
		if (P1->IsFree)
		{
			P1->StrainRate = P1->StrainRate + mj/dj*StrainRate;
			P1->Rotation = P1->Rotation + mj/dj*Rotation;
		}
    omp_unset_lock(&P1->my_lock);

    Mult( GK*xij , mi * ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) ,temp);
    omp_set_lock(&P2->my_lock);
		P2->a		-= temp;
		P2->dDensity+= dj * (mi/di) * dot( -vij , -GK*xij );
		if (P2->IsFree)
		{
			P2->StrainRate = P2->StrainRate + mi/di*StrainRate;
			P2->Rotation = P2->Rotation + mi/di*Rotation;
		}
    omp_unset_lock(&P2->my_lock);
}

inline void Domain::CalcForce33(Particle * P1, Particle * P2)
{
	double h = (P1->h+P2->h)/2;
	double di,dj;
    double mi = P1->Mass;
    double mj = P2->Mass;
    Mat3_t Sigmaj,Sigmai;

    if (!P1->IsFree) di = P2->Density; else di = P1->Density;
    if (!P2->IsFree) dj = P1->Density; else dj = P2->Density;

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

	Sigmai = P1->Sigma;
	Sigmaj = P2->Sigma;

    //Tensile Instability
    Mat3_t TIij;
    set_to_zero(TIij);

    if (TI > 0.0)
    {
    	if (P1->IsFree*P2->IsFree)
             TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn)*(P1->TIR+P2->TIR);
     	else
    	{
    		if (P1->IsFree) TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn)*(2.0*P1->TIR);
    		if (P2->IsFree) TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),TIn)*(2.0*P2->TIR);
    	}
    }

    Mat3_t StrainRate,Rotation;
    set_to_zero(StrainRate);
    set_to_zero(Rotation);

    Vec3_t vab = 0.0;
	if (P1->IsFree*P2->IsFree)
	{
		vab = vij;
	}
	else
	{
		if (P1->NoSlip || P2->NoSlip)
		{
			// No-Slip velocity correction
			if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->vb); else vab = (2.0*P1->v-P1->vb) - P2->v;
		}
		if (!(P1->NoSlip || P2->NoSlip))
		{
			if (P1->IsFree)
			{
				if (P2->FreeSlip(0) > 0.0)
					vab = 2.0*P1->v(0),0.0,0.0;
				else if (P2->FreeSlip(1) > 0.0)
					vab = 0.0,2.0*P1->v(1),0.0;
				else if (P2->FreeSlip(2) > 0.0)
					vab = 0.0,0.0,2.0*P1->v(2);
				else vab = P1->v;
			}
			else
			{
				if (P1->FreeSlip(0) > 0.0)
					vab = -2.0*P2->v(0),0.0,0.0;
				else if (P1->FreeSlip(1) > 0.0)
					vab = 0.0,-2.0*P2->v(1),0.0;
				else if (P1->FreeSlip(2) > 0.0)
					vab = 0.0,0.0,-2.0*P2->v(2);
				else vab = -P2->v;
			}
		}

	}

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
        P1->VXSPH += XSPH*mj/(0.5*(di+dj))*K*-vij;
        omp_unset_lock(&P1->my_lock);

        omp_set_lock(&P2->my_lock);
		P2->VXSPH += XSPH*mi/(0.5*(di+dj))*K*vij;
		omp_unset_lock(&P2->my_lock);
    }

    Vec3_t temp;
    Mult( GK*xij , mj * ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
    omp_set_lock(&P1->my_lock);
		P1->a += temp;
		if (P1->IsFree)
		{
			P1->StrainRate	= P1->StrainRate + mj/dj*StrainRate;
			P1->Rotation	= P1->Rotation + mj/dj*Rotation;
		}
    omp_unset_lock(&P1->my_lock);

    Mult( GK*xij , mi * ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) ,temp);
    omp_set_lock(&P2->my_lock);
		P2->a -= temp;
		if (P2->IsFree)
		{
			P2->StrainRate	= P2->StrainRate + mi/di*StrainRate;
			P2->Rotation	= P2->Rotation + mi/di*Rotation;
		}
    omp_unset_lock(&P2->my_lock);
}

}; // namespace SPH
