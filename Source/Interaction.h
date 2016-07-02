/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

// Std Lib
#include <stdio.h>		/// for NULL
#include <algorithm>		/// for min,max

#include "Particle.h"
#include "Functions.h"
#include "Domain.h"

namespace SPH {

inline void Domain::CalcForce11(Particle * P1, Particle * P2)
{
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	// Correction of xij for Periodic BC
	if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;
		Vec3_t vij	= P1->v - P2->v;


		if (!P1->IsFree)
		{
			di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
			mi = P2->Mass;
		}
		else
		{
			di = P1->Density;
			mi = P1->Mass;
		}

		if (!P2->IsFree)
		{
			dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
			mj = P1->Mass;
		}
		else
		{
			dj = P2->Density;
			mj = P2->Mass;
		}

		double GK	= GradKernel(Dimension, KernelType, rij, h);
		double K	= Kernel(Dimension, KernelType, rij, h);

		// Artificial Viscosity
		double PIij = 0.0;
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);						///<(2.75) Li, Liu Book
			if (dot(vij,xij)<0) PIij = (-Alpha*0.5*(Ci+Cj)*MUij+Beta*MUij*MUij)/(0.5*(di+dj));		///<(2.74) Li, Liu Book
		}

		// Tensile Instability
		double TIij = 0.0;
		if (P1->TI > 0.0 || P2->TI > 0.0)
		{
			double Ri,Rj;
			Ri = 0.0;
			Rj = 0.0;
			if ((P1->Pressure+P2->Pressure) < 0.0)
			{
				if (P1->Pressure < 0.0) Ri = -P1->Pressure/(di*di);
				if (P2->Pressure < 0.0) Rj = -P2->Pressure/(dj*dj);
			}
			TIij = (P1->TI*Ri + P2->TI*Rj)*pow((K/Kernel(Dimension, KernelType, InitialDist, h)),(P1->TIn+P2->TIn)/2.0);
		}

		// Real Viscosity
		Mat3_t StrainRate;
		set_to_zero(StrainRate);
		Vec3_t VI = 0.0;
			Vec3_t vab=0.0;
		if ((P1->NoSlip || P2->NoSlip) || (P1->IsFree*P2->IsFree))
		{
			double Mu = 0.0;
			if (P1->IsFree*P2->IsFree)
			{
				vab	= vij;
				Mu	= 2.0*P1->Mu*P2->Mu/(P1->Mu+P2->Mu);
			}
			else
			{
				// No-Slip velocity correction
				if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->vb); else vab = (2.0*P1->v-P1->vb) - P2->v;
				if (!P1->IsFree) Mu = P2->Mu;
				if (!P2->IsFree) Mu = P1->Mu;
			}

			StrainRate =	2.0*vab(0)*xij(0)            , vab(0)*xij(1)+vab(1)*xij(0) , vab(0)*xij(2)+vab(2)*xij(0) ,
					vab(0)*xij(1)+vab(1)*xij(0)  , 2.0*vab(1)*xij(1)           , vab(1)*xij(2)+vab(2)*xij(1) ,
					 vab(0)*xij(2)+vab(2)*xij(0) , vab(1)*xij(2)+vab(2)*xij(1) , 2.0*vab(2)*xij(2)           ;
			StrainRate = -GK * StrainRate;

			if (VisEq==0) VI =  2.0*Mu / (di*dj)          * GK*vab;								//Morris et al 1997
			if (VisEq==1) VI =  8.0*Mu / ((di+dj)*(di+dj))* GK*vab;								//Shao et al 2003
			if (VisEq==2) VI = -Mu     / (di*dj)          * LaplaceKernel(Dimension, KernelType, rij, h)*vab;		//Real Viscosity (considering incompressible fluid)
			if (VisEq==3) VI = -Mu     / (di*dj)          * ( LaplaceKernel(Dimension, KernelType, rij, h)*vab +
						1.0/3.0*(GK*vij + dot(vij,xij) * xij / (rij*rij) * 
						(-GK+SecDerivativeKernel(Dimension, KernelType, rij, h) ) ) );				//Takeda et al 1994
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

		// Calculating the forces for the particle 1 & 2				
		Vec3_t temp	= 0.0;
		double temp1	= 0.0;
//		temp		= ( P1->Pressure/(di*di) + P2->Pressure/(dj*dj) + PIij + TIij ) * GK*xij + VI;
		temp		= ( (P1->Pressure + P2->Pressure)/(di*dj)       + PIij + TIij ) * GK*xij + VI;
		if (Dimension == 2) temp(2) = 0.0;
		temp1		= dot( vij , GK*xij );

		omp_set_lock(&P1->my_lock);
			P1->a		+= -mj * temp;
			P1->dDensity	+=  mj * (di/dj) * temp1;
			if (P1->IsFree) P1->StrainRate = P1->StrainRate + mj/dj*StrainRate; else P1->StrainRate = 0.0;
			if (P1->IsFree && SWIType ==1) P1->S = P1->S + mj/dj*vab(0)*xij(1)*-GK; else P1->S = 0.0;
			if (P1->Shepard)
				if (P1->ShepardCounter == P1->ShepardStep && (P1->IsFree*P2->IsFree))
				{
					P1->SumDen += mj*    K;
					P1->ZWab   += mj/dj* K;
					P1->ShepardNeighbourNo++;
				}
		omp_unset_lock(&P1->my_lock);


		omp_set_lock(&P2->my_lock);
			P2->a		-= -mi * temp;
			P2->dDensity	+=  mi * (dj/di) * temp1;
			if (P2->IsFree) P2->StrainRate = P2->StrainRate + mi/di*StrainRate; else P2->StrainRate = 0.0;
			if (P2->IsFree && SWIType ==1) P2->S = P2->S + mi/di*vab(0)*xij(1)*-GK; else P2->StrainRate = 0.0;
			if (P2->Shepard)
				if (P2->ShepardCounter == P2->ShepardStep && (P1->IsFree*P2->IsFree))
				{
					P2->SumDen += mi*    K;
					P2->ZWab   += mi/di* K;
					P2->ShepardNeighbourNo++;
				}
		omp_unset_lock(&P2->my_lock);
    }
}

inline void Domain::CalcForce2233(Particle * P1, Particle * P2)
{
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	// Correction of xij for Periodic BC
	if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;

		if (P1->Material*P2->Material == 9)
		{
			if (!P1->IsFree) {di = P2->Density;	mi = P2->Mass;} else {di = P1->Density;	mi = P1->Mass;}
			if (!P2->IsFree) {dj = P1->Density;	mj = P1->Mass;} else {dj = P2->Density;	mj = P2->Mass;}
		}
		else
		{
			if (!P1->IsFree)
			{
				di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
				mi = P2->Mass;
			}
			else
			{
				di = P1->Density;
				mi = P1->Mass;
			}
			if (!P2->IsFree)
			{
				dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
				mj = P1->Mass;
			}
			else
			{
				dj = P2->Density;
				mj = P2->Mass;
			}
		}

		Vec3_t vij	= P1->v - P2->v;
		double GK	= GradKernel(Dimension, KernelType, rij, h);
		double K	= Kernel(Dimension, KernelType, rij, h);

		// Artificial Viscosity
		Mat3_t PIij;
		set_to_zero(PIij);
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
			double Cij;
			if (P1->Material*P2->Material == 9)
				Cij = 0.5*(P1->Cs+P2->Cs);
			else
			{
				double Ci,Cj;
				if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
				if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
				Cij = 0.5*(Ci+Cj);
			}
			if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
		}

		Mat3_t Sigmaj,Sigmai;
		set_to_zero(Sigmaj);
		set_to_zero(Sigmai);
		Sigmai = P1->Sigma;
		Sigmaj = P2->Sigma;

//		if (P1->IsFree) Sigmai = P1->Sigma; else  Sigmai = P2->Sigma;
//		if (P2->IsFree) Sigmaj = P2->Sigma; else  Sigmaj = P1->Sigma;

		// Tensile Instability
		Mat3_t TIij;
		set_to_zero(TIij);
		if (P1->TI > 0.0 || P2->TI > 0.0) TIij = pow((K/Kernel(Dimension, KernelType, InitialDist, h)),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);

		// NoSlip BC velocity correction
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
			// Please check
			if (!(P1->NoSlip || P2->NoSlip))
			{
				if (P1->IsFree) vab = P1->v - P2->vb; else vab = P1->vb - P2->v;
//				if (P1->IsFree) vab(0) = P1->v(0) + P2->vb(0); else vab(0) = -P1->vb(0) - P2->v(0);
			}
		}

		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);

		// Calculation strain rate tensor
		StrainRate(0,0) = 2.0*vab(0)*xij(0);
		StrainRate(0,1) = vab(0)*xij(1)+vab(1)*xij(0);
		StrainRate(0,2) = vab(0)*xij(2)+vab(2)*xij(0);
		StrainRate(1,0) = StrainRate(0,1);
		StrainRate(1,1) = 2.0*vab(1)*xij(1);
		StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
		StrainRate(2,0) = StrainRate(0,2);
		StrainRate(2,1) = StrainRate(1,2);
		StrainRate(2,2) = 2.0*vab(2)*xij(2);
		StrainRate	= -0.5 * GK * StrainRate;

		// Calculation rotation rate tensor
		RotationRate(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
		RotationRate(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
		RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
		RotationRate(1,0) = -RotationRate(0,1);
		RotationRate(2,0) = -RotationRate(0,2);
		RotationRate(2,1) = -RotationRate(1,2);
		RotationRate	  = -0.5 * GK * RotationRate;

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

		// Calculating the forces for the particle 1 & 2				
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
//		Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
		Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj) + PIij + TIij ) , temp);
		if (Dimension == 2) temp(2) = 0.0;
		temp1 = dot( vij , GK*xij );

		// Locking the particle 1 for updating the properties
		omp_set_lock(&P1->my_lock);
			P1->a		+= mj * temp;
			P1->dDensity	+= mj * (di/dj) * temp1;
			if (P1->Shepard)
				if (P1->ShepardCounter == P1->ShepardStep && (P1->IsFree*P2->IsFree))
				{
					P1->SumDen += mj*    K;
					P1->ZWab   += mj/dj* K;
					P1->ShepardNeighbourNo++;
				}
			if (P1->IsFree && SWIType ==1) P1->S = P1->S + mj/dj*vab(0)*xij(1)*-GK; else P1->S = 0.0;
			if (P1->IsFree)
			{
				P1->StrainRate	= P1->StrainRate + mj/dj*StrainRate;
				P1->RotationRate= P1->RotationRate + mj/dj*RotationRate;
			}
		omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		omp_set_lock(&P2->my_lock);
			P2->a		-= mi * temp;
			P2->dDensity	+= mi * (dj/di) * temp1;
			if (P2->Shepard)
				if (P2->ShepardCounter == P2->ShepardStep && (P1->IsFree*P2->IsFree))
				{
					P2->SumDen += mi*    K;
					P2->ZWab   += mi/di* K;
					P2->ShepardNeighbourNo++;
				}
			if (P2->IsFree  && SWIType ==1) P2->S = P2->S + mi/di*vab(0)*xij(1)*-GK; else P2->StrainRate = 0.0;
			if (P2->IsFree)
			{
				P2->StrainRate	= P2->StrainRate + mi/di*StrainRate;
				P2->RotationRate= P2->RotationRate + mi/di*RotationRate;
			}
		omp_unset_lock(&P2->my_lock);
	}
}

inline void Domain::CalcForce13(Particle * P1, Particle * P2)
{
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	// Correction of xij for Periodic BC
	if (DomSize(0)>0.0) {if (xij(0)>2*Cellfac*h || xij(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? xij(0) -= DomSize(0) : xij(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (xij(1)>2*Cellfac*h || xij(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? xij(1) -= DomSize(1) : xij(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (xij(2)>2*Cellfac*h || xij(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? xij(2) -= DomSize(2) : xij(2) += DomSize(2);}}

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double K	= Kernel(Dimension, KernelType, rij, h);
		double di	= P1->Density;
		double dj	= P2->Density;
		double SF1=0.0,SF2=0.0;
		Vec3_t SFt=0.0,v=0.0;
		switch(SWIType)
		{
			case 0:
				if (P1->Material == 3 )
				{
					v = P2->v-P1->v;
					Seepage(P1->SeepageType, P1->k, P1->k2, P2->Mu, P2->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) / (di*dj)*K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a += P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a -= P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				else
				{
					v = P1->v-P2->v;
					Seepage(P2->SeepageType, P2->k, P2->k2, P1->Mu, P1->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) / (di*dj)*K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a -= P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a += P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				break;
			case 1:
				if (P1->Material == 3 )
				{
					v = P2->v-P1->v;
					if (P1->ZWab<0.6)
					{
						double Cd = 24.0*(P2->Mu/P2->RefDensity)/(P1->d*norm(v)+0.01*h*h) + 2.0;
						SFt = (3.0/(4.0*P1->d)*P2->RefDensity*(1.0-P1->n)*Cd*norm(v)*v) / (di*dj)*K;
						SFt(1) += (P2->RefDensity*(1.0-P1->n)*norm(v)*(P2->S-P1->S)) / (di*dj)*K;
					}
					else
					{
						Seepage(P1->SeepageType, P1->k, P1->k2, P2->Mu, P2->RefDensity, SF1, SF2);
						SFt = (SF1*v + SF2*norm(v)*v) / (di*dj)*K;
					}			
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a += P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a -= P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				else
				{
					v = P1->v-P2->v;
					if (P2->ZWab<0.6)
					{
						double Cd = 24.0*(P1->Mu/P1->RefDensity)/(P2->d*norm(v)+0.01*h*h) + 2.0;
						SFt = (3.0/(4.0*P2->d)*P1->RefDensity*(1.0-P2->n)*Cd*norm(v)*v) / (di*dj)*K;
						SFt(1) += (P1->RefDensity*(1.0-P2->n)*norm(v)*(P1->S-P2->S))/ (di*dj)*K;;
					}
					else
					{
						Seepage(P2->SeepageType, P2->k, P2->k2, P1->Mu, P1->RefDensity, SF1, SF2);
						SFt = (SF1*v + SF2*norm(v)*v) / (di*dj)*K;
					}			
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a -= P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a += P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				break;
			case 2:
				if (P1->Material == 3 )
				{
					double GK	= GradKernel(Dimension, KernelType, rij, h);
					v = P2->v-P1->v;
					Seepage(P1->SeepageType, P1->k, P1->k2, P2->Mu, P2->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) / (di*dj)*K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a += P2->Mass*SFt - P2->Mass*P2->Pressure/(di*dj)*GK*xij;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a -= P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				else
				{
					double GK	= GradKernel(Dimension, KernelType, rij, h);
					v = P1->v-P2->v;
					Seepage(P2->SeepageType, P2->k, P2->k2, P1->Mu, P1->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) / (di*dj)*K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a -= P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a += P1->Mass*SFt + P1->Mass*P1->Pressure/(di*dj)*GK*xij;
					omp_unset_lock(&P2->my_lock);
				}
				break;
			default:
				std::cout << "Soil-Water Interaction Type No is out of range. Please correct it and run again" << std::endl;
				std::cout << "0 => The seepage force + The bouyant unit weight of soil" << std::endl;
				std::cout << "1 => The seepage force + The surface erosion(Lift+Drag)) + The bouyant unit weight of soil" << std::endl;
				std::cout << "2 => The seepage force + The pore water pressure from water particles " << std::endl;
				abort();
				break;
		}
	}
}

}; // namespace SPH
