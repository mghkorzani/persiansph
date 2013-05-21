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

#ifndef MECHSYS_SPH_FUNCTIONS_H
#define MECHSYS_SPH_FUNCTIONS_H

// Std lib
#include <iostream>
#include <cmath>

namespace SPH {

inline double Kernel(double r,double h)
{
    double C = 10.0/(7*h*h*h*M_PI);
    double q = r/h;
    if ((q>=0.0)&&(q<1)) return C*(1-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
    else if (q<=2)       return C*((1.0/4.0)*(2-q)*(2-q)*(2-q));
    else                 return 0.0;
}

inline double GradKernel(double r, double h)
{
    double C = 10.0/(7*h*h*h*M_PI);
    double q = r/h;
    if ((q>=0.0)&&(q<1)) return C*(1-3.0*q+(9.0/4.0)*q*q);
    else if (q<=2)       return C*(-1*(3.0/4.0)*(2-q)*(2-q));
    else                 return 0.0;
}

inline double Pressure(double Density, double Density0)
{
	double v2 = 2*9.81*20;
	return (100*Density0*v2/7)*(pow(Density/Density0,7)-1);
}

inline double SoundSpeed(double Density, double Density0)
{
	double v2 = 2*9.81*20;
	return sqrt(7*(100*Density0*v2/7)*(pow(Density/Density0,6)/Density0));
}

}; // namespace SPH

#endif // MECHSYS_SPH_FUNCTIONS_H
