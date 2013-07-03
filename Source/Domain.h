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

#ifndef MECHSYS_SPH_DOMAIN_H
#define MECHSYS_SPH_DOMAIN_H

// Std Lib
#include <stdio.h>			/// for NULL
#include <algorithm>		/// for min,max


#include <Source/Interaction.h>

// HDF File Output
#include <hdf5.h>
#include <hdf5_hl.h>

#include <mechsys/util/stopwatch.h>

namespace SPH {

class Domain
{
public:

    // Constructor
    Domain();  ///< Constructor with a vector containing the number of divisions per length, and Xmin Xmax defining the limits of the rectangular domain to be plotted

    // Destructor
    ~Domain ();

    // Methods
    void AddBox              (int tag, Vec3_t const & x, size_t nx, size_t ny, size_t nz,
    		                  double R, double Mass, double Density, double h, bool Fixed);                                    ///< Add a box of particles (should specify radius of particles)
    void AddBoxLength        (int tag, Vec3_t const &V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    		                  double Mass, double Density, double h, bool Fixed);                                              ///< Add a box of particles with length (calculate radius of particles)
    void AddRandomBox        (int tag, Vec3_t const &V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    		                  double Mass, double Density, double h, size_t RandomSeed=100);                                   ///< Add box of random positioned particles (calculate radius of particles)
    void StartAcceleration   (Vec3_t const & a = Vec3_t(0.0,0.0,0.0));                                                         ///< Add a fixed acceleration
    void ComputeAcceleration (double dt);                                                                                      ///< Compute the acceleration due to the other particles
    void Move                (double dt);                                                                                      ///< Compute the acceleration due to the other particles
    void ResetInteractions();                                                                                                  ///< Reset the interaction array
    void ResetContacts();                                                                                                      ///< Reset the possible interactions
    void Solve               (double tf, double dt, double dtOut, char const * TheFileKey);                                    ///< The solving function
    void WriteXDMF           (char const * FileKey);                                                                           ///< Save a XDMF file for visualization

    // Data
    Vec3_t                  Gravity;        ///< Gravity acceleration
    Array <Particle*>       Particles;      ///< Array of particles
    Array <Interaction*>    Interactions;   ///< Array of interactions
    Array <Interaction*>    PInteractions;  ///< Array of possible interactions
    double                  Time;           ///< The simulation Time
    size_t                  idx_out;        ///< Index for output purposes
    double 					Dimension;      ///< Dimension of the problem
    double 					Alpha;
    double					Beta;
    double					MaxVel;
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Constructor
inline Domain::Domain ()
{
    Time    = 0.0;
    Gravity = 0.0,0.0,0.0;
    idx_out = 0;

}

inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactions.Size(); ++i) if (Interactions[i]!=NULL) delete Interactions[i];
}


// Methods
inline void Domain::AddBox(int tag, Vec3_t const & V, size_t nx, size_t ny, size_t nz, double R, double Mass, double Density, double h, bool Fixed)
{
    std::cout << "\n--------------Generating particles by AddBox with defined radius-----------------------------------" << std::endl;

    for (size_t i = 0;i<nx;i++)
    for (size_t j = 0;j<ny;j++)
    for (size_t k = 0;k<nz;k++)
    {
        Vec3_t x(2*i*R,2*j*R,2*k*R);
        x +=V;
        if (Mass == 0.0)
        	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),4/3*M_PI*R*R*R*Density,Density,R,h,Fixed));
        else
        	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,R,h,Fixed));
    }
    std::cout << "\n  Total No. of particles   = " << Particles.Size() << std::endl;
}

inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double Mass, double Density, double h, bool Fixed)
{
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length-----------------------------" << std::endl;

	double R = (std::min(Lx/(2*nx),Ly/(2*ny))>0) ? std::min(Lx/(2*nx),Ly/(2*ny)) : std::max(Lx/(2*nx),Ly/(2*ny));
    R = (std::min(R,Lz/(2*nz))>0) ? std::min(R,Lz/(2*nz)) : std::max(R,Lz/(2*nz));

	for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = V(0)+i*Lx/nx;
                double y = V(1)+j*Ly/ny;
                double z = V(2)+k*Lz/nz;
                if (Mass == 0.0)
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),4/3*M_PI*R*R*R*Density,Density,R,h,Fixed));
                else
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),Mass,Density,R,h,Fixed));
            }
        }
    }
    std::cout << "\n  Total No. of particles   = " << Particles.Size() << std::endl;
}

inline void Domain::AddRandomBox(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double Mass, double Density, double h, size_t RandomSeed)
{
    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating random packing of particles by AddRandomBox-------------------------------" << std::endl;

    double R = (std::min(Lx/(2*nx),Ly/(2*ny))>0) ? std::min(Lx/(2*nx),Ly/(2*ny)) : std::max(Lx/(2*nx),Ly/(2*ny));
    R = (std::min(R,Lz/(2*nz))>0) ? std::min(R,Lz/(2*nz)) : std::max(R,Lz/(2*nz));

	std::cout << R << std::endl;

    double qin = 0.95;
    srand(RandomSeed);
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = V(0)+(i+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Lx/nx;
                double y = V(1)+(j+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Ly/ny;
                double z = V(2)+(k+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*Lz/nz;
                if (Mass == 0.0)
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),4/3*M_PI*R*R*R*Density,Density,R,h,false));
                else
                	Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),Mass,Density,R,h,false));
            }
        }
    }
    std::cout << "\n  Total No. of particles   = " << Particles.Size() << std::endl;
}

inline void Domain::StartAcceleration (Vec3_t const & a)
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->a = a;
        Particles[i]->dDensity = 0.0;
    }
}

inline void Domain::ComputeAcceleration (double dt)
{
    for (size_t i=0; i<PInteractions.Size(); i++) PInteractions[i]->CalcForce(dt);
}

inline void Domain::Move (double dt)
{
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Move(dt);
}

inline void Domain::ResetInteractions()
{
    // delete old interactors
    for (size_t i=0; i<Interactions.Size(); ++i)
    {
        if (Interactions[i]!=NULL) delete Interactions[i];
    }

    // new interactors
    Interactions.Resize(0);
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            // if both particles are fixed, don't create any interactor
            if (!Particles[i]->IsFree && !Particles[j]->IsFree) continue;
            else Interactions.Push(new Interaction(Particles[i],Particles[j],Dimension,Alpha,Beta,MaxVel));
        }
    }
}

inline void Domain::ResetContacts()
{
    PInteractions.Resize(0);
    for (size_t i=0; i<Interactions.Size(); i++)
    {
        if(Interactions[i]->UpdateContacts()) PInteractions.Push(Interactions[i]);
    }
}

inline void Domain::WriteXDMF (char const * FileKey)
{

	String fn(FileKey);
    fn.append(".h5");
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    float * Posvec   = new float[3*Particles.Size()];
    float * Velvec   = new float[3*Particles.Size()];
    float * Pressure = new float[  Particles.Size()];
    float * Radius   = new float[  Particles.Size()];
    int   * Tag      = new int  [  Particles.Size()];
    for (size_t i=0;i<Particles.Size();i++)
    {
        Posvec  [3*i  ] = float(Particles[i]->x(0));
        Posvec  [3*i+1] = float(Particles[i]->x(1));
        Posvec  [3*i+2] = float(Particles[i]->x(2));
        Velvec  [3*i  ] = float(Particles[i]->v(0));
        Velvec  [3*i+1] = float(Particles[i]->v(1));
        Velvec  [3*i+2] = float(Particles[i]->v(2));
        Pressure[i    ] = float(Particles[i]->Pressure);
        Radius  [i    ] = float(Particles[i]->R);
        Tag     [i    ] = int  (Particles[i]->ID);
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Velocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dims[0] = Particles.Size();
    dsname.Printf("Pressure");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tag);


    delete [] Posvec;
    delete [] Velvec;
    delete [] Pressure;
    delete [] Radius;
    delete [] Tag;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);


    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"SPHCenter\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pressure \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey)
{
    Util::Stopwatch stopwatch;
    std::cout << "\n--------------Solving------------------------------------------------------------------------------" << std::endl;

    idx_out = 0;
    double tout = Time;


    ResetInteractions();
    ResetContacts();


    while (Time<tf)
    {

    	// Calculate the acceleration for each particle
        StartAcceleration(Gravity);
        ComputeAcceleration(dt);

         // output
        if (Time>=tout)
        {
            if (TheFileKey!=NULL)
            {
                    String fn;
                    fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                    WriteXDMF    (fn.CStr());
            }
            idx_out++;
            tout += dtOut;
        }

         // Move each particle
         Move(dt);


        // next time position
        Time += dt;

        ResetContacts();
        
    }
}

}; // namespace SPH

#endif // MECHSYS_SPH_DOMAIN_H
