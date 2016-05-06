#####################################################################################
# PersianSPH - A C++ library to simulate Mechanical Systems                         #
# Copyright (C) 2016 Maziar Gholami Korzani and Sergio Galindo-Torres               #
#                                                                                   #
# This file is part of PersianSPH                                                   #
#                                                                                   #
# PersianSPH is free software; you can redistribute it and/or modify it under the   #
# terms of the GNU General Public License as published by the Free Software         #
# Foundation; either version 2 of the License, or (at your option) any later        #
# version.                                                                          #
#                                                                                   #
# PersianSPH is distributed in the hope that it will be useful, but WITHOUT ANY     #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A   #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.          #
#                                                                                   #
# You should have received a copy of the GNU General Public License along with      #
# PersianSPH; if not, write to the Free Software Foundation, Inc.,                  #
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA                       #
#####################################################################################

#!/bin/bash

#checking Mercurial
command -v hg >/dev/null && 
	echo "Mercurial has been found." ||
	{
	echo "Mercurial is not installed."; 
	while true; do
		read -p "Do you wish to install this program (sudo access required)?(y/n)" yn;
    		case $yn in
			[Yy]* ) sudo apt-get install mercurial;
				command -v hg >/dev/null || { echo "You do not have sudo access."; exit; }
				break;;
			[Nn]* ) echo "Mercurial is compulsory to download PersianSPH code."; exit;;
			* ) echo "Please answer yes or no.";;
		esac
	done
	}

#checking g++ and cmake
for $Compulsory_Software in g++ cmake
do
	command -v $Compulsory_Software >/dev/null && 
		echo "$Compulsory_Software has been found." ||
		{
		echo "$Compulsory_Software is not installed."; 
		while true; do
			read -p "Do you wish to install this program (sudo access required)?(y/n)" yn;
	    		case $yn in
				[Yy]* ) sudo apt-get install $Compulsory_Software;
					command -v $Compulsory_Software >/dev/null || { echo "You do not have sudo access."; exit; }
					break;;
				[Nn]* ) echo "$Compulsory_Software is compulsory to install PersianSPH code."; exit;;
				* ) echo "Please answer yes or no.";;
			esac
		done
		}
done

#using the current directory as the main installation directory
Current_Address=$(pwd)

#copying the code on this pc 
if [ -d "$Current_Address/sph" ]; then
	echo "PersianSPH code directory already exists at $Current_Address/sph"
	while true; do
		read -p "Do you wish to remove it?(y/n)" yn;
    		case $yn in
			[Yy]* ) rm -rf $Current_Address/sph; break;;
			[Nn]* ) echo "It is probably an older copy of PersianSPH code on this PC!!!"; break;;
			* ) echo "Please answer yes or no.";;
		esac
	done
else
	hg clone https://mghkorzani@bitbucket.org/persiansph/sph
fi

#setting a permanent environmental variable for the main code 
SPH_Address=$SPH
if [ -z "${SPH_Address}" ]; then
	echo "Creating an environment variable for the path of the SPH code as \$SPH"
	echo "export SPH=$Current_Address/sph" >> $HOME/.bashrc
	SPH_Address="$Current_Address/sph"
	source $HOME/.bashrc
fi

#copying packages on this pc 
if [ -d "$Current_Address/pkg" ]; then
	echo "Libraries' directory already exists at $Current_Address/pkg"
	while true; do
		read -p "Do you wish to remove it?(y/n)" yn;
    		case $yn in
			[Yy]* ) rm -rf $Current_Address/pkg; break;;
			[Nn]* ) echo "These are probably older copies of required libraries on this PC!!!"; break;;
			* ) echo "Please answer yes or no.";;
		esac
	done
else
	hg clone https://mghkorzani@bitbucket.org/persiansph/pkg
fi

#setting a permanent environmental variable for the packages 
PKG_Address=$PKG
if [ -z "${PKG_Address}" ]; then
	echo "Creating an environment variable for the path of required packages as \$PKG"
	echo "export PKG=$Current_Address/pkg" >> $HOME/.bashrc
	PKG_Address="$Current_Address/pkg"
	source $HOME/.bashrc
fi

#installing Blitz library
if [ ! -d "$Current_Address/pkg/blitz-0.9" ]; then
	echo "... Installing Blitz library ..."
	cd $PKG_Address
	tar -xzf blitz-0.9.tar.gz
fi

# installing other libraries
while true; do
	read -p "Do you have sudo access to install libraries?(y/n)"  yn;
	case $yn in
		[Yy]* ) echo "... Installing Blas, Lapacke, Gsl, HDF5, g++ and Cmake libraries ..."
			sudo apt-get install libblas-dev liblapack-dev libgsl0-dev libhdf5-serial-dev;
			break;;
		[Nn]* ) echo "... Unpacking libraries ...";
			cd $PKG_Address;
			tar -xzf hdf5-1.8.11.tar.gz;
			tar -xzf gsl-2.1.tar.gz;
			tar -xzf lapack-3.5.0.tgz;
			echo "... Compiling libraries ...";
			read -p "It takes a few minutes, press any key to continue ..."
			cd $PKG_Address/hdf5-1.8.11;
			./configure --prefix=$PKG_Address/hdf5-1.8.11;
			make;
			make install;
			cd $PKG_Address/gsl-2.1;
			./configure --prefix=$PKG_Address/gsl-2.1;
			make;
			cd $PKG_Address/lapack-3.5.0;
			cmake .;
			make;
			break;;
		* ) echo "Please answer yes or no.";;
	esac
done

echo ""
echo "Installation is completed."
echo "Close all terminal windows and open a new one to take effect the defined environment variables" 
echo "Please refer to the tutorial to run a simulation."



