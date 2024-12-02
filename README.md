# pyANTI-FASc

## 1. General description

Python version of ANTI-FASc

pyANTI-FASc, acronym for python Automatic Numerical Tsunami Initial conditions: on-the-Fly rupture Areas and earthquake Scenarios, is a software allowing the fast computation of large ensembles of slip distributions on complex non-planar fault interfaces (Herrero & Murphy 2018 GJI, Maesano et al. 2017, Sci. Rep.; Tonini et al. 2021, GJI) like the subducting plates. These slip models can be promptly used as initial conditions for the computation of tsunami scenarios in the framework of both Seismic-Probabilistic Tsunami Harzard Assessment (see Scala et al. 2020 PAGEOPH - Basili et al. 2021 Frontiers) and for real-time Probabilisitic Tsunami Forecasting (see Selva et al. 2021 Nature). The newest release (version v1.0.0) is also available at the following DOI: *https://zenodo.org/doi/10.5281/zenodo.13614657* along with all the previous realeses. 

IMPORTANT: Please refer the repository and cite the DOI as indicated at the zenodo webpage, in your pubblications, if you use the software for your research studies. Cite it as:

Antonio Scala, Manuel Mojica, & rissclab-tester. (2024). antonioscalaunina/pyANTI-FASc: pyANTIFASc_1.0.0 (v1.0.0). Zenodo. *https://doi.org/10.5281/zenodo.13614657*


A wiki documentation (currently under construction) will be soon available at this [link](https://github.com/antonioscalaunina/pyANTI-FASc/wiki)

The software is composed by three modules managed by a single standalone python executable that runs the module sequentially, also managing the transfer of the input/output files between the three modules. The software can be run also through the use of a Jupyter Notebook showing the intermediate step processes and outputs
A simple postprocessing tool is also provided to convert the output files in a standard georeference format and plot the slip distributions.

The three modules can be summarized as follows

1- **Preprocess module** 

This module includes a mesh generator that generates input mesh file from a nodes and faces discretization. A set of pre-computed discretization of the main subducting slabs worldwide is also provided. These mesh discretizations are based on [Slab 2.0](https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467). These meshes are composed by triangular elements having sides in a range between 10 and 15 km.  The available modelled slabs, released with this version of the software is shown [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/utils/sz_slabs/map_of_slabs.png).

On the selected meshes some preliminary processes are performed during the execution of this module, like the computation of the area of the elements, the definition of element-to-element connection and a preliminary estimation of the interdistance between the mesh nodes
    

2 - **Rupture areas computation**:
    
This module computes a set of possible rupture geometry on the selected fault mesh. It has two different use mode:
         
   - **Hazard**: it computes a large set of possible different rupture areas in all the prescribed magnitude bins to cover in a homogeneous way the whole provided meshed zone
         
   - **PTF**: it computes all the scenarios “compatible” with estimation and uncertainty of magnitude and location for a given earthquake

The application selection, as well as the setting of the other parameters will be shown in more details in the proposed example and in the wiki documentation (currently under construction)


3 - **k223d:** 

This module, based on the original software presented in Herrero & Murphy (2018, GJI) and available at *https://github.com/s-murfy/k223d* performs:

   - A refined computation of inter-distance between nodes to be used in k-square slip distribution computation. It is based on the lateration algorithm presented in Herrero & Murphy (2018, GJI)
   
   - The computation of ensembles of stochastic k-square slip distributions for all the previously selected areas also accounting for other conditions (e.g. homogeneous or variable rigidity, surface slip amplification, variable stress-drop etc.)
  
**IMPORTANT:** More details about the k223d module and its original sources can be found in the repository of ANTI-FASc at the following link *https://github.com/antonioscalaunina/ANTI-FASc/blob/main/src/k223d/README.md*. The use of this module is shown in the examples in the main folder.

Below are the instructions for installing the software dependencies. Please refer to the examples in the main folder and to the wiki documentation **(BOTH UNDER CONSTRUCTION)** for further details regarding the code functionality, the configuration of the input files and the database of precomputed mesh discretizations

## 2 Installation

### 2.1 Linux & Windows WSL

#### 2.1.1 Linux environment

The most practical way to run the code in a Linux environment is to create a conda environment where all the necessary libraries and dependencies are available. If you do not have a conda distributions already installed you might install miniconda, a minimal distribution of conda, through the following commands

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    ./Miniconda3-latest-Linux-x86_64.sh
    conda update

The first command will download the installer while the other ones are used to run the installer and install conda to the latest version. During the installation the user will be asked to decide if the conda execution must be automatic every time a terminal is open. Such a choice is optional.

Once Miniconda is installed, simply navigate to the main folder and create the Conda environment antifasc by entering the following commands:
    
    cd pyANTI-FASc
    conda env create -f ANTIFASc.yml

at the end of installation you may enter into the environment typing the command:

    conda activate antifasc

Alternatively, you can create the conda environment by using the following commands:

    conda create --name antifasc python=3.9
    conda activate antifasc
    pip install -r requirements.txt

Once the environment has been created and activated the fortran executable *k223d.x* available in the folder [pyANTI-FASc/bin](https://github.com/antonioscalaunina/pyANTI-FASc/bin/) must be copied among the executable scripts of the environment. It might be necessary to give the execution permission through the following commands:

    cd ~/pyANTI-FASc/bin
    cp k223d.x ~/miniconda3/envs/antifasc/bin/
    chmod +x ~/miniconda3/envs/antifasc/bin/
    
Please note that if you are using a different conda distribution the folder of the environment might have a different name

    sudo apt-get install gfortran

Then for compiling:

    cd pyANTIFASc/src/k223d
    make
    cp k223d.x ../../bin

If everything worked you are now ready to have fun with pyANTIFASc (😉). Try and run the [Example1_Tohoku](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md) to start!

#### 2.1.2 Windows wsl

If you use the wsl distribution on Windows to work in a virtual Ubuntu environment you can use the same guide for installation just presented in the [section 2.1.1](#211-linux-environment). 
To install the wsl distributions you might follow the instructions at this [link](https://learn.microsoft.com/en-us/windows/wsl/install). The easiest way is to open the Windows PowerShell and digit:

    wsl --install

This command will enable the features necessary to run WSL and install the Ubuntu distribution of Linux. Once the wsl is installed all the instructions for [Linux distributions](#211-linux-environment) can be used


### 2.2 Windows through conda GUI

 This installation is a bit more complicated than the one described for Linux and WSL environments but it will allow to run the code in a fully Windows environment, e.g. using a conda GUI, like Anaconda Navigator. This also allows to run the code through the Jupyter Notebook available [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb). To do that you should follow the steps outlined in the next subsections.

#### 2.2.1 Create the conda environment

1 - Download the repository at the main page [https://github.com/antonioscalaunina/pyANTI-FASc/tree/main](https://github.com/antonioscalaunina/pyANTI-FASc/tree/main) or with the direct link [https://github.com/antonioscalaunina/pyANTI-FASc/archive/refs/heads/main.zip](https://github.com/antonioscalaunina/pyANTI-FASc/archive/refs/heads/main.zip) and unzip it.

2 - Open a PowerShell, within the Conda GUI you are using, and within the main folder of the repository type the following command:

    conda env create -f ANTIFASc.yml

this command creates the conda environment *antifasc* (installing Python 3.9.16 within it) and installs all the needed libraries and dependencies within it. Then enter into the environment by typing:

    conda activate antifasc

or searching for the environment antifasc through the menu *Environments* of the Conda GUI.

Alternatively, you can create the conda environment by typing the following commands in Conda GUI PowerShell:

    conda create --name antifasc python=3.9
    conda activate antifasc
    pip install -r requirements.txt

#### 2.2.2 Install the fortran compiler through MinGW

The following steps will allow to download a fortran compiler for the Windows environment. This is necessary to compile the fortran module of the software.

1 - At the webpage [https://sourceforge.net/projects/mingw/](https://sourceforge.net/projects/mingw/) download the MinGW - Minimalist GNU for Windows installer

2 - Run the installer

3 - If the MinGW Installation Manager automatically opens at the end of the installation, close it.

4 - Search for the MinGW Installation Manager among the available apps (e.g. through the Windows menu) and run it.

5 - From the available menu select all the packages having **gcc-fortran** in the name

6 - From the menu on the top select Installation > Apply changes

7 - Wait for the end of the process installing the gfortran compiler and the needed libraries

#### 2.2.3 Add the fortran compiler to the enviromental variables

1 - Enter in Advanced system settings from This Pc menu: on *"This PC"* menu right-click button select *Properties* and then *Advanced system settings* 

2 - Within the pop-up window that appears select *Enviromental Variables* and on the next window double-click on *Path* in the *System variables* menu 

3 - Select *new* and add the folder C:\MinGW\bin (**WARNING**: C:\ should be the standard position where MinGW is automatically installed, but it could be elsewhere!)

4 - Repeat the step 3 to add the folders C:\MinGW and C:\MinGW\lib\gcc\mingw32\6.3.0\ (**WARNING**: 6.3.0 should be the installed version. However if you install a different version, the name of the last folder will change accordingly)

5 - Click OK to close all the windows and restart your computer to make this modifications effective.

#### 2.2.4 Compile the k223d module

1 - After restarting, first check if your fortran compiler is properly installed opening the Windows prompt and typing:

    gfortran --version

2 - If it is installed and correctly configured you should see a message like this:

    GNU Fortran (MinGW.org GCC-6.3.0-1) 6.3.0
    Copyright (C) 2016 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

3 - Enter in the main folder of pyANTI-FASc and then move into the *src\k223d* folder:

    cd src\k223d\

4 - Compile the code and move the executable in the bin folder with the following commands:

    gfortran -c forparse.f90 utils.f90 lateration.f90 typedef_Win.f90 makepdf.f90 k223d.f90
    gfortran -o k223d.x forparse.o utils.o lateration.o typedef_Win.o makepdf.o k223d.o
    copy k223d.x  ..\..\bin

If all these steps went fine, now you are ready to have fun with ANTI-FASc (😉). Try and run the [Example1_Tohoku](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md) to start!

### 2.3 Mac

If you use a macOS you might create a terminal similar to the one that you can use on Linux through different package managers for macOS, which allow users to easily install and manage software packages from the command line of the Mac terminal. The two most popular package managers are Homebrew and MacPorts

#### 2.3.1 Conda and gfortran with Homebrew

To firstly install Homebrew open the terminal and run the following command:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew update

Through brew, you might install gfortran and conda (miniconda distribution) with the following commands:

    brew install gcc
    brew install --cask miniconda

you might verify the installed version of gfortran with the following command:

    gfortran --version

If it does not work it might be necessary to restart the terminal.
After installation of miniconda follow the instructions to initialize conda (you may need to restart the shell):

    conda init

and verify the installation checking the version:

    conda --version

See now section [2.3.3](#233-antifasc-on-mac-terminal) to complete the ANTIFASc installation

#### 2.3.2 Conda and gfortran with MacPorts

To firstly install MacPorts go the MacPorts [webpage](https://www.macports.org/install.php) and download the right version of installer for your OS. Then launch the *.pkg* and complete the installation following the instructions

MacPorts requires Xcode and its Command Line Tools. You can install them via:

   xcode-select --install

After installation you might update MacPorts with the command:

    sudo port selfupdate

Through port you might install gfortran and verify the version with the following command:

    sudo port install gcc12
    gfortran-mp-12 --version

This will install GCC 12, which includes gfortran. You may replace gcc12 with another version if necessary (e.g., gcc11, gcc10, etc.). In MacPorts, gfortran is typically suffixed with the GCC version, like gfortran-mp-12. We suggest to create an alias for easier usage. editing the ~/.bash_profile or ~/.bashrc adding the following line:

    alias gfortran='/usr/local/bin/gfortran-12'

Alternatively, you could modify the [makefile](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/src/k223d/makefile) on the line 4 replacing gfortran with gfortran-mp-12 (or the other installed version) 

MacPorts does not have a direct package for conda, but you can download and install it manually. Download the right version at this [webpage](https://www.macports.org/install.php)
Run the installer script:
    
    bash Miniconda3-latest-MacOSX-x86_64.sh

After installation, follow the on-screen prompts and initialize conda

    conda init

and verify the installation checking the version:

    conda --version

See now section [2.3.3](#233-antifasc-on-mac-terminal) to complete the ANTIFASc installation

#### 2.3.3 ANTIFASc on MAC terminal

After the installations described in the previous sub-sections the steps are very similar to the ones for linux. Enter into the main folder and create the conda environment:

    cd pyANTIFASc
    conda env create -f ANTIFASc.yml

at the end of installation you may enter into the environment by typing:

    conda activate antifasc

Alternatively, you can create the conda environment by typing:

    conda create --name antifasc python=3.9
    conda activate antifasc
    pip install -r requirements.txt

Then compile k223d and copy the executable to the bin folder:

    cd pyANTIFASc/src/k223d
    make
    cp k223d.x ../../bin

If everything worked you are now ready to have fun with pyANTIFASc (😉). Try and run the [Example1_Tohoku](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/Example1_Tohoku.md) to start!

### 2.4 Run the Jupyter Notebook

Regardless of the operating system, the software can always be run through the Jupyter Notebook available [here](https://github.com/antonioscalaunina/pyANTI-FASc/blob/main/bin/antifasc_main.ipynb). For example, you can install an IDE like [Visual Studio Code](https://code.visualstudio.com/download)  to run the Jupyter Notebook on Ubuntu, macOS, or Windows. 
**IMPORTANT**: To ensure everything works correctly, it's essential to compile the Fortran code *k223d* using the appropriate compiler for your operating system, as outlined in the previous sections.


## 3 ACKNOWLEDGEMENTS

A special thanks to Stefano Lorito, Fabrizio Romano, Manuela Volpe, Hafize Basak Bayraktar, Jacopo Selva, Gaetano Festa and Antonio Giovanni Iaccarino for participating at the different phases of conceiving, revising, developing and testing of the current version of the platform.

Thanks to Roberto Basili, Francesco Emanuele Maesano, Mara Monica Tiberti and Gareth Davies for their precious contribution in providing detailed geometrical models for some of the slabs included in the database distributed along with this platform.

Thanks to Shane Murphy and Andre Herrero for their valuable contribution in developing and training to use the softwares composing the third module of the platform  


## 4 BIBLIOGRAPHY

Basili R. et al. (2021), The Making of the NEAM Tsunami Hazard Model 2018 (NEAMTHM18), Frontiers in Earth Science, DOI: 10.3389/feart.2020.616594 

Herrero and Murphy (2018), 	Self-similar slip distributions on irregular shaped faults, Geophysical Journal International, DOI: 10.1093/gji/ggy104

Maesano, F.E., Tiberti, M.M. and Basili, R. (2017), The Calabrian Arc: Three-dimensional modelling of the subduction interface, Scientific Reports DOI: 10.1038/s41598-017-09074-8

Murotani, S., Satake, K., and Fujii, Y. (2013), Scaling relations of seismic moment, rupture area, average slip, and asperity size for M~9 subduction-zone earthquakes, Geophysical Research Letters, 40, DOI: 10.1002/grl.50976.

Scala A. et al. (2020), Effect of Shallow Slip Amplification Uncertainty on Probabilistic Tsunami Hazard Analysis in Subduction Zones: Use of Long-Term Balanced Stochastic Slip Models, Pure and Applied Geophysics, DOI: 10.1007/s00024-019-02260-x

Selva, J., Lorito, S., Volpe, M. et al. (2021). Probabilistic tsunami forecasting for early warning. Nat Commun 12, 5677 (2021). DOI: 10.1038/s41467-021-25815-w

Strasser, F. O., Arango, M. C., & Bommer, J. J. (2010), Scaling of the source dimensions of interface and intraslab subduction-zone earthquakes with moment magnitude. Seismological Research Letters, DOI: 10.1785/gssrl.81.6.941.

Tonini et al. (2020), Importance of earthquake rupture geometry on tsunami modelling: The Calabrian Arc subduction interface (Italy) case study, Geophysical Journal International, DOI: 10.1093/gji/ggaa409
