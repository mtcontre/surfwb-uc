# SurfWB-UC #

A two-dimensional shallow water model that accurately reproduces shock-waves and steady states through a shock-capturing well-balanced scheme (Guerra et al., 2014). It was initially developed by Maricarmen Guerra, but has had several additions from Leandro Suarez, José Galaz and María Teresa Contreras. This is the main repository for the code.

## System requirements
Under linux the requirements are

* gfortran
* python 2.7
* numpy
* matplotlib

## Installation instructions
First clone the repository in your home directory with

    git clone https://JoseGalazM@bitbucket.org/JoseGalazM/surfwb-uc.git

then add to your .bashrc

    export SURF=/home/username/surfwb-uc/

and start a new terminal or update your environment variables with

    source .bashrc

## Running an example

If everything is fine you should go to the tests directory with just

    cd $SURF/tests/

and now pick one of them, for example the two dimensional partial dambreak

    cd test2_db2d

if you have the gnu fortran compiler gfortran, then try compiling with

     make clean
     make

now you need all python 2.7, numpy and matplotlib to run

     python setrun.py

this will create directories data/ and results/, and inside data/ you can see the bathymetry and initial condition in .png files. Now we are ready to run

    ./xsurf

and to see some results, you can try

     python ver.py

which should produce many png files with frames of the simulation in vis/*.png.



### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact