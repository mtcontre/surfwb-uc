# SurfWB-UC #

A two-dimensional shallow water model that accurately reproduces shock-waves and steady states through a shock-capturing well-balanced scheme using boundary fitted curvilinear grids (Guerra et al., 2014). Its core was written by Maricarmen Guerra during her Master's thesis, and the mpi-parallelized version was written by José Galaz. It has had several contributions by Leandro Suarez, and María Teresa Contreras. 

This is the main repository for the code and is currently administrated by José Galaz.

## System requirements
Under linux the requirements are (

* gfortran
* openmpi 1.6.0 or greater (or some mpich)

and for pre and post processing tools

* python 2.7
* numpy
* matplotlib

these can be easily installed with the anaconda-python distribution.

## Installation instructions
First download the software from

    https://bitbucket.org/JoseGalazM/surfwb-uc/downloads/surfwb-uc.zip

an extract all files in a folder, say surfwb-uc in your home directory and then add to your .bashrc

    export SURF=/home/username/surfwb-uc/

and start a new terminal or update your environment variables with

    source .bashrc

## Running an example

If everything is fine you should go to the tests directory with just

    cd $SURF/tests/

and now pick one of them, for example the two dimensional partial dambreak

    cd test2_db2d

if you have all dependencies then you should compile with

     make clean
     make

and for the data you need all python 2.7, numpy and matplotlib to run

     python setrun.py

this will create directories data/ and results/, and inside data/ you can see the bathymetry and initial condition in .png files. Tell the computer where the input data directory is by typing

    export INDIR=data

now we are ready to run the model, say with 4 cores, through

    mpirun -np 4 xsurf

or with just one core just run
    
    ./xsurf

To visualize your results just execute

     python ver.py

which should produce many png files with frames of the simulation in ver/*.png.

## Contribution guidelines

### Test cases

Every new feature must be able to pass each one of the tests in the tests folder, and before merging, each test case must be reconfigured so it can run under the latest version of the code.

### Pull requests

Only shared repository pull requests are allowed. Each new feature must be reflected in a new branch. For more information see [this tutorial](https://es.atlassian.com/git/tutorials/making-a-pull-request/).

<!--
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
* Other community or team contact-->