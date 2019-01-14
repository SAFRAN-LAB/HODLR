#!/usr/bin/env bash

##############################################################
# This is the install script for HODLRlib                    #
                                                             #
# It allows an easy installation by getting all dependencies #
# and setting all the enviroment variables needed            #
##############################################################

# Setting the environment variable HODLR_PATH:
export HODLR_PATH=$PWD
# Initializing a dependencies folder:
export DEPS_DIR=$HODLR_PATH/deps
if [ ! -d "${DEPS_DIR}" ]; then
    mkdir ${DEPS_DIR}
fi
cd ${DEPS_DIR}

# Getting CMake:
if [ ! -d "cmake/" ]; then
    echo "Setting up the necessary dependencies. Hold on tight! (this may take a few minutes)"
    printf "\nDOWNLOADING CMAKE...\n"
    CMAKE_URL="https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz"
    mkdir cmake && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake
    export PATH=${DEPS_DIR}/cmake/bin:${PATH}
    printf "\nCOMPLETED CMAKE INSTALLATION.\n"
fi

# Getting Eigen:
if [ ! -d "eigen/" ]; then
    printf "\nDOWNLOADING EIGEN...\n"
    wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
    tar xjf 3.3.7.tar.bz2
    mv eigen-eigen-323c052e1731 eigen
    export EIGEN_PATH=$PWD/eigen/
    printf "\nCOMPLETED EIGEN INSTALLATION.\n"
fi

# Returning the the home folder:
cd ..

echo "Writing the necessary enviroment variables to .bash_profile"
echo "# Added by HODLRlib:" >> ~/.bash_profile
to_add="export HODLR_PATH=$HODLR_PATH"
echo $to_add >> ~/.bash_profile
to_add="export PATH=$PATH:${DEPS_DIR}/cmake/bin/"
echo $to_add >> ~/.bash_profile
to_add="export EIGEN_PATH=$EIGEN_PATH"
echo $to_add >> ~/.bash_profile
source ~/.bash_profile
