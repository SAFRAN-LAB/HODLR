#!/usr/bin/env bash

##############################################################
# This is the install script for HODLRlib                    #
                                                             #
# It allows an easy installation by getting all dependencies #
# and setting all the enviroment variables needed            #
##############################################################

# Setting the environment variable HODLR_PATH:
HODLR_PATH=$PWD
# Initializing a dependencies folder:
DEPS_DIR=$HODLR_PATH/deps
if [ ! -d "${DEPS_DIR}" ]; then
    mkdir ${DEPS_DIR}
fi
cd ${DEPS_DIR}

# Getting CMake:
# In our experience, building CMake on a Mac and then using this CMake seems to throw an error
# We find that everything works as expected when using the system CMake installed using brew
# TODO: Look into further
CMAKE_INSTALLED_NOW=0
if [ "$(uname)" == "Darwin" ]; then
    printf "\nPlease ensure that you have installed CMake > 3.12 (preferrably from brew)...\n"
else
    if [ ! -d "cmake/" ]; then
        CMAKE_INSTALLED_NOW=1
        echo "Setting up the necessary dependencies. Hold on tight! (this may take a few minutes)"
        printf "\nDOWNLOADING CMAKE...\n"
        CMAKE_URL="https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz"
        mkdir cmake && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake
        PATH=${DEPS_DIR}/cmake/bin:${PATH}
        printf "\nCOMPLETED CMAKE INSTALLATION.\n"
    fi
fi

# Getting Eigen:
if [ ! -d "eigen/" ]; then
    printf "\nDOWNLOADING EIGEN...\n"
    wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
    tar xjf 3.3.7.tar.bz2
    mv eigen-eigen-323c052e1731 eigen
    EIGEN_PATH=$PWD/eigen/
    printf "\nCOMPLETED EIGEN INSTALLATION.\n"
fi

# Returning the the home folder:
cd ..

echo "Writing the necessary enviroment variables to .bash_profile, if not already defined..."

[ -z "$HODLR_PATH" ]
if env | grep -q ^HODLR_PATH=
then
  echo HODLR_PATH is already set
else
  echo HODLR_PATH updated in .bash_profile
  echo "# Added by HODLRlib:" >> ~/.bash_profile
  to_add="export HODLR_PATH=$HODLR_PATH"
  echo $to_add >> ~/.bash_profile
fi

if [ "$CMAKE_INSTALLED_NOW" == "1" ]; then
    echo PATH updated in .bash_profile to point to newly installed CMake
    to_add="export PATH=$PATH:${DEPS_DIR}/cmake/bin/"
    echo $to_add >> ~/.bash_profile
fi

[ -z "$EIGEN_PATH" ]
if env | grep -q ^EIGEN_PATH=
then
  echo EIGEN_PATH is already set
else
  echo EIGEN_PATH updated in .bash_profile
  to_add="export EIGEN_PATH=$EIGEN_PATH"
  echo $to_add >> ~/.bash_profile
fi

source ~/.bash_profile
