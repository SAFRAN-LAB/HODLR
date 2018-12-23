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
mkdir ${DEPS_DIR} && cd ${DEPS_DIR} 

# Getting CMake:
CMAKE_URL="http://www.cmake.org/files/v3.3/cmake-3.3.2-Linux-x86_64.tar.gz"
mkdir cmake && wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C cmake
export PATH=${DEPS_DIR}/cmake/bin:${PATH}

# Getting Eigen:
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
tar xjf 3.3.7.tar.bz2
mv eigen-eigen-323c052e1731 eigen
export EIGEN_PATH=$PWD/eigen/
# Returning the the home folder:
cd ..
