#!/bin/bash -x

BUILD_DIR=$(pwd)
echo $BUILD_DIR

cat > $BUILD_TAG.pbs << EOF
#PBS -A MAT151
#PBS -N $BUILD_TAG
#PBS -j oe
#PBS -l walltime=1:00:00,nodes=1
#PBS -d $BUILD_DIR
#PBS -l partition=rhea

cd $BUILD_DIR

source /sw/rhea/environment-modules/3.2.10/rhel6.7_gnu4.4.7/init/bash

module swap PE-intel PE-gnu/5.3.0-1.10.2
module load git

env

module list

cd build 

cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..

make

echo 
echo checking J1
echo ----------------------------------------------------
echo 

./bin/diff_j1

echo 
echo checking J2
echo ----------------------------------------------------
echo 

./bin/diff_j2

echo 
echo checking JeeI
echo ----------------------------------------------------
echo 

./bin/diff_jeeI

EOF

cp $BUILD_TAG.pbs $BUILD_DIR

cd $BUILD_DIR

source scl_source enable rh-python35
which python 

$BUILD_DIR/../../../scripts/blocking_qsub.py $BUILD_DIR $BUILD_TAG.pbs

grep 'Fail' ../$BUILD_TAG.o*
