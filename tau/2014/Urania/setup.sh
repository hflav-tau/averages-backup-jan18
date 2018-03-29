source /cvmfs/sft.cern.ch/lcg/views/LCG_91/x86_64-slc6-gcc62-opt/setup.sh
#rm -rf install build
#mkdir -vp install build
#cd build
#cmake .. -DCMAKE_INSTALL_PREFIX=$(pwd)/../install
#cmake --build . --target install
#cd -
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(pwd)/install/lib64
export PYTHONPATH=${PYTHONPATH}:$(pwd)/install/python
#cd Phys/Tau23Mu/math/3fb
