cd build
rm -rf *
cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX ..
make
