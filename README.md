# Finimizers
Missing definition

# Building
First, pull the submodules with:
```
git submodule init
git submodule update
```
Then, go the SBWT submodule and build it using the instructions in the submodule. Then, compile the experiments with:
```
cd SBWT/build

cmake .. -DCMAKE_C_COMPILER=$(which gcc-10) -DCMAKE_CXX_COMPILER=$(which g++-10)
make -j4

cd ../..
make benchmark --always-make CXX=g++-10
```
# Index construction

The code takes a plain-matrix sbwt file as input. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 30
```

Then, you can build the Finimizers index with:

```
./benchmark build-fmin -u unitigs.fna -i index.sbwt [--lcs LCS.sdsl] -o out-file -t freq --type shortest
```
and query the data with
