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
./benchmark build-fmin -o out-file -u unitigs.fna -i index.sbwt [--lcs LCS.sdsl] -t freq [--type shortest]
```
and query the data with
```
./benchmark search-fmin -o out-file  -q query-file.fa -i index.sbwt [--lcs LCS.sdsl] -f finimizers_bv --unitigs-v finimizers-unitigs_v -t freq [--type shortest]

```
type has to be the same for both commands. The default type is "rarest", t=1. Selecting shortest the shortest finimizers is selected among those with frequency smaller than t. If t=1 then the two types are equivalent. 
