# Finimizers

# Building
```
cd SBWT/build

cmake .. -DCMAKE_C_COMPILER=$(which gcc-10) -DCMAKE_CXX_COMPILER=$(which g++-10)
make -j4

cd ../..
make benchmark --always-make CXX=g++-10
```
# Running

The code takes a plain-matrix sbwt file as input. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 30
```

Then, you can run the benchmark with:

```
./benchmark index.sbwt .....
```
