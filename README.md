# Finimizers
Missing definition



## Building
First, pull the submodules with:
```
git submodule update --init --recursive
```
Then, go the SBWT submodule and build it using the instructions in the submodule. And compile the experiments with:
```
cd SBWT/build

cmake .. -DCMAKE_C_COMPILER=$(which gcc-10) -DCMAKE_CXX_COMPILER=$(which g++-10)
make -j4

cd ../..
make benchmark --always-make CXX=g++-10
```
## Index construction

The code takes a plain-matrix sbwt file as input generated from canonical unitigs. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 30
```

Then, you can build the Finimizers index with:
```
build-fmin [OPTION...]

  -o, --out-file arg            Output filename.
  -i, --index-file arg          Index input file. This has to be a binary matrix.
  -u, --in-file arg             The query in FASTA or FASTQ format, 
                                possibly gzipped. Multi-line FASTQ is not 
                                supported. If the file extension is .txt, 
                                this is interpreted as a list of query 
                                files, one per line. In this case, 
                                --out-file is also interpreted as a list of 
                                output files in the same manner, one line 
                                for each input file.
  -z, --gzip-output             Writes output in gzipped form. This can 
                                shrink the output files by an order of 
                                magnitude.
      --type arg                Decide which streaming search type you 
                                prefer. Available types:  rarest shortest 
                                optimal verify (default: rarest)
  -t arg                        Maximum finimizer frequency
      --lcs arg                 Provide in input the LCS file if available. 
                                (default: "")
  -h, --help                    Print usage
```

** Modify the example **

```
./benchmark build-fmin -o out-file -u unitigs.fna -i index.sbwt [--lcs LCS.sdsl] -t freq [--type shortest]
```
and query the data with:

```
Query all Finimizers of all input reads.
Usage:
  search-fmin [OPTION...]

  -o, --out-file arg    Output filename.
  -i, --index-file arg  Index input file. This has to be a binary matrix.
  -q, --query-file arg  The query in FASTA or FASTQ format, possibly 
                        gzipped. Multi-line FASTQ is not supported. If the 
                        file extension is .txt, this is interpreted as a 
                        list of query files, one per line. In this case, 
                        --out-file is also interpreted as a list of output 
                        files in the same manner, one line for each input 
                        file.
  -z, --gzip-output     Writes output in gzipped form. This can shrink the 
                        output files by an order of magnitude.
  -t arg                Maximum finimizer frequency
      --type arg        Decide which streaming search type you prefer. 
                        Available types:  rarest shortest optimal verify 
                        (default: rarest)
      --lcs arg         Provide in input the LCS file if available. 
                        (default: "")
  -f, --fmin_bv arg     Provide in input the finimizers binary kmers 
                        vector. (default: "")
  -e, --endpoints arg       Provide in input the endpoints of the 
                            concatenated unitigs. (default: "")
  -g, --global-offsets arg  Provide in input the global offsets of 
                            finimizers in the concatenated unitigs. 
                            (default: "")
  -h, --help            Print usage
```
** Modify the example **
```
./benchmark search-fmin -o out-file  -q query-file.fa -i index.sbwt [--lcs LCS.sdsl] -f fmin_bv --unitigs-v fmin_unitigs -t freq [--type shortest]

```
type has to be the same for both commands. The default type is "rarest", t=1. Selecting shortest the shortest finimizers is selected among those with frequency smaller than t. If t=1 then the two types are equivalent. 

### Canonical unitigs
Canonical unitigs are required as input to build the SBWT index. You can obtain canonical unitigs using [ggcat](https://github.com/algbio/ggcat).

```
ggcat build --min-multiplicity 1 -k 31 --output-file out.fna --threads-count 48 input.fna
```

