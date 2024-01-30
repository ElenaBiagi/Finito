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

cmake .. -DCMAKE_C_COMPILER=$(which gcc-10) -DCMAKE_CXX_COMPILER=$(which g++-10) -D MAX_KMER_LENGTH=250
make -j4

cd ../..
make benchmark --always-make CXX=g++-10
```
## Index construction

The code takes a plain-matrix sbwt file as input generated from canonical unitigs. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 30 --add-reverse-complements
```

You also need to create a file with the reverse complement of the unitigs.

```
./benchmark reverse -i unitigs.fna -o rev_unitigs.fna
```

```
  reverse -i <input> -o <output>

  -i, --in-file arg   The SPSS in FASTA or FASTQ format, possibly gzipped. 
                      Multi-line FASTQ is not supported. If the file 
                      extension is .txt, this is interpreted as a list of 
                      query files, one per line. In this case, --out-file 
                      is also interpreted as a list of output files in the 
                      same manner, one line for each input file.
  -o, --out-file arg  Reverse complement files. (default: out.fna)
  -h, --help          Print usage
```

Then, you can build the Finimizers index with:
```
build-fmin [OPTION...]

  -o, --out-file arg    Output index filename prefix.
  -i, --index-file arg  SBWT file. This has to be a binary matrix.
  -f, --f-file arg      The unitigs in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the 
                        file extension is .txt, this is interpreted as a list of query files, one per line. In this case, 
                        --out-file is also interpreted as a list of output files in the same manner, one line for each input 
                        file.
  -r, --r-file arg      reverse complement of f-file
      --type arg        Available types:  rarest
                        (default: rarest)
  -t arg                Maximum finimizer frequency
      --lcs arg         Provide in input the LCS file if available. 
                        (default: "")
  -h, --help            Print usage
```

```
./benchmark build-fmin -o out-file -u unitigs.fna -i index.sbwt [--lcs LCS.sdsl] [-t 1] [--type rarest]
```
and query the data with:

```
Query all Finimizers of all input reads.
Usage:
  search-fmin [OPTION...]

  -o, --out-file arg    Output filename, or stdout if not given.
  -i, --index-file arg  Index filename prefix.
  -q, --query-file arg  The query in FASTA or FASTQ format, possibly 
                        gzipped. Multi-line FASTQ is not supported. If the 
                        file extension is .txt, this is interpreted as a 
                        list of query files, one per line. In this case, 
                        --out-file is also interpreted as a list of output 
                        files in the same manner, one line for each input 
                        file.

  -h, --help            Print usage
```
** Modify the example **
```
./benchmark search-fmin -o <outfile>  -q <query-file.fa> -i index.sbwt [--lcs LCS.sdsl] -f fmin_bv --unitigs-v fmin_unitigs -t freq 

```
type has to be the same for both commands. The default type is "rarest", t=1. Selecting shortest the shortest finimizers is selected among those with frequency smaller than t. If t=1 then the two types are equivalent.
Support for lookup queries is currently available only for "rarest".

### Spectrum Preserving String Set (SPSS)
A SPSS is required as input to build the SBWT index. You can obtain canonical unitigs using [ggcat](https://github.com/algbio/ggcat).

```
ggcat build --min-multiplicity 1 -k <k> --output-file unitigs.fna --threads-count 48 input.fna
```

### Unitigs flipping
To reduce the space usage it is advisable to flip the unitigs with [unitig-flipper](https://github.com/jnalanko/unitig_flipper).

```
unitig_flipper --input unitigs.fna --output flipped_unitigs.fna -k <k>

```
