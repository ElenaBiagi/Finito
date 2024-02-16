# Finimizers
This is the code for the paper **Finimizers: Variable-length bounded-frequency minimizers for $k$-mer sets** by J. N. Alanko, E. Biagi,  S. J. Puglisi. 

### Shortest Unique Finimizers
Let $G$ be the de Bruijn graph of a set of $k$-mers $R$, $t \geq 1$ be an integer, $X$ be a $k$-mer, and $Y$ be a substring of $X$. We say $Y$ is a **shortest $t$-finimizer** of $X$ with respect to the $k$-mer set $R$ if $Y$ has at most $t$ occurrences in $G$ and there does not exist a shorter substring of $X$ with at most $t$ occurrences in $G$. If $t = 1$ then we say $Y$ is a _shortest-unique finimizer_ of $X$.


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
./SBWT/build/bin/sbwt build -i SBWT/example_data/coli3.fna -o index.sbwt -k 31 --add-reverse-complements
```

You also need to create a file with the reverse complement of the unitigs.

```
./benchmark reverse -i unitigs.fna -o rev_unitigs.fna
```

```
Usage:
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
./benchmark build-fmin -o out-file -f unitigs.fna -r rev_unitigs.fna -i index.sbwt [--lcs LCS.sdsl] [-t 1] [--type rarest] 
```
```
Usage:
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

## k-mer lookup queries

You can query $k$-mer in the unitigs with:
```
./benchmark search-fmin -o <outfile>  -q <query-file.fa> -i index.sbwt [--lcs LCS.sdsl] -f fmin_bv --unitigs-v fmin_unitigs -t freq 
```
```
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
type has to be the same for both commands.
Support for lookup queries is currently available only for "rarest".

### Disjoint Spectrum Preserving String Set (DSPSS)
A DSPSS is required as input to build the SBWT index. You can obtain canonical unitigs or eulertigs using [ggcat](https://github.com/algbio/ggcat).

```
ggcat build --min-multiplicity 1 -k <k> --output-file unitigs.fna --threads-count 48 input.fna
```

### Unitigs flipping
To reduce the space usage it is advisable to flip the unitigs with [unitig-flipper](https://github.com/jnalanko/unitig_flipper).

```
unitig_flipper --input unitigs.fna --output flipped_unitigs.fna -k <k>

```
