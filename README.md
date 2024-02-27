# Finito
This is the code for the paper [**Finimizers: Variable-length bounded-frequency minimizers for $k$-mer sets**](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1) by J. N. Alanko, E. Biagi,  S. J. Puglisi. 

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
```
Select the desired branch: *main*: sigle index, *double*: double index (+reverse complements).
The following instructions are for the main branch.
```
make benchmark --always-make CXX=g++-10
```
## Single index construction

The code takes a plain-matrix SBWT file as input generated from canonical unitigs. You can generate one by running:

```
./SBWT/build/bin/sbwt build -i <unitigs.fna> -o <index.sbwt> -k <31> 
```

Then, you can build the Finimizers index with:

```
./benchmark build-fmin -o <finimizer-index>  -i <index.sbwt> -u <unitigs.fna> [--lcs LCS.sdsl] [-t 1] [--type rarest] 
```
```
Usage:
build-fmin [OPTION...]

  -o, --out-file arg    Output index filename prefix.
  -i, --index-file arg  SBWT file. This has to be a binary matrix.
  -u, --in-file arg     The unitigs in FASTA or FASTQ format, possibly gzipped.
                        Multi-line FASTQ is not supported.
      --type arg        Decide which streaming search type you prefer. 
                        Available types:  rarest shortest verify.
                        The latter two only provide some stats. (default: rarest)
  -t arg                Maximum finimizer frequency (default: 1)
      --lcs arg         Provide in input the LCS file if available. 
                        (default: "")
  -h, --help            Print usage
```

## k-mer localization queries

You can query $k$-mer in the unitigs with:
```
./benchmark search-fmin -o <out-file>  -i <finimizer-index> -q <query-file.fa> 
```
```
Usage:
  search-fmin [OPTION...]

  -o, --out-file arg    Output filename, or stdout if not given.
  -i, --index-file arg  Index filename prefix.
  -q, --query-file arg  The query in FASTA or FASTQ format, possibly gzipped.
                        Multi-line FASTQ is not supported.
  -h, --help            Print usage
```
Support for localization queries is currently available only for "rarest".
The output for each kmer is a pair (unitig id, index) or (-1,-1) if not found.


### Disjoint Spectrum Preserving String Set (DSPSS)
A DSPSS is required as input to build the SBWT index. You can obtain canonical unitigs or eulertigs using [ggcat](https://github.com/algbio/ggcat).

```
ggcat build --min-multiplicity 1 -k <k> --output-file <unitigs.fna> --threads-count 48 <input.fna>
```

### Unitigs flipping
To reduce the space usage it is advisable to flip the unitigs with [unitig-flipper](https://github.com/jnalanko/unitig_flipper).

```
unitig_flipper --input <unitigs.fna> --output <flipped_unitigs.fna> -k <k>

```
