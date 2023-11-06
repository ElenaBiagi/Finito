.PHONY: benchmark
benchmark:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o benchmark -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp -lz -D MAX_KMER_LENGTH=250

benchmark_debug:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o benchmark -Wno-deprecated-declarations -march=native -fopenmp -lz -D MAX_KMER_LENGTH=250
