build_fmin:
	g++ build_fmin.cpp src/SeqIO.cpp src/globals.cpp SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I ./SBWT/include -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o build_fmin -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp

