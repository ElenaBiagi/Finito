.PHONY: build_fmin search_fmin benchmark
benchmark:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o benchmark -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp -lz

build_fmin:
	$(CXX) src/build_fmin.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o build_fmin -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp -lz
search_fmin:
	$(CXX) src/search_fmin.cpp SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -O3 -I ./SBWT/include -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o search_fmin -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp
