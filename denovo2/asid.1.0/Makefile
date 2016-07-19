all:
	g++ -O3 -march=nocona -pipe -ffast-math   -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -c src/util.cxx -o src/util.o
	g++ -O3 -march=nocona -pipe -ffast-math   -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -c src/zutil.cxx -o src/zutil.o
	g++ -O3 -march=nocona -pipe -ffast-math   -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE src/asid_light.cpp src/util.o src/zutil.o -o asid_light

clean: 
	rm -f src/*.o
	rm -f asid_light
