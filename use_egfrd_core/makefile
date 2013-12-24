epdp_path=/home/sakamoto/local/include/egfrd_core
gsl_path=/usr/lib

static_libs=/usr/lib/libgslcblas.a /usr/lib/libgsl.a

build:
	g++ -DHAVE_CONFIG_H -DHAVE_INLINE -I${epdp_path}  -lm -legfrd_core ${static_libs} mymapk.cpp -o mymapk

clean:
	rm mymapk
