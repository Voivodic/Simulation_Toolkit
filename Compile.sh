gcc Translate.c -o Translate

gcc Split.c -o Split

$HOME/Documentos/CGAL/cgal-releases-CGAL-4.13.1/Scripts/scripts/cgal_create_CMakeLists -s Critical
cmake -DCGAL_DIR=$HOME/Documentos/CGAL/cgal-releases-CGAL-4.13.1/ .
make

gcc HVFinder.c -o HVFinder -lm -fopenmp

gcc Profiles.c -o Profiles -lm -fopenmp

gcc Stacking.c -o Stacking -lm

gcc PowerSpectrum.c -o PowerSpectrum -I$HOME/local/include -IFFTLog-master/include/ -L$HOME/local/lib -lm -lfftw3


