# Simulation_Toolkit

Some useful codes in the study of halos and voids in N-body simulations. The codes here can find halo and voids in dark matter simulations, compute their density and velocity profiles, stack them, and compute the linear biases using the halo/void-matter and the auto power spectrum.

This folder contains the following codes:

-Compile.sh: This is a simple bash file used to compile the other codes;

-Translate.c: This code translates the file with the positions and velocities of the particles to the format accepted by the other codes;

-Split.c: This code splits the particle catalogue to use less RAM and run some computations in parallel;

-Critical.cpp: This code makes the Delaunay triangulation to estimate a local density for each particle and to find the local minima and maxima;

-HVFinder.c: This code grows spheres around the density maxima and minima to find the halos and voids;

-Profiles.c: This code computes the density and velocity profile for each halo and void given by the halo/void finder;

-Stacking.c: This code makes the stack of the profiles compute above for a given binning;

-PowerSpectrum.c: This code computes the cross and auto power spectra for the halos, voids, and matter.

The complete list of input options for each code is given by the code when started without any input.

The codes here needs the library CGAL to work.

If you use any of the codes present here, please cite at least one of the papers: https://arxiv.org/abs/1609.02544, https://arxiv.org/abs/2003.06411.

Please, contact me if you have any question: rodrigo.voivodic@usp.br
