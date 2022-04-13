# parmetis-ghostlayers

The source files contained in this repository are provied in the context of the publication of a research paper entitled "Domain decomposition with multiple layers of ghost cells using the CFD General Notation System : Application to the resolution of the Shallow-water equations using a multi-GPU solver" by Delmas, Vincent and Soulaimani, Azzeddine from *Ecole de Technologie Supérieure*, Montréal in the journal *Computer Communications in Physics*.

A mockup of a "
CMakeLists.txt" is provided for compilation with CMake. The required libraries with the specific version used are :
- openmpi/4.0.3
- metis/5.1.0
- parmetis/4.0.3
- hdf5-mpi/1.10.6
- cgns/4.1.2

It has to be noted that the CGNS library needs to be compiled with parallel support (PCGNS). More information on this library can be found on https://github.com/CGNS/
