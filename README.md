# parmetis-ghostlayers

The source files contained in this repository are provied in the context of the publication of a research paper entitled "Domain decomposition with multiple layers of ghost cells using the CFD General Notation System : Application to the resolution of the Shallow-water equations using a multi-GPU solver" by Delmas, Vincent and Soulaimani, Azzeddine from *Ecole de Technologie Supérieure*, Montréal in the journal *Computer Physics Communications*.

A mockup of a "
CMakeLists.txt" is provided for compilation with CMake. The required libraries with the specific versions used are :
- openmpi/4.0.3
- metis/5.1.0
- parmetis/4.0.3
- hdf5-mpi/1.10.6
- cgns/4.1.2

It should be noted that the CGNS library needs to be compiled with parallel support (PCGNS). More information on this library can be found on https://github.com/CGNS/

Once the code has been compiled it can be launched using mpirun to partition a mesh using ParMETIS and add as many layers of ghost cells as demanded with

mpirun -n NP main NSD MESH_FILE.cgns METHOD NLAYERS

where, NP is the number of MPI processes (>=2), NSD is the number of sub-domains wanted, MESH_FILE.cgns is the mesh file, METHOD is either 0 to look for adjacent cells with an edge in common and 1 to look for adjacent cells with a vertex in common (finite volumes is fine with 0), and NLAYERS is the number of layer of ghost cells wanted.

A typical output is (here for the command  mpirun -n 2 ./main 16 mille_48.cgns 0 2) :

Info : Np ParMETIS GLAS Reading CGNS_Single PCGNS_Single PCGNS_Hybrid      
Times (s) : 2 0.388324 0.55677 0.517044 0 0 0.657898        
Memory usage (MB) : 2 320016 160008         
        

Where the times for different sections of the code are reported as well as the memory usage (first is total memory usage (RSS) and second is the average memory usage per process). In the last version of the code the CGNS_Single and PCGNS_Single output have been deactivated resulting in times of 0s. The only mesh file produced by the program is named "Mesh_Output_pcgns_ch.cgns" and contains the whole domain decomposition and the send/recv information as desribed in the paper.
