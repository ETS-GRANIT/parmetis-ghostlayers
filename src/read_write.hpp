#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "mpi.h"

#include <metis.h>
#include <parmetis.h>

#include "cgnslib.h"
#include "pcgnslib.h"

#include "module_parmetis.hpp"

#ifndef READ_WRITE
#define READ_WRITE

void ParallelReadMeshCGNS(idx_t*& elmdist, idx_t*& eptr, idx_t*& eind, idx_t*& part, idx_t& esize, idx_t& dim, const idx_t numberingstart, std::string filename, MPI_Comm comm);

//ParallelReadMesh refactored from ParallelReadMesh in parmetis/programs/io.c
void ParallelReadMesh(idx_t*& elmdist, idx_t*& eptr, idx_t*& eind, idx_t*& part, idx_t& esize, idx_t& dim, const idx_t numberingstart, std::string filename, MPI_Comm comm);

void boundaryConditionsNodes(std::vector<submesh> &submeshesowned, std::string filename);

void boundaryConditionsCute(std::vector<submesh> &submeshesowned, std::string file);

void writeneighbors(std::vector<submesh> &submeshesowned, idx_t esize);

void writeCute(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void writeMeshCGNS1(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);

void writeMeshCGNS2(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);

void readArrays(std::vector<submesh> &submeshesowned, std::string filename, MPI_Comm comm);

void writeworecvVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void writeVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void writesendrecvCute(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void writerecvVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void writesendVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void writeboundaryVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim);

void updateNodes(std::vector<submesh> &submeshesowned, std::string nodesFilename, MPI_Comm comm);

void updateNodesCGNS(std::vector<submesh> &submeshesowned, std::string filename, MPI_Comm comm);

void writeMeshPCGNS_wos(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);
void writeMeshPCGNS_ch(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);
void writeMeshPCGNS(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);

#endif
