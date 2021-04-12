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


#ifndef MODULE_PARMETIS_GHOST
#define MODULE_PARMETIS_GHOST

// Simple pair hashing function allow to use unordered set of int pair
struct pair_idx_t_hash { 
  size_t operator()(const std::pair<idx_t, idx_t>& p) const
  { 
    return(p.first);
  } 
}; 

// Simple vector hashing function allow to use unordered set of vector<idx_t>
struct vector_idx_t_hash { 
  size_t operator()(const std::vector<idx_t>& p) const
  { 
    return(p[0]);
  } 
}; 


// Submesh struct so store a sub-domain
struct submesh{
  int submeshid; // Number of the sub-domain
  /* idx_t nelems, nnodes; */ 
  int esize; // Element size (number of vertices of the elements)

  std::vector<std::vector<real_t> > extents; // Spacial extents
  std::set<idx_t> potentialneighbors; // Set of potential neighbors id

  // Main 3 vectors of the sub-domain
  std::vector<real_t> nodes;
  std::vector<idx_t> elems;
  std::vector<idx_t> neighbors;

  idx_t get_nnodes(){return nodes.size()/3;};
  idx_t get_nelems(){return elems.size()/esize;};

  real_t& get_nodes(idx_t i, idx_t j){return nodes[i*3+j];};
  idx_t& get_elems(idx_t i, idx_t j){return elems[i*esize+j];};
  idx_t& get_neighbors(idx_t i, idx_t j){return neighbors[i*esize+j];};

  // Sets of boundary elems and nodes
  std::set<idx_t> boundaryelems;
  std::set<idx_t> boundarynodes;
  std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> boundary_edges;
  std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> boundary_faces;

  // Renumbering vectors of the elements
  std::vector<idx_t> renumber_otn; // otn = old_to_new
  std::vector<idx_t> renumber_nto; // nto = new_to_old

  // Numbering maps for nodes and elements
  std::map<idx_t,idx_t> elems_gtl;
  std::map<idx_t,idx_t> elems_ltg;
  std::map<idx_t,idx_t> nodes_gtl;
  std::map<idx_t,idx_t> nodes_ltg;
  std::unordered_set<idx_t> nodesindex; 

  // Set of elements to send and recv from each potential neighbour
  std::map<idx_t, std::set<idx_t> > elemstosend;
  std::map<idx_t, std::set<idx_t> > elemstorecv;
  std::map<idx_t, std::set<idx_t> > nodestosend;
  std::map<idx_t, std::set<idx_t> > nodestorecv;

  // Set of nodes for boundary conditions
  std::vector<std::set<idx_t> > boundary_conditions;
  std::vector<std::string> boundary_conditions_names;
  std::vector<CGNS_ENUMV(BCType_t)> boundary_conditions_types;

  // Userdefined arrays in cgnsfile
  std::vector<std::string>  ud_names;
  std::vector<std::vector<std::string> > ar_names;
  std::vector<std::vector<std::vector<double> > > arrays;



  // Function to add a set of elements to the sub-domain
  void Addelems(idx_t *elems, idx_t offset, idx_t nelems, idx_t esize);

  // Function to know if an element is part of the boundary of the sub-domain
  bool isBoundary(idx_t ielem);

  // Function to compute the sub-domain spacial extents
  void computemyextents();
};

void Computesubmeshownership(idx_t nsubmeshes, idx_t &nsubmeshesowned, std::vector<submesh> &submeshesowned, std::vector<idx_t> &ownerofsubmesh, MPI_Comm comm);

void Gathersubmeshes(idx_t*& elmdist, idx_t*& eind, idx_t*& part, const idx_t esize, std::vector<submesh> &submeshesowned, std::vector<idx_t> const ownerofsubmesh, MPI_Comm comm);

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

void buildBoundaryEdges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::set<idx_t> &elems_set);

void buildEdges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges);

void buildBoundaryFaces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::set<idx_t> &elems_set);

void buildFaces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces);

void buildElemConnectivity2D(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::vector<idx_t> &neighbors);

void buildElemConnectivity3D(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::vector<idx_t> &neighbors);

void Buildconnectivity(std::vector<submesh> &submeshesowned, idx_t dimension);

void Findboundaryfromconnectivity(std::vector<submesh> &submeshesowned, idx_t method, idx_t numlayers);

void FindNodesElemsSendRecv(std::vector<submesh> &submeshesowned, idx_t dimension, idx_t method, idx_t numlayers);

void Shareboundary(std::vector<submesh> &submeshesowned, std::vector<idx_t> &ownerofsubmesh, MPI_Comm comm);

void Computepotentialneighbors(idx_t nsubmeshes, std::vector<submesh> &submeshesowned, MPI_Comm comm);

void AddElemsAndRenumber(std::vector<submesh> &submeshesowned);

void writeMeshPCGNS_wos(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);
void writeMeshPCGNS_ch(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);
void writeMeshPCGNS(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh);

#endif
