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

struct potential_neighbors_boundary{
  int esize;

  std::vector<real_t> nodes;
  std::vector<idx_t> elems;

  idx_t get_nnodes(){return nodes.size()/3;};
  idx_t get_nelems(){return elems.size()/esize;};

  real_t& get_nodes(idx_t i, idx_t j){return nodes[i*3+j];};
  idx_t& get_elems(idx_t i, idx_t j){return elems[i*esize+j];};

  //Numbering tables : gtl = global to local, ltg = local to global
  std::map<idx_t,idx_t> elems_gtl;
  std::map<idx_t,idx_t> elems_ltg;
  std::map<idx_t,idx_t> nodes_gtl;
  std::map<idx_t,idx_t> nodes_ltg;

  //Edges/faces and neighbors of the boundary
  std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> edges;
  std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> faces;
  std::vector<idx_t> neighbors;

  bool already_computed(){if(neighbors.size() > 0){return 1;}else{return 0;}};
};

//Global variable of potential neighbors inner boundary to reduce memory usage
extern std::map<idx_t, potential_neighbors_boundary> g_potentialneighbors;

// partition struct so store a sub-domain
struct partition {
  int partitionid;

  int esize; // Element size (number of vertices of the elements)

  std::vector<std::vector<real_t> > extents; // Spacial extents
  std::set<idx_t> potentialneighbors_extents; // Set of extents computed potential neighbors id
  std::set<idx_t> potentialneighbors; // Set of potential neighbors id fixed by looking for same nodes

  // Main 3 vectors of the partition
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
  std::set<idx_t> l1_boundarynodes;
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

  // Userdefined arrays in cgnsfile
  std::vector<std::string>  ud_names;
  std::vector<std::vector<std::string> > ar_names;
  std::vector<std::vector<std::vector<double> > > arrays;

  // Function to add a set of elements to the sub-domain
  void add_elems(idx_t *elems, idx_t offset, idx_t nelems, idx_t esize);

  // Function to know if an element is epart of the boundary of the sub-domain
  bool is_boundary(idx_t ielem);

  // Function to compute the sub-domain spacial extents
  void compute_extents();
};

void compute_partition_ownership(idx_t ntotparts, idx_t &ntotpartsowned, std::vector<partition> &parts, std::vector<idx_t> &ownerofpartition, MPI_Comm comm);

void gather_partitions(idx_t*& elmdist, idx_t*& eind, idx_t*& epart, const idx_t esize, std::vector<partition> &parts, std::vector<idx_t> const ownerofpartition, MPI_Comm comm);

void build_boundary_edges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::set<idx_t> &elems_set);

void build_edges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges);

void buildBoundaryFaces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::set<idx_t> &elems_set);

void build_faces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces);

void build_elems_connectivity_2d(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::vector<idx_t> &neighbors);

void build_elems_connectivity_3d(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::vector<idx_t> &neighbors);

void build_connectivity(std::vector<partition> &parts, idx_t dimension);

void find_boundary_from_connectivity(std::vector<partition> &parts, idx_t method, idx_t numlayers);

void find_nodes_send_recv(std::vector<partition> &parts, idx_t dimension, idx_t method, idx_t numlayers);

void share_boundary(std::vector<partition> &parts, std::vector<idx_t> &ownerofpartition, MPI_Comm comm);

void share_boundaryFixPotentialNeighbors(std::vector<partition> &parts, std::vector<idx_t> &ownerofpartition, MPI_Comm comm);

void compute_potential_neighbors_from_extents(idx_t ntotparts, std::vector<partition> &parts, MPI_Comm comm);

void add_elemsAndRenumber(std::vector<partition> &parts);

#endif
