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

void parallel_read_mesh_cgns(idx_t*& elmdist, idx_t*& eptr, idx_t*& eind, idx_t*& epart, idx_t& esize, idx_t& dim, const idx_t numberingstart, std::string filename, MPI_Comm comm);

void read_boundary_conditions(std::vector<partition> &parts, std::string filename);

void write_cgns_separate(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition);

void write_cgns_single(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition);

void read_arrays(std::vector<partition> &parts, std::string filename, MPI_Comm comm);

void write_wo_recv_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim);

void write_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim);

void write_recv_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim);

void write_send_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim);

void read_nodes_cgns(std::vector<partition> &parts, std::string filename, MPI_Comm comm);

void write_pcgns_without_send_recv_info(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition);

void write_pcgns_hybird_with_send_recv_info(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition);

#endif
