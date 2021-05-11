#include <iostream>
#include <iomanip>
#include <chrono>
#include <assert.h>

#include <metis.h>
#include <parmetis.h>

#include "module_parmetis.hpp"
#include "read_write.hpp"

extern "C" int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb);


#if CGNSVERSION < 3100
# define cgsizet int
#else
# if CGBUILDSCOPE
# error enumeration scoping needs to be off
# endif
#endif

int main(int argc, char *argv[]) {
  int me, nprocs;
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm = MPI_COMM_WORLD;

  //Memory usage variables
  long vmrss, vmsize;

  //User inputs
  assert(argc >= 5);
  idx_t ntotparts = atoi(argv[1]);
  std::string basename = argv[2];
  idx_t method = atoi(argv[3]);
  idx_t numlayers = atoi(argv[4]);

  idx_t numberingstart = 1;

  //ParMetis variables
  idx_t *elmdist, *eptr, *eind, *epart;
  idx_t edgecut(-1), ncon(1), ncommonnodes(2), wgtflag(0), numflag(0);
  idx_t esize, dim; //elemsize (number of nodes) and dimension (2D or 3D) to be read
  real_t *tpwgts, *ubvec;

  tpwgts = new real_t[ntotparts]; //ncon=1
  for(int i=0;i<ntotparts;i++){
    tpwgts[i] = 1./ntotparts;
  }

  ubvec = new real_t[1]; //ncon=1
  ubvec[0] = 1.05;

  //Read the element mesh file and construct all the vectors needed by parmetis
  /* ParallelReadMesh(elmdist, eptr, eind, epart, esize, dim, numberingstart, elemsfile, comm); */
  MPI_Barrier(comm);
  auto tio01 = std::chrono::high_resolution_clock::now();
  parallel_read_mesh_cgns(elmdist, eptr, eind, epart, esize, dim, numberingstart, basename, comm);
  auto tio02 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  assert(dim==2 or dim==3);
  if(dim==2){
    assert(esize==3 or esize==4); //Triangles or quadrangles for 2D (could be extented)
  }
  else if(dim==3){
    assert(esize==4); //Only tet for 3D
  }

  MPI_Barrier(comm);
  idx_t options[METIS_NOPTIONS];
  options[METIS_OPTION_CONTIG] = 1;
  auto t1 = std::chrono::high_resolution_clock::now();
  ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &ntotparts, tpwgts, ubvec, options, &edgecut, epart, &comm);
  auto t2 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);


  idx_t ntotpartsowned;
  std::vector<partition> parts;
  std::vector<idx_t> ownerofpartition;

  delete [] eptr;

  MPI_Barrier(comm);
  auto t3 = std::chrono::high_resolution_clock::now();
  compute_partition_ownership(ntotparts, ntotpartsowned, parts, ownerofpartition, comm);
  gather_partitions(elmdist, eind, epart, esize, parts, ownerofpartition, comm);
  auto t4 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  delete [] eind;
  delete [] epart;

  MPI_Barrier(comm);
  auto tio1 = std::chrono::high_resolution_clock::now();
  read_nodes_cgns(parts, basename, comm);
  auto tio2 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  MPI_Barrier(comm);
  auto t5 = std::chrono::high_resolution_clock::now();
  build_connectivity(parts, dim);
  find_boundary_from_connectivity(parts, method, numlayers);
  compute_potential_neighbors_from_extents(ntotparts, parts, comm);
  share_boundaryFixPotentialNeighbors(parts, ownerofpartition, comm);
  share_boundary(parts, ownerofpartition, comm);
  find_nodes_send_recv(parts, dim, method, numlayers);
  add_elemsAndRenumber(parts);
  auto t6 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  auto tio3 = std::chrono::high_resolution_clock::now();
  read_boundary_conditions(parts, basename);
  read_arrays(parts, basename, comm);
  auto tio4 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  get_memory_usage_kb(&vmrss, &vmsize);

  MPI_Barrier(comm);
  auto tio7 = std::chrono::high_resolution_clock::now();
  write_cgns_separate(parts, esize, dim, ownerofpartition);
  auto tio8 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  MPI_Barrier(comm);
  auto tio5 = std::chrono::high_resolution_clock::now();
  /* write_cgns_single(parts, esize, dim, ownerofpartition); */
  auto tio6 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  MPI_Barrier(comm);
  auto tio9 = std::chrono::high_resolution_clock::now();
  write_pcgns_without_send_recv_info(parts, esize, dim, ownerofpartition);
  auto tio10 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  MPI_Barrier(comm);
  auto tio11 = std::chrono::high_resolution_clock::now();
  write_pcgns_hybird_with_send_recv_info(parts, esize, dim, ownerofpartition);
  auto tio12 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  double duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()/1.0e6;
  double duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count()/1.0e6;
  double duration3 = std::chrono::duration_cast<std::chrono::microseconds>(t6-t5).count()/1.0e6;
  double duration4 = std::chrono::duration_cast<std::chrono::microseconds>(tio2-tio1).count()/1.0e6;
  double duration5 = std::chrono::duration_cast<std::chrono::microseconds>(tio4-tio3).count()/1.0e6;
  double duration6 = std::chrono::duration_cast<std::chrono::microseconds>(tio6-tio5).count()/1.0e6;
  double duration7 = std::chrono::duration_cast<std::chrono::microseconds>(tio8-tio7).count()/1.0e6;
  double duration8 = std::chrono::duration_cast<std::chrono::microseconds>(tio10-tio9).count()/1.0e6;
  double duration04 = std::chrono::duration_cast<std::chrono::microseconds>(tio02-tio01).count()/1.0e6;
  double duration9 = std::chrono::duration_cast<std::chrono::microseconds>(tio12-tio11).count()/1.0e6;
  
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, &duration1, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration2, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration3, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration4, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration5, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration6, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration7, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration8, 1, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &duration9, 1, MPI_DOUBLE, MPI_SUM, comm);

  MPI_Allreduce(MPI_IN_PLACE, &vmrss, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &vmsize, 1, MPI_LONG, MPI_SUM, comm);

  MPI_Barrier(MPI_COMM_WORLD);

  if(me==0){
    std::cout << " Info : ParMETIS GLAS Reading CGNS_Single  CGNS_Multi_files PCGNS_Single Hybrid" << std::endl;
    std::cout << std::setfill(' ') << std::setw(5) << "Temps : " << nprocs << " " << duration1/nprocs << " " << (duration2+duration3)/nprocs << " " << (duration04+duration4+duration5)/nprocs << " " << duration6/nprocs << " " << duration7/nprocs <<  " "  << duration8/nprocs << " " << duration9/nprocs << std::endl;
    std::cout << "Memory usage : " << nprocs << " " << vmrss << " " << vmrss/nprocs << std::endl;
  }

  MPI_Finalize();
}
