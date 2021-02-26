#include "cgnslib.h"
#include <iostream>
#include <parmetis.h>
#include <metis.h>
#include <chrono>
#include <assert.h>
#include <unistd.h> 
#include <sstream>


#include "module_parmetis.hpp"

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

  //User inputs
  assert(argc >= 5);
  idx_t nsubmeshes = atoi(argv[1]);
  std::string basename = argv[2];
  idx_t method = atoi(argv[3]);
  idx_t numlayers = atoi(argv[4]);

  idx_t numberingstart = 1;

  //ParMetis variables
  idx_t *elmdist, *eptr, *eind, *part;
  idx_t edgecut(-1), ncon(1), ncommonnodes(2), wgtflag(0), numflag(0), options(0);
  idx_t esize, dim; //elemsize (number of nodes) and dimension (2D or 3D) to be read
  real_t *tpwgts, *ubvec;

  tpwgts = new real_t[nsubmeshes]; //ncon=1
  for(int i=0;i<nsubmeshes;i++){
    tpwgts[i] = 1./nsubmeshes;
  }

  ubvec = new real_t[1]; //ncon=1
  ubvec[0] = 1.05;

  //Read the element mesh file and construct all the vectors needed by parmetis
  /* ParallelReadMesh(elmdist, eptr, eind, part, esize, dim, numberingstart, elemsfile, comm); */
  auto tio01 = std::chrono::high_resolution_clock::now();
  ParallelReadMeshCGNS(elmdist, eptr, eind, part, esize, dim, numberingstart, basename, comm);
  auto tio02 = std::chrono::high_resolution_clock::now();

  assert(dim==2 or dim==3);
  if(dim==2){
    assert(esize==3 or esize==4); //Triangles or quadrangles for 2D (could be extented)
  }
  else if(dim==3){
    assert(esize==4); //Only tet for 3D
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nsubmeshes, tpwgts, ubvec, &options, &edgecut, part, &comm);
  auto t2 = std::chrono::high_resolution_clock::now();

  MPI_Barrier(comm);

  idx_t nsubmeshesowned;
  std::vector<submesh> submeshesowned;
  std::vector<idx_t> ownerofsubmesh;

  if(eptr) delete [] eptr;

  auto t3 = std::chrono::high_resolution_clock::now();
  Computesubmeshownership(nsubmeshes, nsubmeshesowned, submeshesowned, ownerofsubmesh, comm);
  Gathersubmeshes(elmdist, eind, part, esize, submeshesowned, ownerofsubmesh, comm);
  auto t4 = std::chrono::high_resolution_clock::now();

  MPI_Barrier(comm);
  if(eind) delete [] eind;
  if(part) delete [] part;

  auto tio1 = std::chrono::high_resolution_clock::now();
  updateNodesCGNS(submeshesowned, basename, comm);
  auto tio2 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  auto t5 = std::chrono::high_resolution_clock::now();
  Buildconnectivity(submeshesowned, dim);
  Findboundaryfromconnectivity(submeshesowned, method, numlayers);
  Computepotentialneighbors(nsubmeshes, submeshesowned, comm);
  Shareboundary(submeshesowned, ownerofsubmesh, comm);
  FindNodesElemsSendRecv(submeshesowned, dim, method, numlayers);
  AddElemsAndRenumber(submeshesowned);
  auto t6 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  auto tio3 = std::chrono::high_resolution_clock::now();
  boundaryConditionsNodes(submeshesowned, basename);
  readArrays(submeshesowned, basename, comm);
  auto tio4 = std::chrono::high_resolution_clock::now();
  MPI_Barrier(comm);

  auto tio5 = std::chrono::high_resolution_clock::now();
  /* writeMeshCGNS1(submeshesowned, esize, dim, ownerofsubmesh); */
  writeMeshCGNS2(submeshesowned, esize, dim, ownerofsubmesh);
  auto tio6 = std::chrono::high_resolution_clock::now();

  /* writeMeshCGNS(submeshesowned, esize, dim, ownerofsubmesh); */
  /* writeVTK(submeshesowned, esize, dim); */
  /* writeCute(submeshesowned, esize, dim); */
  /* writesendrecvCute(submeshesowned, esize, dim); */
  /* writesendVTK(submeshesowned, esize, dim); */
  /* writerecvVTK(submeshesowned, esize, dim); */
  /* writeneighbors(submeshesowned, esize); */

  auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
  auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count();
  auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(t6-t5).count();
  auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(tio2-tio1).count();
  auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(tio4-tio3).count();
  auto duration6 = std::chrono::duration_cast<std::chrono::microseconds>(tio6-tio5).count();
  auto duration04 = std::chrono::duration_cast<std::chrono::microseconds>(tio02-tio01).count();
  std::cout << std::setfill(' ') << std::setw(5) << me << "   ParMetis : "  << duration1/1.0e6 << "   GhostLayers : " << (duration2+duration3)/1.0e6 << "   FS Read : " << (duration04+duration4+duration5)/1.0e6 << "   FS Write : " << duration6/1.0e6 << "\n";

  MPI_Finalize();
}
