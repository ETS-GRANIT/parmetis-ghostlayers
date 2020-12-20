#include <iostream>
#include <parmetis.h>
#include <metis.h>
#include <chrono>
#include <assert.h>

#include "module_parmetis.hpp"

int main(int argc, char *argv[]) {
  int me, nprocs;
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm = MPI_COMM_WORLD;

  //User inputs
  assert(argc >= 3);
  idx_t nsubmeshes = atoi(argv[1]);
  std::string elemsfile = argv[2];
  std::string nodesfile = argv[3];

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
  ParallelReadMesh(elmdist, eptr, eind, part, esize, dim, numberingstart, elemsfile, comm);

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

  idx_t nsubmeshesowned;
  std::vector<submesh> submeshesowned;
  std::vector<idx_t> ownerofsubmesh;

  auto t3 = std::chrono::high_resolution_clock::now();
  Computesubmeshownership(nsubmeshes, nsubmeshesowned, submeshesowned, ownerofsubmesh, comm);
  Gathersubmeshes(elmdist, eind, part, esize, submeshesowned, ownerofsubmesh, comm);
  auto t4 = std::chrono::high_resolution_clock::now();

  //We do not need eind, eptr, part anymore
  if(eind) delete [] eind;
  if(eptr) delete [] eptr;
  if(part) delete [] part;

  /* updateNodes(submeshesowned, nodesfile, comm); */

  auto t5 = std::chrono::high_resolution_clock::now();
  Buildconnectivity(submeshesowned, dim);
  Findboundaryfromconnectivity(submeshesowned);
  auto t6 = std::chrono::high_resolution_clock::now();

  /* writeVTK(submeshesowned, esize, dim); */
  /* writeboundaryVTK(submeshesowned, esize, dim); */

  auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
  auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count();
  auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(t6-t5).count();
  std::cout << me << " "  << duration1/100000 << " " << duration2/100000 << " " << duration3/100000 << "\n";

  MPI_Finalize();
}
