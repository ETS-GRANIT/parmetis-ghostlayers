#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <parmetis.h>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <assert.h>
#include "pcgnslib.h"

#include "mpi.h"

struct pair_idx_t_hash { 
  size_t operator()(const std::pair<idx_t, idx_t>& p) const
  { 
    return(p.first);
  } 
}; 

struct vector_idx_t_hash { 
  size_t operator()(const std::vector<idx_t>& p) const
  { 
    return(p[0]);
  } 
}; 

struct potential_neighbors_boundary{
  std::vector<std::vector<real_t> > nodes;
  std::vector<std::vector<idx_t> > elems;

  //Numbering of boundary nodes
  std::map<idx_t,idx_t> elems_gtl;
  std::map<idx_t,idx_t> elems_ltg;
  std::map<idx_t,idx_t> nodes_gtl;
  std::map<idx_t,idx_t> nodes_ltg;
};

//Global variable of potential nieghbos to reduce memory usage
std::map<idx_t, potential_neighbors_boundary> g_potentialneighbors;

struct submesh{
  idx_t submeshid;
  idx_t nelems, nnodes, nboundaryelems;
  idx_t esize;

  std::vector<std::vector<real_t> > extents;
  std::set<idx_t> potentialneighbors;

  std::vector<std::vector<idx_t> > elems;
  std::vector<std::vector<idx_t> > neighbors;
  std::vector<std::vector<idx_t> > boundary_neighbors;
  std::set<idx_t> boundaryelems;
  std::set<idx_t> boundarynodes;
  std::vector<std::vector<real_t> > nodes;
  std::vector<idx_t> erenumber;
  std::vector<idx_t> erenumberr;

  //Numbering maps for nodes and elements
  std::map<idx_t,idx_t> elems_gtl;
  std::map<idx_t,idx_t> elems_ltg;
  std::map<idx_t,idx_t> nodes_gtl;
  std::map<idx_t,idx_t> nodes_ltg;
  std::unordered_set<idx_t> nodesindex; 

  //Set of elements to send and recv from each potential neighbour
  std::map<idx_t, std::set<idx_t> > elemstosend;
  std::map<idx_t, std::set<idx_t> > elemstorecv;
  std::map<idx_t, std::set<idx_t> > nodestosend;
  std::map<idx_t, std::set<idx_t> > nodestorecv;

  //Edges and faces
  std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> boundary_edges;
  std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> boundary_faces;

  //Set of nodes for boundary conditions
  std::vector<std::set<idx_t> > boundary_conditions;
  std::vector<std::string> boundary_conditions_names;
  std::vector<CGNS_ENUMV(BCType_t)> boundary_conditions_types;

  //Userdefines arrays in cgnsfile
  std::vector<std::string>  ud_names;
  std::vector<std::vector<std::string> > ar_names;
  std::vector<std::vector<std::vector<double> > > arrays;

  void Addelems(idx_t *elems, idx_t offset, idx_t nelems, idx_t esize);
  bool isBoundary(idx_t ielem);
  void computemyextents(idx_t nsubmeshes);
};

bool submesh::isBoundary(idx_t ielem){
  for(idx_t i=0; i<esize; i++){
    if(neighbors[ielem][i] == -1) return(1);
  }
  return(0);
}


void submesh::Addelems(idx_t *newelems, idx_t offset, idx_t nnewelems, idx_t newesize){
  nelems += nnewelems;
  esize = newesize; //should'nt change
  std::vector<std::vector<idx_t> > tmpelems(nnewelems, std::vector<idx_t>(esize, 0));
  for(idx_t i=0;i<nnewelems;i++){
    elems_gtl.insert(std::make_pair(newelems[offset+i*(esize+1)],i+nelems-nnewelems));
    elems_ltg.insert(std::make_pair(i+nelems-nnewelems,newelems[offset+i*(esize+1)]));

    for(idx_t k=0; k<esize; k++){
      tmpelems[i][k] = newelems[offset+i*(esize+1)+k+1];
      nodesindex.insert(tmpelems[i][k]);
    }
  }
  elems.insert(elems.end(), tmpelems.begin(), tmpelems.end());
}

void Computesubmeshownership(idx_t nsubmeshes, idx_t &nsubmeshesowned, std::vector<submesh> &submeshesowned, std::vector<idx_t> &ownerofsubmesh, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  //Compute the list of partition to ba handled by me
  idx_t firstsubmeshid;
  idx_t q=nsubmeshes/nprocs;
  idx_t r=nsubmeshes-q*nprocs;

  ownerofsubmesh.resize(nsubmeshes);

  idx_t nsubmeshesofp;
  for(idx_t p=0; p<nprocs; p++){
    if(p<r){
      nsubmeshesofp = q+1;
      firstsubmeshid = p*(q+1);
    }
    else{
      nsubmeshesofp = q;
      firstsubmeshid = p*q + r;
    }

    for(idx_t k=0;k<nsubmeshesofp;k++){
      ownerofsubmesh[k+firstsubmeshid] = p;
    }
    if(me==p){
      nsubmeshesowned = nsubmeshesofp;
      submeshesowned.resize(nsubmeshesowned);
      for(idx_t p=0;p<nsubmeshesowned;p++){
        submeshesowned[p].submeshid = p+firstsubmeshid; 
      }
    }
  }
}

void Gathersubmeshes(idx_t*& elmdist, idx_t*& eind, idx_t*& part, const idx_t esize, std::vector<submesh> &submeshesowned, std::vector<idx_t> const ownerofsubmesh, MPI_Comm comm){

  std::stringstream line;
  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;
  MPI_Request *requestSendPtr, *requestRecvPtr;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  idx_t nelems = elmdist[me+1]-elmdist[me];

  //Generate submeshes
  std::map<idx_t, idx_t> iter;
  std::map<idx_t, idx_t*> submesheselems; //Needs to be properly freed
  std::map<idx_t, idx_t> submeshessizes;

  for(idx_t i=0;i<nelems;i++){
    if(submeshessizes.find(part[i]) == submeshessizes.end()){
      submeshessizes.insert(std::make_pair(part[i], 0));
    }
    submeshessizes[part[i]] += 1;
  }

  for(idx_t i=0;i<nelems;i++){
    if(submesheselems.find(part[i]) == submesheselems.end()){
      submesheselems.insert(std::make_pair(part[i], new idx_t[((esize+1)*submeshessizes[part[i]])]));
      iter.insert(std::make_pair(part[i],0));
    }
    submesheselems[part[i]][iter[part[i]]*(esize+1)] = elmdist[me]+i;
    for(idx_t j=0;j<esize;j++){
      submesheselems[part[i]][iter[part[i]]*(esize+1)+j+1] = eind[esize*i+j];
    }
    iter[part[i]] += 1;
  }

  //Communicate the size of each submesh owned
  idx_t ntotalsubmeshes = ownerofsubmesh.size();
  idx_t *gsubmeshessizes = new idx_t[(ntotalsubmeshes)*nprocs];

  for(idx_t i=0;i<ntotalsubmeshes*nprocs;i++){
    gsubmeshessizes[i] = 0;
  }

  for(idx_t i=0;i<ntotalsubmeshes;i++){
    if(submesheselems.find(i) == submesheselems.end()){
      gsubmeshessizes[me*ntotalsubmeshes+i] = 0;
    }
    else{
      gsubmeshessizes[me*ntotalsubmeshes+i] = submeshessizes[i];
    }
  } 

  MPI_Allgather((void *)(gsubmeshessizes+me*ntotalsubmeshes), ntotalsubmeshes, IDX_T, (void *)gsubmeshessizes, ntotalsubmeshes, IDX_T, comm);

  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  idx_t nSendRequest(0);
  for(idx_t i=0;i<ntotalsubmeshes;i++){
    if(gsubmeshessizes[me*ntotalsubmeshes+i] > 0){
      if(ownerofsubmesh[i]!=me){
        nSendRequest+=1;
      }
    }
  } 

  idx_t nRecvRequest(0), maxRecv(0);
  for(idx_t i=0;i<ntotalsubmeshes;i++){
    if(ownerofsubmesh[i]==me){
      for(idx_t p=0;p<nprocs;p++){
        if(gsubmeshessizes[p*ntotalsubmeshes+i]>0 and p!=me){
          nRecvRequest+=1;
          maxRecv=std::max((esize+1)*gsubmeshessizes[p*ntotalsubmeshes+i],maxRecv);
        }
      }
    }
  }

  requestSendPtr = new MPI_Request[nSendRequest];
  requestRecvPtr = new MPI_Request[nRecvRequest];
  statSend = new MPI_Status[nSendRequest];
  statRecv = new MPI_Status[nRecvRequest];

  idx_t k=0;
  for(idx_t i=0;i<ntotalsubmeshes;i++){
    if(gsubmeshessizes[me*ntotalsubmeshes+i] > 0){
      if(ownerofsubmesh[i]!=me){//Send to ownerofsubmesh[i]
        MPI_Isend((void *)submesheselems[i], (esize+1)*submeshessizes[i], IDX_T, ownerofsubmesh[i], me*ntotalsubmeshes+i, comm, &requestSendPtr[k]);
        k+=1;
      }
      else{//Add to my partition
        submeshesowned[submeshlocid[i]].Addelems(submesheselems[i], 0, submeshessizes[i], esize);
      }
    }
  } 

  k=0;
  idx_t *tmpElems = new idx_t[(esize+1)*nRecvRequest*maxRecv];
  for(idx_t i=0;i<ntotalsubmeshes;i++){
    if(ownerofsubmesh[i]==me){
      for(idx_t p=0;p<nprocs;p++){
        if(gsubmeshessizes[p*ntotalsubmeshes+i]>0 and p!=me){//Recv from p and add to my partition
          MPI_Irecv((void *)(tmpElems+k*maxRecv*(esize+1)), (esize+1)*gsubmeshessizes[p*ntotalsubmeshes+i], IDX_T, p, p*ntotalsubmeshes+i, comm, &requestRecvPtr[k]);
          k+=1;
        }
      }
    }
  }

  MPI_Waitall(nSendRequest, requestSendPtr, statSend);
  MPI_Waitall(nRecvRequest, requestRecvPtr, statRecv);

  k=0;
  for(idx_t i=0;i<ntotalsubmeshes;i++){
    if(ownerofsubmesh[i]==me){
      for(idx_t p=0;p<nprocs;p++){
        if(gsubmeshessizes[p*ntotalsubmeshes+i]>0 and p!=me){//Recv from p and add to my partition
          submeshesowned[submeshlocid[i]].Addelems(tmpElems, k*maxRecv*(esize+1), gsubmeshessizes[p*ntotalsubmeshes+i], esize);
          k+=1;
        }
      }
    }
  }

  //Free of submesheselems
  std::map<idx_t, idx_t*>::iterator it;
  for(it=submesheselems.begin();it!=submesheselems.end();it++){
    delete [] it->second;
  }
}

void ParallelReadMeshCGNS(idx_t*& elmdist, idx_t*& eptr, idx_t*& eind, idx_t*& part, idx_t& esize, idx_t& dim, const idx_t numberingstart, std::string filename, MPI_Comm comm){

  int me, nprocs;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  int index_file, index_base, n_bases, base, physDim, cellDim, nZones, zone;
  int gnelems, nelems;
  int nSections;
  cgsize_t gnnodes;
  int nCoords;
  char basename[40];
  char zonename[40];
  char name[40];
  char secname[40];
  cgsize_t sizes[2];
  CGNS_ENUMV(ZoneType_t) zoneType;
  CGNS_ENUMV(DataType_t) dataType;
  CGNS_ENUMV(ElementType_t) type;

  if (cg_open(filename.c_str(),CG_MODE_READ,&index_file)) cg_error_exit();
  if(cg_nbases(index_file, &n_bases)!= CG_OK) cg_get_error();
  if(n_bases != 1) cg_get_error(); 
  base=1;
  if(cg_base_read(index_file, base, basename, &cellDim, &physDim) != CG_OK) cg_get_error();
  if(cg_nzones(index_file, base, &nZones) != CG_OK) cg_get_error();
  zone = 1;
  if(cg_zone_type(index_file, base, zone, &zoneType) != CG_OK) cg_get_error();
  assert(zoneType == CGNS_ENUMV(Unstructured));
  if(cg_zone_read(index_file, base, zone, zonename, sizes) != CG_OK) cg_get_error();
  gnnodes = sizes[0];
  gnelems = sizes[1];

  //Create elmdist
  elmdist = new idx_t[nprocs+1];
  elmdist[0] = 0;
  idx_t j=gnelems;idx_t k;
  for(idx_t i=0; i<nprocs; i++){
    k = j/(nprocs-i);
    elmdist[i+1] = elmdist[i]+k;
    j -= k;
  }
  nelems = elmdist[me+1]-elmdist[me];

  if(cg_nsections(index_file, base, zone, &nSections) != CG_OK) cg_get_error();
  int sec = 1;
  int nBdry;
  cgsize_t eBeg, eEnd, *conn;
  int parentFlag;
  if(cg_section_read(index_file, base, zone, sec, secname, &type,
        &eBeg, &eEnd, &nBdry, &parentFlag) != CG_OK) cg_get_error();
  switch (type)
  {
    case CGNS_ENUMV(TETRA_4):
      esize = 4;
      dim = 3;
      conn = new cgsize_t[nelems*esize]; break;
    case CGNS_ENUMV(TRI_3):
      esize = 3;
      dim = 2;
      conn = new cgsize_t[nelems*esize]; break;
    case CGNS_ENUMV(QUAD_4):
      esize = 4;
      dim = 2;
      conn = new cgsize_t[nelems*esize]; break;
    default:
      exit(0);break;
  }


  // Conectivity starts numbering at 1 !!!
  eBeg=elmdist[me]+1;
  eEnd=elmdist[me+1];
  if(cg_elements_partial_read(index_file, base, zone, sec, eBeg, eEnd, conn, NULL) != CG_OK) cg_get_error();
  eind = new idx_t[esize*nelems];
  part = new idx_t[nelems];
  eptr = new idx_t[nelems+1];

  for(idx_t i=0;i<nelems;i++){
    for(idx_t j=0;j<esize;j++){
      /* eind[i*esize+j] = (idx_t) (conn[i*esize+j]-1); */
      eind[i*esize+j] = (idx_t) (conn[i*esize+j]-1);
    }
  }
  for(idx_t i=0;i<=nelems;i++){
    eptr[i] = esize*i;
  }

  if (cg_close(index_file)) cg_error_exit();

}

//ParallelReadMesh refactored from ParallelReadMesh in parmetis/programs/io.c
void ParallelReadMesh(idx_t*& elmdist, idx_t*& eptr, idx_t*& eind, idx_t*& part, idx_t& esize, idx_t& dim, const idx_t numberingstart, std::string filename, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  elmdist = new idx_t[nprocs+1];

  // Construction of the elmdist vector (see parmetis manual)
  std::stringstream fFilename;
  fFilename << filename;
  std::ifstream meshFile;

  if(me==(nprocs-1)){
    meshFile.open(fFilename.str());
    idx_t gnelems;
    meshFile >> gnelems >> dim >> esize;

    elmdist[0] = 0;
    idx_t j=gnelems;idx_t k;
    for(idx_t i=0; i<nprocs; i++){
      k = j/(nprocs-i);
      elmdist[i+1] = elmdist[i]+k;
      j -= k;
    }
    MPI_Bcast((void *)elmdist, nprocs+1, IDX_T, nprocs-1, comm);
    MPI_Bcast(&esize, 1, IDX_T, nprocs-1, comm);
    MPI_Bcast(&dim, 1, IDX_T, nprocs-1, comm);
  }
  else{
    MPI_Bcast((void *)elmdist, nprocs+1, IDX_T, nprocs-1, comm);
    MPI_Bcast(&esize, 1, IDX_T, nprocs-1, comm);
    MPI_Bcast(&dim, 1, IDX_T, nprocs-1, comm);
  }

  idx_t nelems = elmdist[me+1]-elmdist[me];

  eind = new idx_t[esize*nelems];
  part = new idx_t[nelems];
  eptr = new idx_t[nelems+1];

  if(me==(nprocs-1)){
    idx_t maxnelems = 0;
    for(idx_t i=0;i<nprocs;i++){
      if(maxnelems < (elmdist[i+1]-elmdist[i])){
        maxnelems = elmdist[i+1]-elmdist[i];
      }
    }

    idx_t *tmpElems = new idx_t[esize*maxnelems];

    for(idx_t proc=0; proc<nprocs; proc++){
      idx_t nTmpElems=(elmdist[proc+1]-elmdist[proc]);
      for(idx_t i=0;i<nTmpElems;i++){
        for(idx_t j=0;j<esize;j++){
          meshFile >> tmpElems[i*esize+j];
        }
      }

      if(proc!=nprocs-1){
        MPI_Send((void *)tmpElems, esize*nTmpElems, IDX_T, proc, 0, comm);
      }
      else{
        for(idx_t i=0;i<nelems;i++){
          for(idx_t j=0;j<esize;j++){
            eind[i*esize+j] = tmpElems[i*esize+j];
          }
        }
      }
    }
    delete [] tmpElems;
  }
  else{
    MPI_Recv((void *)eind, esize*nelems, IDX_T, nprocs-1, 0, comm, &stat);
  }

  //Correction of the numbering according to numberingstart
  if(numberingstart != 0){
    for(idx_t i=0;i<=nelems*esize;i++){
      eind[i] -= numberingstart;
    }
  }

  for(idx_t i=0;i<=nelems;i++){
    eptr[i] = esize*i;
  }
  meshFile.close();
}

void boundaryConditionsNodes(std::vector<submesh> &submeshesowned, std::string filename){
  int index_file, index_base, n_bases, base, physDim, cellDim, nZones, zone, nBocos;
  int normalIndex, nDataSet;
  int gnelems, nelems;
  int nSections;
  cgsize_t gnnodes;
  cgsize_t normListFlag;
  int nCoords;
  char basename[40];
  char zonename[40];
  char name[40];
  char secname[40];
  char boconame[40];
  cgsize_t sizes[2];
  CGNS_ENUMV(ZoneType_t) zoneType;
  CGNS_ENUMV(DataType_t) dataType, normDataType;
  CGNS_ENUMV(ElementType_t) type;
  CGNS_ENUMV(BCType_t) bocoType;
  CGNS_ENUMV(PointSetType_t) ptsetType;
  cgsize_t nBCNodes;
  cgsize_t *bcnodes;

  if (cg_open(filename.c_str(),CG_MODE_READ,&index_file)) cg_error_exit();
  if(cg_nbases(index_file, &n_bases)!= CG_OK) cg_get_error();
  if(n_bases != 1) cg_get_error(); 
  base=1;
  if(cg_base_read(index_file, base, basename, &cellDim, &physDim) != CG_OK) cg_get_error();
  if(cg_nzones(index_file, base, &nZones) != CG_OK) cg_get_error();
  zone = 1;
  if(cg_zone_type(index_file, base, zone, &zoneType) != CG_OK) cg_get_error();
  assert(zoneType == CGNS_ENUMV(Unstructured));
  if(cg_zone_read(index_file, base, zone, zonename, sizes) != CG_OK) cg_get_error();
  gnnodes = sizes[0];
  gnelems = sizes[1];

  if(cg_nbocos(index_file, base, zone, &nBocos) != CG_OK) cg_get_error();
  for(int k=0; k<submeshesowned.size(); k++){
    submeshesowned[k].boundary_conditions.resize(nBocos);
    submeshesowned[k].boundary_conditions_names.resize(nBocos);
    submeshesowned[k].boundary_conditions_types.resize(nBocos);
  }
  for(int boco=1; boco<=nBocos; boco++)
  {
    if(cg_boco_info(index_file, base, zone, boco, boconame, &bocoType,
          &ptsetType, &nBCNodes, &normalIndex,
          &normListFlag, &normDataType, &nDataSet) != CG_OK) cg_get_error();
    bcnodes = new cgsize_t[nBCNodes];
    CGNS_ENUMV(GridLocation_t) location;
    if(cg_boco_gridlocation_read(index_file, base, zone, boco, &location) != CG_OK) cg_get_error();
    assert(location==CGNS_ENUMV(Vertex));
    if(cg_boco_read(index_file, base, zone, boco, bcnodes,
          NULL) != CG_OK) cg_get_error();

    for(int k=0; k<submeshesowned.size(); k++){
      submeshesowned[k].boundary_conditions_names[boco-1] = boconame;
      submeshesowned[k].boundary_conditions_types[boco-1] = bocoType;
    }
    for(cgsize_t i=0; i<nBCNodes; i++){
      for(int k=0; k<submeshesowned.size(); k++){
        if(submeshesowned[k].nodes_gtl.count(bcnodes[i]) != 0){
          submeshesowned[k].boundary_conditions[boco-1].insert(bcnodes[i]);
        }
      }
    }

    delete [] bcnodes;
  }
}

void boundaryConditionsCute(std::vector<submesh> &submeshesowned, std::string file){
  idx_t nnodesin;
  idx_t *nnodesout = new idx_t[submeshesowned.size()];
  for(idx_t k=0; k<submeshesowned.size(); k++){
    nnodesout[k] = 0;
  }
  idx_t dummy, dummy2;
  std::ifstream inputfile;
  std::ofstream *outputfiles;
  outputfiles = new std::ofstream[submeshesowned.size()];
  std::stringstream name;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    name.str(" ");
    name << submeshesowned[k].submeshid << "_" << file;
    outputfiles[k].open(name.str());
  }

  inputfile.open(file);
  inputfile >> nnodesin;
  for(idx_t i=0; i<nnodesin; i++){
    inputfile >> dummy;
    for(idx_t k=0; k<submeshesowned.size(); k++){
      if(submeshesowned[k].nodes_gtl.find(dummy) != submeshesowned[k].nodes_gtl.end()){
        nnodesout[k]++;
      }
    }
  }
  inputfile.close();

  for(idx_t k=0; k<submeshesowned.size(); k++){
    outputfiles[k] << nnodesout[k] << std::endl;
  }

  inputfile.open(file);
  inputfile >> nnodesin;
  for(idx_t i=0; i<nnodesin; i++){
    inputfile >> dummy;
    for(idx_t k=0; k<submeshesowned.size(); k++){
      if(submeshesowned[k].nodes_gtl.find(dummy) != submeshesowned[k].nodes_gtl.end()){
        outputfiles[k] << submeshesowned[k].nodes_gtl[dummy] << std::endl;
      }
    }
  }
  inputfile.close();
}

void writeneighbors(std::vector<submesh> &submeshesowned, idx_t esize){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_neighbors.dat";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    for(idx_t i=0; i<submeshesowned[k].nelems; i++){
      outFile << submeshesowned[k].elems_ltg[i] << " ";
      for(idx_t j=0; j<submeshesowned[k].esize; j++){
        idx_t value = submeshesowned[k].neighbors[i][j];
        if(value!=-1){
          outFile << submeshesowned[k].elems_ltg[submeshesowned[k].neighbors[i][j]] << " ";
        }
        else{
          outFile << "-1 ";
        }
      }
      outFile << std::endl;
    }
  }
}

void writeCute(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t indrenum;
    idx_t nNodes = submeshesowned[k].nnodes;
    idx_t nElems = submeshesowned[k].nelems;
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_cuteflow.txt";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "Table de coordonnees" << std::endl;
    outFile << nNodes << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << i+1 << " " << submeshesowned[k].nodes[i][0] << " " << submeshesowned[k].nodes[i][1] << " " << submeshesowned[k].nodes[i][2] << std::endl;
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems << " " << nElems*(esize+1) << std::endl;
    outFile << "Table de connectivité" << std::endl;
    outFile << nElems << std::endl;
    for(idx_t i=0; i<nElems; i++){
      outFile << i+1 << " ";
      indrenum = submeshesowned[k].erenumber[i];
      for(idx_t p=0; p<esize; p++){
        outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].elems[indrenum][p]] << " ";
      }
      outFile << std::endl;
    }
  }
}

void writeMeshCGNS1(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm(MPI_COMM_WORLD);

  int index_file, index_base, index_zone, n_bases, base, physDim, cellDim, nZones, zone, index_coord, index_section, index_bc;
  int gnelems, nelems;
  int nSections;
  cgsize_t gnnodes;
  int nCoords;
  const char *zonename;
  char name[40];
  char secname[40];
  cgsize_t sizes[2];
  CGNS_ENUMV(ZoneType_t) zoneType;
  CGNS_ENUMV(DataType_t) dataType;
  CGNS_ENUMV(ElementType_t) elementType;
  double *coord;
  cgsize_t *elems, *elemstosend, *elemstorecv;
  int Cx, Cy, Cz;

  int icelldim;
  int iphysdim;
  if(dim == 2 and esize==3){
    icelldim = 2;
    iphysdim = 3;
    elementType = CGNS_ENUMV(TRI_3);
  }
  else if(dim==2 and esize==4){
    icelldim = 2;
    iphysdim = 3;
    elementType = CGNS_ENUMV(QUAD_4);
  }
  else if(dim==3 and esize==4){
    icelldim = 3;
    iphysdim = 3;
    elementType = CGNS_ENUMV(TETRA_4);
  }

  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  std::vector<std::string> zms(ownerofsubmesh.size());
  for(int k=0; k<ownerofsubmesh.size(); k++){
    std::stringstream zss;
    zss << "River_" << k;
    zms[k] = zss.str();
  }

  cgsize_t cgzones[ownerofsubmesh.size()][4];

  index_base=1;
  for(int k=0; k<submeshesowned.size(); k++){
    std::stringstream fname;
    fname << submeshesowned[k].submeshid << "_Mesh_Output.cgns" << std::endl;
    if(cg_open(fname.str().c_str(),CG_MODE_WRITE,&index_file)) cg_error_exit();
    if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base)) cg_error_exit();
    /* cgsize_t estart, eend; */
    cgsize_t isize[1][3];
    isize[0][0]=submeshesowned[k].nodes.size();
    isize[0][1]=submeshesowned[k].elems.size();
    isize[0][2]=0;
    std::stringstream zss;
    zss << "River_" << submeshesowned[k].submeshid;
    if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&index_zone) ) cg_error_exit();
    coord = new double[submeshesowned[k].nodes.size()];
    for(idx_t i=0; i<submeshesowned[k].nodes.size(); i++){
      coord[i] = submeshesowned[k].nodes[i][0];
    }
    if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateX", coord, &Cx)) cg_error_exit();
    for(idx_t i=0; i<submeshesowned[k].nodes.size(); i++){
      coord[i] = submeshesowned[k].nodes[i][1];
    }
    if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateY", coord, &Cy)) cg_error_exit();
    for(idx_t i=0; i<submeshesowned[k].nodes.size(); i++){
      coord[i] = submeshesowned[k].nodes[i][2];
    }
    if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateZ", coord, &Cz)) cg_error_exit();
    delete [] coord;
    int nstart=1, nend=submeshesowned[k].elems.size();
    elems = new cgsize_t[esize*submeshesowned[k].elems.size()];
    for(idx_t i=0; i<submeshesowned[k].elems.size(); i++){
      for(idx_t j=0; j<esize; j++){
        elems[esize*i + j] = submeshesowned[k].erenumber[submeshesowned[k].nodes_gtl[submeshesowned[k].elems[i][j]]]+1;
      }
    }
    if(cg_section_write(index_file,index_base,index_zone,"Elements",elementType,nstart,nend,0,elems,&index_section)) cg_error_exit();
    delete [] elems;
    std::string name="ElemsToSend";
    if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
    if(cg_user_data_write(name.c_str())) cg_error_exit();
    for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[k].elemstosend.begin(); it!=submeshesowned[k].elemstosend.end(); it++){
      if(it->second.size()>0){
        elemstosend = new cgsize_t[it->second.size()];
        cgsize_t i=0;
        for(std::set<idx_t>::iterator iter=it->second.begin(); iter!=it->second.end(); iter++){
          elemstosend[i] = submeshesowned[k].erenumber[*iter]+1;
          i++;
        }
        std::stringstream ssname;
        ssname << it->first;
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit();
        if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToSend", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
        /* std::cout << it->second.size() << std::endl; */
        if(cg_ptset_write(CGNS_ENUMV(PointList), it->second.size(), elemstosend)) cg_error_exit();
        if(cg_ordinal_write(it->first));
        delete [] elemstosend;
      }
    }
    name="ElemsToRecv";
    if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
    if(cg_user_data_write(name.c_str())) cg_error_exit();
    for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[k].elemstorecv.begin(); it!=submeshesowned[k].elemstorecv.end(); it++){
      if(it->second.size()>0){
        std::stringstream ssname;
        ssname << it->first;
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToRecv", 0, "end")) cg_error_exit();
        if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToRecv", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
        elemstorecv = new cgsize_t[2];
        std::set<idx_t>::iterator iter=it->second.begin();
        idx_t itglobloc = submeshesowned[k].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*iter]];
        elemstorecv[0] = submeshesowned[k].erenumber[itglobloc]+1;
        elemstorecv[1] = it->second.size();
        if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
        if(cg_ordinal_write(it->first));
        delete [] elemstorecv;
      }
    }

    //Write boundary conditions
    for(int bc=0; bc<submeshesowned[k].boundary_conditions.size(); bc++){
      if(submeshesowned[k].boundary_conditions[bc].size()>0){
        cgsize_t bcnodes[submeshesowned[k].boundary_conditions[bc].size()];
        cgsize_t nbc=0;
        for(std::set<idx_t>::iterator it=submeshesowned[k].boundary_conditions[bc].begin();
            it!=submeshesowned[k].boundary_conditions[bc].end();it++){
          bcnodes[nbc] = submeshesowned[k].nodes_gtl[*it] ;
          nbc++;
        }
        if(cg_boco_write(index_file,index_base,index_zone,submeshesowned[k].boundary_conditions_names[bc].c_str(),submeshesowned[k].boundary_conditions_types[bc],CGNS_ENUMV(PointList),submeshesowned[k].boundary_conditions[bc].size(),bcnodes,&index_bc)) cg_error_exit();
        cg_boco_gridlocation_write(index_file,index_base,index_zone,index_bc,CGNS_ENUMV(Vertex));
      }
    }

    int nuserdata = submeshesowned[k].ud_names.size();
    for(int nud=1; nud<=nuserdata; nud++){
      int narrays = submeshesowned[k].ar_names[nud-1].size();
      if(narrays>0){
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
        /* std::cout << submeshesowned[k].submeshid << " " << submeshesowned[k].ud_names[nud-1] << std::endl; */
        if(cg_user_data_write(submeshesowned[k].ud_names[nud-1].c_str())) cg_error_exit();
        cgsize_t dimensions=submeshesowned[k].arrays[nud-1][0].size();
        CGNS_ENUMV(GridLocation_t) location;
        if(dimensions==submeshesowned[k].nodes.size()){
          location = CGNS_ENUMV(Vertex);
        }
        if(dimensions==submeshesowned[k].elems.size()){
          location = CGNS_ENUMV(CellCenter);
        }
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, submeshesowned[k].ud_names[nud-1].c_str(), 0, "end")) cg_error_exit();
        if(cg_gridlocation_write(location)) cg_error_exit();
        for(int na=1; na<=narrays; na++){
          int rank=1;
          if(cg_array_write(submeshesowned[k].ar_names[nud-1][na-1].c_str(), CGNS_ENUMV(RealDouble), 1, &dimensions, submeshesowned[k].arrays[nud-1][na-1].data())) cg_error_exit();

        }
      }

    }
    if(cg_close(index_file)) cg_error_exit();
  }
}

void writeMeshCGNS2(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm(MPI_COMM_WORLD);

  int index_file, index_base, index_zone, n_bases, base, physDim, cellDim, nZones, zone, index_coord, index_section, index_bc;
  int gnelems, nelems;
  int nSections;
  cgsize_t gnnodes;
  int nCoords;
  const char *zonename;
  char name[40];
  char secname[40];
  cgsize_t sizes[2];
  CGNS_ENUMV(ZoneType_t) zoneType;
  CGNS_ENUMV(DataType_t) dataType;
  CGNS_ENUMV(ElementType_t) elementType;
  double *coord;
  cgsize_t *elems, *elemstosend, *elemstorecv;
  int Cx, Cy, Cz;

  int icelldim;
  int iphysdim;
  if(dim == 2 and esize==3){
    icelldim = 2;
    iphysdim = 3;
    elementType = CGNS_ENUMV(TRI_3);
  }
  else if(dim==2 and esize==4){
    icelldim = 2;
    iphysdim = 3;
    elementType = CGNS_ENUMV(QUAD_4);
  }
  else if(dim==3 and esize==4){
    icelldim = 3;
    iphysdim = 3;
    elementType = CGNS_ENUMV(TETRA_4);
  }

  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  std::vector<std::string> zms(ownerofsubmesh.size());
  for(int k=0; k<ownerofsubmesh.size(); k++){
    std::stringstream zss;
    zss << "River_" << k;
    zms[k] = zss.str();
  }

  if(me==0){
    if(cg_open("Mesh_Output.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
    if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base)) cg_error_exit();
    if(cg_close(index_file)) cg_error_exit();
  }
  MPI_Barrier(comm);

  cgsize_t cgzones[ownerofsubmesh.size()][4];

  for(int p=0; p<nprocs; p++){
    if(me==p){
      if(cg_open("Mesh_Output.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();
      index_base=1;
      for(int k=0; k<submeshesowned.size(); k++){
        /* cgsize_t estart, eend; */
        cgsize_t isize[1][3];
        isize[0][0]=submeshesowned[k].nodes.size();
        isize[0][1]=submeshesowned[k].elems.size();
        isize[0][2]=0;
        std::stringstream zss;
        zss << "River_" << submeshesowned[k].submeshid;
        if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&index_zone) ) cg_error_exit();
        coord = new double[submeshesowned[k].nodes.size()];
        for(idx_t i=0; i<submeshesowned[k].nodes.size(); i++){
          coord[i] = submeshesowned[k].nodes[i][0];
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateX", coord, &Cx)) cg_error_exit();
        for(idx_t i=0; i<submeshesowned[k].nodes.size(); i++){
          coord[i] = submeshesowned[k].nodes[i][1];
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateY", coord, &Cy)) cg_error_exit();
        for(idx_t i=0; i<submeshesowned[k].nodes.size(); i++){
          coord[i] = submeshesowned[k].nodes[i][2];
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateZ", coord, &Cz)) cg_error_exit();
        delete [] coord;
        int nstart=1, nend=submeshesowned[k].elems.size();
        elems = new cgsize_t[esize*submeshesowned[k].elems.size()];
        for(idx_t i=0; i<submeshesowned[k].elems.size(); i++){
          for(idx_t j=0; j<esize; j++){
            elems[esize*i + j] = submeshesowned[k].erenumber[submeshesowned[k].nodes_gtl[submeshesowned[k].elems[i][j]]]+1;
          }
        }
        if(cg_section_write(index_file,index_base,index_zone,"Elements",elementType,nstart,nend,0,elems,&index_section)) cg_error_exit();
        delete [] elems;
        std::string name="ElemsToSend";
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_user_data_write(name.c_str())) cg_error_exit();
        for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[k].elemstosend.begin(); it!=submeshesowned[k].elemstosend.end(); it++){
          if(it->second.size()>0){
            elemstosend = new cgsize_t[it->second.size()];
            cgsize_t i=0;
            for(std::set<idx_t>::iterator iter=it->second.begin(); iter!=it->second.end(); iter++){
              elemstosend[i] = submeshesowned[k].erenumber[*iter]+1;
              i++;
            }
            std::stringstream ssname;
            ssname << it->first;
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit();
            if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToSend", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
            if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
            /* std::cout << it->second.size() << std::endl; */
            if(cg_ptset_write(CGNS_ENUMV(PointList), it->second.size(), elemstosend)) cg_error_exit();
            if(cg_ordinal_write(it->first));
            delete [] elemstosend;
          }
        }
        name="ElemsToRecv";
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_user_data_write(name.c_str())) cg_error_exit();
        for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[k].elemstorecv.begin(); it!=submeshesowned[k].elemstorecv.end(); it++){
          if(it->second.size()>0){
            std::stringstream ssname;
            ssname << it->first;
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToRecv", 0, "end")) cg_error_exit();
            if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToRecv", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
            if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
            elemstorecv = new cgsize_t[2];
            std::set<idx_t>::iterator iter=it->second.begin();
            idx_t itglobloc = submeshesowned[k].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*iter]];
            elemstorecv[0] = submeshesowned[k].erenumber[itglobloc]+1;
            elemstorecv[1] = it->second.size();
            if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
            if(cg_ordinal_write(it->first));
            delete [] elemstorecv;
          }
        }

        //Write boundary conditions
        for(int bc=0; bc<submeshesowned[k].boundary_conditions.size(); bc++){
          if(submeshesowned[k].boundary_conditions[bc].size()>0){
            cgsize_t bcnodes[submeshesowned[k].boundary_conditions[bc].size()];
            cgsize_t nbc=0;
            for(std::set<idx_t>::iterator it=submeshesowned[k].boundary_conditions[bc].begin();
              it!=submeshesowned[k].boundary_conditions[bc].end();it++){
              bcnodes[nbc] = submeshesowned[k].nodes_gtl[*it] ;
              nbc++;
            }
            if(cg_boco_write(index_file,index_base,index_zone,submeshesowned[k].boundary_conditions_names[bc].c_str(),submeshesowned[k].boundary_conditions_types[bc],CGNS_ENUMV(PointList),submeshesowned[k].boundary_conditions[bc].size(),bcnodes,&index_bc)) cg_error_exit();
            cg_boco_gridlocation_write(index_file,index_base,index_zone,index_bc,CGNS_ENUMV(Vertex));
          }
        }

        int nuserdata = submeshesowned[k].ud_names.size();
        for(int nud=1; nud<=nuserdata; nud++){
          int narrays = submeshesowned[k].ar_names[nud-1].size();
          if(narrays>0){
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
            /* std::cout << submeshesowned[k].submeshid << " " << submeshesowned[k].ud_names[nud-1] << std::endl; */
            if(cg_user_data_write(submeshesowned[k].ud_names[nud-1].c_str())) cg_error_exit();
            cgsize_t dimensions=submeshesowned[k].arrays[nud-1][0].size();
            CGNS_ENUMV(GridLocation_t) location;
            if(dimensions==submeshesowned[k].nodes.size()){
              location = CGNS_ENUMV(Vertex);
            }
            if(dimensions==submeshesowned[k].elems.size()){
              location = CGNS_ENUMV(CellCenter);
            }
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, submeshesowned[k].ud_names[nud-1].c_str(), 0, "end")) cg_error_exit();
            if(cg_gridlocation_write(location)) cg_error_exit();
            for(int na=1; na<=narrays; na++){
              int rank=1;
              if(cg_array_write(submeshesowned[k].ar_names[nud-1][na-1].c_str(), CGNS_ENUMV(RealDouble), 1, &dimensions, submeshesowned[k].arrays[nud-1][na-1].data())) cg_error_exit();

            }
          }
        }

      }
      if(cg_close(index_file)) cg_error_exit();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void readArrays(std::vector<submesh> &submeshesowned, std::string filename, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  int index_file, index_base, n_bases, physDim, cellDim, nZones, index_zone;
  int gnelems, nelems;
  int nSections;
  cgsize_t gnnodes;
  int nCoords;
  char basename[40];
  char zonename[40];
  char name[40];
  char secname[40];
  cgsize_t sizes[2];
  CGNS_ENUMV(ZoneType_t) zoneType;
  CGNS_ENUMV(DataType_t) dataType;
  CGNS_ENUMV(ElementType_t) type;


  if (cg_open(filename.c_str(),CG_MODE_READ,&index_file)) cg_error_exit();
  if(cg_nbases(index_file, &n_bases)!= CG_OK) cg_get_error();
  if(n_bases != 1) cg_get_error(); 
  index_base=1;
  if(cg_base_read(index_file, index_base, basename, &cellDim, &physDim) != CG_OK) cg_get_error();
  if(cg_nzones(index_file, index_base, &nZones) != CG_OK) cg_get_error();
  index_zone = 1;
  if(cg_zone_type(index_file, index_base, index_zone, &zoneType) != CG_OK) cg_get_error();
  assert(zoneType == CGNS_ENUMV(Unstructured));
  if(cg_zone_read(index_file, index_base, index_zone, zonename, sizes) != CG_OK) cg_get_error();
  gnnodes = sizes[0];
  gnelems = sizes[1];

  int narrays, nuserdata;
  if(cg_goto(index_file, index_base, zonename, 0, "end")) cg_error_exit();
  if(cg_nuser_data(&nuserdata)) cg_error_exit();
  /* std::cout << nuserdata << std::endl; */
  for(int nud=1; nud<=nuserdata; nud++){
    char udname[40];
    if(cg_user_data_read(nud, udname)) cg_error_exit();
    if(cg_goto(index_file, index_base, zonename, 0, udname, 0, "end")) cg_error_exit();
    if(cg_narrays(&narrays)) cg_error_exit();
    /* std::cout << narrays << std::endl; */
    for(int na=1; na<=narrays; na++){
      char aname[40];
      CGNS_ENUMV(DataType_t) dataType;
      int rank;
      cgsize_t dimensions;
      /* assert(dataType=CGNS_ENUMV(RealDouble)); */
      CGNS_ENUMV(GridLocation_t) location;
      if(cg_gridlocation_read(&location)) cg_error_exit();
      assert(location==CGNS_ENUMV(CellCenter) or location==CGNS_ENUMV(Vertex));
      if(cg_array_info(na,aname,&dataType,&rank,&dimensions)) cg_error_exit();
      /* std::cout << dataType << std::endl; */
      if(dataType == CGNS_ENUMV(RealDouble)){
        double *array = new double[dimensions];
        if(cg_array_read(na, array)) cg_error_exit();
        /* std::cout << array[0] << std::endl; */
        std::vector<idx_t> nloc(submeshesowned.size(), 0);
        //Resize if first time
        if(nud==1 and na==1){
          for(idx_t k=0; k<submeshesowned.size(); k++){
            submeshesowned[k].ud_names.resize(nuserdata);
            submeshesowned[k].ar_names.resize(nuserdata, std::vector<std::string>(narrays));
            submeshesowned[k].arrays.resize(nuserdata, std::vector<std::vector<double> >(narrays));
          }
        }
        for(idx_t i=0; i<dimensions; i++){
          for(idx_t k=0; k<submeshesowned.size(); k++){
            if(location=CGNS_ENUMV(CellCenter)){
              if(submeshesowned[k].elems_gtl.count(i)!=0){
                if(submeshesowned[k].arrays[nud-1][na-1].size() == 0) {
                  submeshesowned[k].ud_names[nud-1] = std::string(udname);
                  submeshesowned[k].ar_names[nud-1][na-1] = std::string(aname);
                  /* std::cout << udname << " " << aname << std::endl; */
                  /* std::cout << submeshesowned[k].submeshid << " " << submeshesowned[k].ud_names[nud-1] << std::endl; */
                  submeshesowned[k].arrays[nud-1][na-1].resize(submeshesowned[k].elems.size());
                }
                int iloc = submeshesowned[k].elems_gtl[i];
                submeshesowned[k].arrays[nud-1][na-1][iloc] = array[i];
              }  
            }
          }
        }
        delete [] array;
      }
    }
  }

  if (cg_close(index_file)) cg_error_exit();
}


void writeVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t indrenum;
    idx_t nNodes = submeshesowned[k].nnodes;
    idx_t nElems = submeshesowned[k].nelems;
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][0] << " " << submeshesowned[k].nodes[i][1] << " " << submeshesowned[k].nodes[i][2] << std::endl;
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems << " " << nElems*(esize+1) << std::endl;
    for(idx_t i=0; i<nElems; i++){
      outFile << esize << " ";
      indrenum = submeshesowned[k].erenumberr[i];
      /* indrenum = i; */
      for(idx_t p=0; p<esize; p++){
        outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].elems[indrenum][p]] << " ";
      }
      outFile << std::endl;
    }

    outFile << std::endl;
    outFile << "CELL_TYPES " << nElems << std::endl;
    idx_t cellType;
    if(dim==3) cellType=10; //tetrahedrons
    if(dim==2 and esize==3) cellType=5; //triangles
    if(dim==2 and esize==4) cellType=9; //quadrangles
    for(idx_t i=0; i<nElems; i++){
      outFile << cellType << std::endl;
    }
  }
}

void writesendrecvCute(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){

  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].nnodes;
    idx_t nboundaryelems = submeshesowned[k].nboundaryelems;

    std::map<idx_t, idx_t> trueneighborssend;
    std::map<idx_t, idx_t> trueneighborsrecv;
    for(std::set<idx_t>::iterator it=submeshesowned[k].potentialneighbors.begin();
        it!=submeshesowned[k].potentialneighbors.end();it++){
      if(submeshesowned[k].elemstorecv[*it].size() > 0){
        trueneighborssend.insert(std::make_pair(*it,submeshesowned[k].elemstorecv[*it].size()));
        trueneighborsrecv.insert(std::make_pair(*it,submeshesowned[k].elemstosend[*it].size()));
      }
    }

    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_info_send_recv.txt";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << trueneighborssend.size() << std::endl;
    for(std::map<idx_t,idx_t>::iterator it=trueneighborssend.begin();
        it!=trueneighborssend.end();it++){
      outFile << it->first << " " << it->second << std::endl;
    }

    for(std::map<idx_t,idx_t>::iterator it=trueneighborssend.begin();
        it!=trueneighborssend.end();it++){
      for(std::set<idx_t>::iterator it1=submeshesowned[k].elemstosend[it->first].begin();
          it1!=submeshesowned[k].elemstosend[it->first].end();it1++){
        outFile << submeshesowned[k].erenumber[*it1] << std::endl;
      }
    }

    for(std::map<idx_t,idx_t>::iterator it=trueneighborsrecv.begin();
        it!=trueneighborsrecv.end();it++){
      for(std::set<idx_t>::iterator it1=submeshesowned[k].elemstorecv[it->first].begin();
          it1!=submeshesowned[k].elemstorecv[it->first].end();it1++){
        idx_t itglobloc = submeshesowned[k].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*it1]];
        outFile << it->first << " " << it->second << " " << submeshesowned[k].erenumber[itglobloc] << std::endl;
        break;
      }
    }


  }

}
void writerecvVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){

  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].nnodes;
    idx_t nboundaryelems = submeshesowned[k].nboundaryelems;
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_recv_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][0] << " " << submeshesowned[k].nodes[i][1] << " " << submeshesowned[k].nodes[i][2] << std::endl;
    }



    idx_t nbelemsrecv(0);
    std::set<idx_t>::iterator iter;
    for(iter=submeshesowned[k].potentialneighbors.begin();iter!=submeshesowned[k].potentialneighbors.end();iter++){
      nbelemsrecv +=  submeshesowned[k].elemstorecv[*iter].size();
    }

    outFile << std::endl;
    outFile << "CELLS " << nbelemsrecv << " " << nbelemsrecv*(esize+1) << std::endl;
    std::set<idx_t>::iterator it;
    for(iter=submeshesowned[k].potentialneighbors.begin();iter!=submeshesowned[k].potentialneighbors.end();iter++){
      for(it=submeshesowned[k].elemstorecv[*iter].begin(); it!=submeshesowned[k].elemstorecv[*iter].end(); it++){
        idx_t itloc = submeshesowned[k].elems_gtl[g_potentialneighbors[*iter].elems_ltg[*it]];
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].elems[itloc][p]] << " ";
        }
        outFile << std::endl;
      }
    }

    outFile << std::endl;
    outFile << "CELL_TYPES " << nbelemsrecv << std::endl;
    idx_t cellType;
    if(dim==3) cellType=10; //tetrahedrons
    if(dim==2 and esize==3) cellType=5; //triangles
    if(dim==2 and esize==4) cellType=9; //quadrangles
    for(idx_t i=0; i<nbelemsrecv; i++){
      outFile << cellType << std::endl;
    }

    outFile << std::endl;
    outFile << "POINT_DATA " << nNodes << std::endl;
    outFile << "SCALARS b float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][2] << std::endl;
    }
  }
}

void writesendVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].nnodes;
    idx_t nboundaryelems = submeshesowned[k].nboundaryelems;
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_send_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][0] << " " << submeshesowned[k].nodes[i][1] << " " << submeshesowned[k].nodes[i][2] << std::endl;
    }

    idx_t nbelemssend(0);
    std::set<idx_t>::iterator iter;
    for(iter=submeshesowned[k].potentialneighbors.begin();iter!=submeshesowned[k].potentialneighbors.end();iter++){
      nbelemssend +=  submeshesowned[k].elemstosend[*iter].size();
    }

    outFile << std::endl;
    outFile << "CELLS " << nbelemssend << " " << nbelemssend*(esize+1) << std::endl;
    std::set<idx_t>::iterator it;
    for(iter=submeshesowned[k].potentialneighbors.begin();iter!=submeshesowned[k].potentialneighbors.end();iter++){
      for(it=submeshesowned[k].elemstosend[*iter].begin(); it!=submeshesowned[k].elemstosend[*iter].end(); it++){
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].elems[*it][p]] << " ";
        }
        outFile << std::endl;
      }
    }

    outFile << std::endl;
    outFile << "CELL_TYPES " << nbelemssend << std::endl;
    idx_t cellType;
    if(dim==3) cellType=10; //tetrahedrons
    if(dim==2 and esize==3) cellType=5; //triangles
    if(dim==2 and esize==4) cellType=9; //quadrangles
    for(idx_t i=0; i<nbelemssend; i++){
      outFile << cellType << std::endl;
    }

    outFile << std::endl;
    outFile << "POINT_DATA " << nNodes << std::endl;
    outFile << "SCALARS b float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][2] << std::endl;
    }
  }
}
void writeboundaryVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].nnodes;
    idx_t nboundaryelems = submeshesowned[k].nboundaryelems;
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_boundary_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][0] << " " << submeshesowned[k].nodes[i][1] << " " << submeshesowned[k].nodes[i][2] << std::endl;
    }

    idx_t nbelemssend(0);
    std::set<idx_t>::iterator iter;
    for(iter=submeshesowned[k].potentialneighbors.begin();iter!=submeshesowned[k].potentialneighbors.end();iter++){
      nbelemssend +=  submeshesowned[k].elemstosend[*iter].size();
    }

    outFile << std::endl;
    outFile << "CELLS " << nbelemssend << " " << nbelemssend*(esize+1) << std::endl;
    std::set<idx_t>::iterator it;
    for(iter=submeshesowned[k].potentialneighbors.begin();iter!=submeshesowned[k].potentialneighbors.end();iter++){
      for(it=submeshesowned[k].elemstosend[*iter].begin(); it!=submeshesowned[k].elemstosend[*iter].end(); it++){
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].elems[*it][p]] << " ";
        }
        outFile << std::endl;
      }
    }

    outFile << std::endl;
    outFile << "CELL_TYPES " << nbelemssend << std::endl;
    idx_t cellType;
    if(dim==3) cellType=10; //tetrahedrons
    if(dim==2 and esize==3) cellType=5; //triangles
    if(dim==2 and esize==4) cellType=9; //quadrangles
    for(idx_t i=0; i<nbelemssend; i++){
      outFile << cellType << std::endl;
    }

    outFile << std::endl;
    outFile << "POINT_DATA " << nNodes << std::endl;
    outFile << "SCALARS b float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][2] << std::endl;
    }
  }
}

void updateNodes(std::vector<submesh> &submeshesowned, std::string nodesFilename, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  for(idx_t k=0; k<submeshesowned.size(); k++){
    submeshesowned[k].nnodes = submeshesowned[k].nodesindex.size();
    submeshesowned[k].nodes.resize(submeshesowned[k].nnodes, std::vector<real_t>(3, 0.));
  }

  std::ifstream nodesFile(nodesFilename);
  idx_t nNodesTot;
  nodesFile >> nNodesTot;

  std::vector<idx_t> nloc(submeshesowned.size(), 0);
  real_t sx, sy, sz;
  for(idx_t i=0; i<nNodesTot; i++){
    nodesFile >> sx >> sy >> sz;
    for(idx_t k=0; k<submeshesowned.size(); k++){
      if(submeshesowned[k].nodesindex.find(i)!=submeshesowned[k].nodesindex.end()){
        submeshesowned[k].nodes[nloc[k]][0] = sx;
        submeshesowned[k].nodes[nloc[k]][1] = sy;
        submeshesowned[k].nodes[nloc[k]][2] = sz;
        submeshesowned[k].nodes_ltg.insert({nloc[k], i});
        submeshesowned[k].nodes_gtl.insert({i, nloc[k]});
        nloc[k] += 1;
      }
    }
  }
  nodesFile.close();
}

void updateNodesCGNS(std::vector<submesh> &submeshesowned, std::string filename, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  int index_file, index_base, n_bases, base, physDim, cellDim, nZones, zone;
  int gnelems, nelems;
  int nSections;
  cgsize_t gnnodes;
  int nCoords;
  char basename[40];
  char zonename[40];
  char name[40];
  char secname[40];
  cgsize_t sizes[2];
  CGNS_ENUMV(ZoneType_t) zoneType;
  CGNS_ENUMV(DataType_t) dataType;
  CGNS_ENUMV(ElementType_t) type;

  if (cg_open(filename.c_str(),CG_MODE_READ,&index_file)) cg_error_exit();
  if(cg_nbases(index_file, &n_bases)!= CG_OK) cg_get_error();
  if(n_bases != 1) cg_get_error(); 
  base=1;
  if(cg_base_read(index_file, base, basename, &cellDim, &physDim) != CG_OK) cg_get_error();
  if(cg_nzones(index_file, base, &nZones) != CG_OK) cg_get_error();
  zone = 1;
  if(cg_zone_type(index_file, base, zone, &zoneType) != CG_OK) cg_get_error();
  assert(zoneType == CGNS_ENUMV(Unstructured));
  if(cg_zone_read(index_file, base, zone, zonename, sizes) != CG_OK) cg_get_error();
  gnnodes = sizes[0];
  gnelems = sizes[1];

  if(cg_ncoords(index_file, base, zone, &nCoords) != CG_OK) cg_get_error();
  cgsize_t one = 1;
  double *x, *y, *z;

  if(cg_coord_info(index_file, base, zone, 1, &dataType, zonename) != CG_OK) cg_get_error();
  x = new double[gnnodes];
  if(cg_coord_read(index_file, base, zone, zonename, CGNS_ENUMV(RealDouble),
        &one, &gnnodes, x) != CG_OK) cg_get_error();

  if(cg_coord_info(index_file, base, zone, 2, &dataType, zonename) != CG_OK) cg_get_error();
  y = new double[gnnodes];
  if(cg_coord_read(index_file, base, zone, zonename, CGNS_ENUMV(RealDouble),
        &one, &gnnodes, y) != CG_OK) cg_get_error();

  if(cg_coord_info(index_file, base, zone, 3, &dataType, zonename) != CG_OK) cg_get_error();
  z = new double[gnnodes];
  if(cg_coord_read(index_file, base, zone, zonename, CGNS_ENUMV(RealDouble),
        &one, &gnnodes, z) != CG_OK) cg_get_error();

  for(idx_t k=0; k<submeshesowned.size(); k++){
    submeshesowned[k].nnodes = submeshesowned[k].nodesindex.size();
    submeshesowned[k].nodes.resize(submeshesowned[k].nnodes, std::vector<real_t>(3, 0.));
  }

  std::map<idx_t, std::set<idx_t>> global_nodes;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(std::unordered_set<idx_t>::iterator it=submeshesowned[k].nodesindex.begin();
        it!=submeshesowned[k].nodesindex.end(); it++){
      global_nodes[*it].insert(k);
    }
  }


  std::vector<idx_t> nloc(submeshesowned.size(), 0);
  real_t sx, sy, sz;
  for(idx_t i=0; i<gnnodes; i++){
    if(global_nodes.count(i)!=0){
      for(std::set<idx_t>::iterator it=global_nodes[i].begin();
          it!=global_nodes[i].end(); it++){
        /* for(idx_t k=0; k<submeshesowned.size(); k++){ */
        /*   if(submeshesowned[k].nodesindex.find(i)!=submeshesowned[k].nodesindex.end()){ */
        int k=*it;
        submeshesowned[k].nodes[nloc[k]][0] = x[i];
        submeshesowned[k].nodes[nloc[k]][1] = y[i];
        submeshesowned[k].nodes[nloc[k]][2] = z[i];
        submeshesowned[k].nodes_ltg.insert({nloc[k], i});
        submeshesowned[k].nodes_gtl.insert({i, nloc[k]});
        nloc[k] += 1;
        /* } */
        /* } */
      }
    }
  }
  if (cg_close(index_file)) cg_error_exit();
}

void buildBoundaryEdges(std::vector<std::vector<idx_t> > &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::set<idx_t> &elems_set){

  idx_t nelems=elems.size();
  idx_t esize=elems[0].size();

  std::set<idx_t>::iterator it=elems_set.begin();

  std::pair<idx_t, idx_t> dummyPair;
  idx_t sp;
  for(it=elems_set.begin();it!=elems_set.end();it++){
    for(idx_t s=0; s<esize; s++){
      sp = (s+1)%esize;
      dummyPair = std::make_pair(std::min(elems[*it][s], elems[*it][sp]), std::max(elems[*it][s],elems[*it][sp]));
      if(edges.find(dummyPair) != edges.end()){
        edges[dummyPair] = std::make_pair(edges[dummyPair].first, *it);
      }
      else{
        edges.insert({dummyPair, std::make_pair(*it, -1)});
      }
    }
  }
}

void buildEdges(std::vector<std::vector<idx_t> > &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges){

  idx_t nelems=elems.size();
  idx_t esize=elems[0].size();

  std::pair<idx_t, idx_t> dummyPair;
  idx_t sp;
  for(idx_t i=0;i<nelems;i++){
    for(idx_t s=0; s<esize; s++){
      sp = (s+1)%esize;
      dummyPair = std::make_pair(std::min(elems[i][s], elems[i][sp]), std::max(elems[i][s],elems[i][sp]));
      if(edges.find(dummyPair) != edges.end()){
        edges[dummyPair] = std::make_pair(edges[dummyPair].first, i);
      }
      else{
        edges.insert({dummyPair, std::make_pair(i, -1)});
      }
    }
  }
}

void buildBoundaryFaces(std::vector<std::vector<idx_t> > &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::set<idx_t> &elems_set){
  std::set<idx_t>::iterator it;
  idx_t nelems=elems.size();
  idx_t esize=elems[0].size();
  idx_t fsize=3; //tetrahedron
  std::vector<idx_t> dummy(fsize,-1);
  std::vector<idx_t> sp(fsize,-1);
  for(it=elems_set.begin();it!=elems_set.end();it++){
    for(idx_t s=0; s<esize; s++){
      sp[0] = s;
      for(idx_t k=1;k<fsize;k++){
        sp[k] = (s+k)%esize;
      }
      for(idx_t k=0;k<fsize;k++){
        dummy[k] = elems[*it][sp[k]];
      }
      std::sort(dummy.begin(), dummy.end());
      if(faces.find(dummy) != faces.end()){
        faces[dummy] = std::make_pair(faces[dummy].first, *it);
      }
      else{
        faces.insert({dummy, std::make_pair(*it, -1)});
      }
    }
  }
}

void buildFaces(std::vector<std::vector<idx_t> > &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces){
  idx_t nelems=elems.size();
  idx_t esize=elems[0].size();
  idx_t fsize=3; //tetrahedron
  std::vector<idx_t> dummy(fsize,-1);
  std::vector<idx_t> sp(fsize,-1);
  for(idx_t i=0;i<nelems;i++){
    for(idx_t s=0; s<esize; s++){
      sp[0] = s;
      for(idx_t k=1;k<fsize;k++){
        sp[k] = (s+k)%esize;
      }
      for(idx_t k=0;k<fsize;k++){
        dummy[k] = elems[i][sp[k]];
      }
      std::sort(dummy.begin(), dummy.end());
      if(faces.find(dummy) != faces.end()){
        faces[dummy] = std::make_pair(faces[dummy].first, i);
      }
      else{
        faces.insert({dummy, std::make_pair(i, -1)});
      }
    }
  }
}

void buildElemConnectivity2D(std::vector<std::vector<idx_t> > &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::vector<std::vector<idx_t> > &neighbors){
  std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash>:: iterator it = edges.begin();

  idx_t esize=elems[0].size();

  idx_t s1,s2,elem1,elem2,sp;
  while(it != edges.end())
  {
    s1 = it->first.first;
    s2 = it->first.second;
    elem1 = it->second.first;
    elem2 = it->second.second;

    if(elem2!=-1){
      for(idx_t s=0; s<esize; s++){
        sp = (s+1)%esize;
        if(((elems[elem1][s] == s1) and (elems[elem1][sp] == s2)) or 
            ((elems[elem1][s] == s2) and (elems[elem1][sp] == s1))){
          neighbors[elem1][s] = elem2;
        }
        if(((elems[elem2][s] == s1) and (elems[elem2][sp] == s2)) or 
            ((elems[elem2][s] == s2) and (elems[elem2][sp] == s1))){
          neighbors[elem2][s] = elem1;
        }
      }
    }
    it++;
  }
}

void buildElemConnectivity3D(std::vector<std::vector<idx_t> > &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::vector<std::vector<idx_t> > &neighbors){
  idx_t esize=elems[0].size();
  std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash>:: iterator it = faces.begin();

  idx_t faceSize=3;//Tetrahedron
  std::vector<idx_t> faceVert(faceSize,-1);
  idx_t elem1,elem2;
  bool cond1;
  while(it != faces.end())
  {
    elem1 = it->second.first;
    elem2 = it->second.second;

    if(elem2!=-1){
      for(idx_t s=0; s<esize; s++){ //Pour chaque face

        faceVert[0] = elems[elem1][s];
        for(idx_t k=1; k<faceSize; k++){
          faceVert[k] = elems[elem1][(s+k)%esize];
        }

        cond1 = true;
        for(idx_t k=0;k<faceSize;k++){
          cond1 = ((std::find(faceVert.begin(), faceVert.end(), it->first[k])!=faceVert.end()) and cond1);
        }
        if(cond1) neighbors[elem1][s] = elem2;

        faceVert[0] = elems[elem2][s];
        for(idx_t k=1; k<faceSize; k++){
          faceVert[k] = elems[elem2][(s+k)%esize];
        }

        cond1 = true;
        for(idx_t k=0;k<faceSize;k++){
          cond1 = ((std::find(faceVert.begin(), faceVert.end(), it->first[k])!=faceVert.end()) and cond1);
        }
        if(cond1) neighbors[elem2][s] = elem1;
      }
    }
    it++;
  }
}

void Buildconnectivity(std::vector<submesh> &submeshesowned, idx_t dimension){
  if(dimension==2){
    for(idx_t k=0;k<submeshesowned.size();k++){
      std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> edges;
      submeshesowned[k].neighbors.resize(submeshesowned[k].nelems, std::vector<idx_t>(submeshesowned[k].esize,-1));
      buildEdges(submeshesowned[k].elems, edges);
      buildElemConnectivity2D(submeshesowned[k].elems, edges, submeshesowned[k].neighbors);
    }
  }
  if(dimension==3){
    for(idx_t k=0;k<submeshesowned.size();k++){
      std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> faces;
      submeshesowned[k].neighbors.resize(submeshesowned[k].nelems, std::vector<idx_t>(submeshesowned[k].esize,-1));
      buildFaces(submeshesowned[k].elems, faces);
      buildElemConnectivity3D(submeshesowned[k].elems, faces, submeshesowned[k].neighbors);
    }
  }
}




void Findboundaryfromconnectivity(std::vector<submesh> &submeshesowned, idx_t method, idx_t numlayers){

  for(idx_t k=0;k<submeshesowned.size();k++){
    for(idx_t i=0;i<submeshesowned[k].nelems;i++){
      if(submeshesowned[k].isBoundary(i)){
        submeshesowned[k].boundaryelems.insert(i);

        /* if(method==1){ */
        //Add boundarynodes
        for(idx_t j=0; j<submeshesowned[k].esize; j++){
          if(submeshesowned[k].neighbors[i][j] == -1){
            submeshesowned[k].boundarynodes.insert(submeshesowned[k].elems[i][j]);
            submeshesowned[k].boundarynodes.insert(submeshesowned[k].elems[i][(j+1)%submeshesowned[k].esize]);
          }
        }
        /* } */
      }
    }  

    //If method=1 add cells that have only a point in the boundary points
    if(method==1){
      for(idx_t i=0;i<submeshesowned[k].nelems;i++){
        for(idx_t j=0;j<submeshesowned[k].esize;j++){
          if(submeshesowned[k].boundarynodes.find(submeshesowned[k].elems[i][j])!=submeshesowned[k].boundarynodes.end()){
            submeshesowned[k].boundaryelems.insert(i);
          }
        }
      }

      for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
          it!=submeshesowned[k].boundaryelems.end();it++){
        for(idx_t p=0; p<submeshesowned[k].esize; p++){
          submeshesowned[k].boundarynodes.insert(submeshesowned[k].elems[*it][p]);
        }
      }
    }

    for(idx_t l=1; l<numlayers; l++){
      std::set<idx_t> tmpnewboundaryelems = submeshesowned[k].boundaryelems;
      std::set<idx_t>::iterator it;
      for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
        for(idx_t j=0;j<submeshesowned[k].esize;j++){
          if(submeshesowned[k].neighbors[*it][j] != -1){
            submeshesowned[k].boundaryelems.insert(submeshesowned[k].neighbors[*it][j]);

            /* if(method==1){ */
            //Add boundarynodes
            for(idx_t j=0; j<submeshesowned[k].esize; j++){
              submeshesowned[k].boundarynodes.insert(submeshesowned[k].elems[*it][j]);
            }
            /* } */

          }
        } 
      }

      if(method==1){
        for(idx_t i=0;i<submeshesowned[k].nelems;i++){
          for(idx_t j=0;j<submeshesowned[k].esize;j++){
            if(submeshesowned[k].boundarynodes.find(submeshesowned[k].elems[i][j])!=submeshesowned[k].boundarynodes.end()){
              submeshesowned[k].boundaryelems.insert(i);
            }
          }
        }
        for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
            it!=submeshesowned[k].boundaryelems.end();it++){
          for(idx_t p=0; p<submeshesowned[k].esize; p++){
            submeshesowned[k].boundarynodes.insert(submeshesowned[k].elems[*it][p]);
          }
        }
      }
    }

    submeshesowned[k].nboundaryelems = submeshesowned[k].boundaryelems.size();
    //Suite à enlver peut etre
    for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
        it!=submeshesowned[k].boundaryelems.end();it++){
      for(idx_t p=0; p<submeshesowned[k].esize; p++){
        submeshesowned[k].boundarynodes.insert(submeshesowned[k].elems[*it][p]);
      }
    }

  }
}

void submesh::computemyextents(idx_t nsubmeshes){
  extents.resize(3,std::vector<real_t>(2, 0.));

  extents[0][0] = nodes[0][0];
  extents[0][1] = nodes[0][0];

  extents[1][0] = nodes[0][1];
  extents[1][1] = nodes[0][1];

  extents[2][0] = nodes[0][2];
  extents[2][1] = nodes[0][2];

  for(int i=1;i<nnodes;i++){
    if(nodes[i][0] < extents[0][0]) extents[0][0] = nodes[i][0];
    if(nodes[i][0] > extents[0][1]) extents[0][1] = nodes[i][0];

    if(nodes[i][1] < extents[1][0]) extents[1][0] = nodes[i][1];
    if(nodes[i][1] > extents[1][1]) extents[1][1] = nodes[i][1];

    if(nodes[i][2] < extents[2][0]) extents[2][0] = nodes[i][2];
    if(nodes[i][2] > extents[2][1]) extents[2][1] = nodes[i][2];
  }

}

void FindNodesElemsSendRecv(std::vector<submesh> &submeshesowned, idx_t dimension, idx_t method, idx_t numlayers){

  for(idx_t k; k<submeshesowned.size(); k++){
    //Build my boundary connectivity
    submeshesowned[k].boundary_neighbors.resize(submeshesowned[k].nelems, std::vector<idx_t>(submeshesowned[k].esize,-1));
    if(dimension==2){
      buildBoundaryEdges(submeshesowned[k].elems, submeshesowned[k].boundary_edges, submeshesowned[k].boundaryelems);
      buildElemConnectivity2D(submeshesowned[k].elems, submeshesowned[k].boundary_edges, submeshesowned[k].boundary_neighbors);

      std::set<idx_t>::iterator iter;
      for(iter=submeshesowned[k].potentialneighbors.begin(); 
          iter!= submeshesowned[k].potentialneighbors.end(); iter++){

        std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> pn_edges;
        std::vector<std::vector<idx_t> > pn_neighbors(g_potentialneighbors[*iter].elems.size(), std::vector<idx_t>(submeshesowned[k].esize,-1));

        buildEdges(g_potentialneighbors[*iter].elems, pn_edges);
        buildElemConnectivity2D(g_potentialneighbors[*iter].elems, pn_edges, pn_neighbors);

        std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash>::iterator it_edges;
        for(it_edges=pn_edges.begin();it_edges!=pn_edges.end();it_edges++){
          if(submeshesowned[k].boundary_edges.find(it_edges->first)!=submeshesowned[k].boundary_edges.end()){

            //Add elem to recv
            idx_t elem=it_edges->second.first;
            submeshesowned[k].elemstorecv[*iter].insert(elem);
            //Add boundarynodes to recv
            submeshesowned[k].nodestorecv[*iter].insert(it_edges->first.first);
            submeshesowned[k].nodestorecv[*iter].insert(it_edges->first.second);

            //Add elem to send
            elem=submeshesowned[k].boundary_edges[it_edges->first].first;
            submeshesowned[k].elemstosend[*iter].insert(elem);
            //Add boundarynodes to send
            submeshesowned[k].nodestosend[*iter].insert(it_edges->first.first);
            submeshesowned[k].nodestosend[*iter].insert(it_edges->first.second);


          }
        }

        //RECV If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(idx_t i=0;i<g_potentialneighbors[*iter].elems.size();i++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].elems[i][j])!=submeshesowned[k].nodestorecv[*iter].end()){
                submeshesowned[k].elemstorecv[*iter].insert(i);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
              it!=submeshesowned[k].elemstorecv[*iter].end();
              it++){
            for(idx_t p=0; p<submeshesowned[k].esize; p++){
              submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][p]);
            }
          }
        }

        //SEND If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
              it!=submeshesowned[k].boundaryelems.end();
              it++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].elems[*it][j])!=submeshesowned[k].nodestosend[*iter].end()){
                submeshesowned[k].elemstosend[*iter].insert(*it);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
              it!=submeshesowned[k].elemstosend[*iter].end();
              it++){
            for(idx_t p=0; p<submeshesowned[k].esize; p++){
              submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].elems[*it][p]);
            }
          }
        }

        //Add all the layers
        for(idx_t l=1; l<numlayers; l++){
          //Add to send list
          std::set<idx_t> tmpnewboundaryelems = submeshesowned[k].elemstosend[*iter];
          std::set<idx_t>::iterator it;
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(submeshesowned[k].neighbors[*it][j] != -1){
                submeshesowned[k].elemstosend[*iter].insert(submeshesowned[k].neighbors[*it][j]);
                //Add boundarynodes
                for(idx_t j=0; j<submeshesowned[k].esize; j++){
                  submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].elems[*it][j]);
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
                it!=submeshesowned[k].boundaryelems.end();
                it++){
              for(idx_t j=0;j<submeshesowned[k].esize;j++){
                if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].elems[*it][j])!=submeshesowned[k].nodestosend[*iter].end()){
                  submeshesowned[k].elemstosend[*iter].insert(*it);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
                it!=submeshesowned[k].elemstosend[*iter].end();
                it++){
              for(idx_t p=0; p<submeshesowned[k].esize; p++){
                submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].elems[*it][p]);
              }
            }
          }

          //Add to recv list
          tmpnewboundaryelems.clear();
          tmpnewboundaryelems = submeshesowned[k].elemstorecv[*iter];
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(pn_neighbors[*it][j] != -1){
                submeshesowned[k].elemstorecv[*iter].insert(pn_neighbors[*it][j]);
                //Add boundarynodes
                for(idx_t j=0; j<submeshesowned[k].esize; j++){
                  submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][j]);
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(idx_t i=0; i<g_potentialneighbors[*iter].elems.size(); i++){
              for(idx_t j=0;j<submeshesowned[k].esize;j++){
                if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].elems[i][j])!=submeshesowned[k].nodestorecv[*iter].end()){
                  submeshesowned[k].elemstorecv[*iter].insert(i);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
                it!=submeshesowned[k].elemstorecv[*iter].end();
                it++){
              for(idx_t p=0; p<submeshesowned[k].esize; p++){
                submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][p]);
              }
            }
          }


        }

        for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
            it!=submeshesowned[k].elemstorecv[*iter].end();it++){
          for(idx_t p=0; p<submeshesowned[k].esize; p++){
            submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][p]);
          }
        }


      }
    }
    if(dimension==3){
      idx_t faceSize=3;//Tetrahedron
      buildBoundaryFaces(submeshesowned[k].elems, submeshesowned[k].boundary_faces, submeshesowned[k].boundaryelems);
      buildElemConnectivity3D(submeshesowned[k].elems, submeshesowned[k].boundary_faces, submeshesowned[k].boundary_neighbors);

      std::set<idx_t>::iterator iter;
      for(iter=submeshesowned[k].potentialneighbors.begin(); 
          iter!= submeshesowned[k].potentialneighbors.end(); iter++){

        std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> pn_faces;
        std::vector<std::vector<idx_t> > pn_neighbors(g_potentialneighbors[*iter].elems.size(), std::vector<idx_t>(submeshesowned[k].esize,-1));

        buildFaces(g_potentialneighbors[*iter].elems, pn_faces);
        buildElemConnectivity3D(g_potentialneighbors[*iter].elems, pn_faces, pn_neighbors);

        std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash>::iterator it_faces;

        for(it_faces=pn_faces.begin();it_faces!=pn_faces.end();it_faces++){
          if(submeshesowned[k].boundary_faces.find(it_faces->first)!=submeshesowned[k].boundary_faces.end()){

            //Add elem to recv
            idx_t elem=it_faces->second.first;
            submeshesowned[k].elemstorecv[*iter].insert(elem);
            //Add boundarynodes to recv
            for(idx_t p=0; p<faceSize; p++){
              submeshesowned[k].nodestorecv[*iter].insert(it_faces->first[p]);
            }

            //Add elem to send
            elem=submeshesowned[k].boundary_faces[it_faces->first].first;
            submeshesowned[k].elemstosend[*iter].insert(elem);
            //Add boundarynodes to send
            for(idx_t p=0; p<faceSize; p++){
              submeshesowned[k].nodestosend[*iter].insert(it_faces->first[p]);
            }


          }
        }
        //RECV If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(idx_t i=0;i<g_potentialneighbors[*iter].elems.size();i++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].elems[i][j])!=submeshesowned[k].nodestorecv[*iter].end()){
                submeshesowned[k].elemstorecv[*iter].insert(i);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
              it!=submeshesowned[k].elemstorecv[*iter].end();
              it++){
            for(idx_t p=0; p<submeshesowned[k].esize; p++){
              submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][p]);
            }
          }
        }

        //SEND If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
              it!=submeshesowned[k].boundaryelems.end();
              it++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].elems[*it][j])!=submeshesowned[k].nodestosend[*iter].end()){
                submeshesowned[k].elemstosend[*iter].insert(*it);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
              it!=submeshesowned[k].elemstosend[*iter].end();
              it++){
            for(idx_t p=0; p<submeshesowned[k].esize; p++){
              submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].elems[*it][p]);
            }
          }
        }

        //Add all the layers
        for(idx_t l=1; l<numlayers; l++){
          //Add to send list
          std::set<idx_t> tmpnewboundaryelems = submeshesowned[k].elemstosend[*iter];
          std::set<idx_t>::iterator it;
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(submeshesowned[k].neighbors[*it][j] != -1){
                submeshesowned[k].elemstosend[*iter].insert(submeshesowned[k].neighbors[*it][j]);
                //Add boundarynodes
                for(idx_t j=0; j<submeshesowned[k].esize; j++){
                  submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].elems[*it][j]);
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
                it!=submeshesowned[k].boundaryelems.end();
                it++){
              for(idx_t j=0;j<submeshesowned[k].esize;j++){
                if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].elems[*it][j])!=submeshesowned[k].nodestosend[*iter].end()){
                  submeshesowned[k].elemstosend[*iter].insert(*it);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
                it!=submeshesowned[k].elemstosend[*iter].end();
                it++){
              for(idx_t p=0; p<submeshesowned[k].esize; p++){
                submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].elems[*it][p]);
              }
            }
          }

          //Add to recv list
          tmpnewboundaryelems.clear();
          tmpnewboundaryelems = submeshesowned[k].elemstorecv[*iter];
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<submeshesowned[k].esize;j++){
              if(pn_neighbors[*it][j] != -1){
                submeshesowned[k].elemstorecv[*iter].insert(pn_neighbors[*it][j]);
                //Add boundarynodes
                for(idx_t j=0; j<submeshesowned[k].esize; j++){
                  submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][j]);
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(idx_t i=0; i<g_potentialneighbors[*iter].elems.size(); i++){
              for(idx_t j=0;j<submeshesowned[k].esize;j++){
                if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].elems[i][j])!=submeshesowned[k].nodestorecv[*iter].end()){
                  submeshesowned[k].elemstorecv[*iter].insert(i);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
                it!=submeshesowned[k].elemstorecv[*iter].end();
                it++){
              for(idx_t p=0; p<submeshesowned[k].esize; p++){
                submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][p]);
              }
            }
          }


        }

        for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
            it!=submeshesowned[k].elemstorecv[*iter].end();it++){
          for(idx_t p=0; p<submeshesowned[k].esize; p++){
            submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].elems[*it][p]);
          }
        }
      }

    }
  }

}

void Shareboundary(std::vector<submesh> &submeshesowned, std::vector<idx_t> &ownerofsubmesh, MPI_Comm comm){

  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;
  MPI_Request *requestSendPtr, *requestRecvPtr;

  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  std::vector<std::set<idx_t> > msgs_send(nprocs);
  std::vector<std::set<idx_t> > msgs_recv(nprocs);

  idx_t esize=submeshesowned[0].esize;


  std::set<idx_t>::iterator it;
  idx_t iloc(0);
  idx_t idum;

  //Construct submeshlocid (reverse lookup of mesh id)
  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  //Find max number of nodes and elems to send
  idx_t maxnodestosend(0), maxelemstosend(0);
  for(idx_t k=0; k<submeshesowned.size(); k++){
    maxnodestosend = std::max(maxnodestosend, (idx_t) submeshesowned[k].boundarynodes.size());
    maxelemstosend = std::max(maxelemstosend, (idx_t) submeshesowned[k].boundaryelems.size());
  }
  real_t *nodestosend = new real_t[maxnodestosend*4*submeshesowned.size()];
  idx_t *elemstosend = new idx_t[maxelemstosend*(1+esize)*submeshesowned.size()];

  //Compute total number of messages/request
  idx_t ntotreq_send=0;
  idx_t ntotreq_recv=0;
  idx_t nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(ownerofsubmesh[*it]!=me && (msgs_send[ownerofsubmesh[*it]].count(submeshesowned[k].submeshid)==0)){
        msgs_send[ownerofsubmesh[*it]].insert(submeshesowned[k].submeshid);
        ntotreq_send++;
      }
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        ntotreq_recv++;
      }
    }
  }

  requestSendPtr = new MPI_Request[ntotreq_send];

  //Construc contiguous vector of nodes  to send
  for(idx_t k=0; k<submeshesowned.size(); k++){
    iloc=0;
    for(it=submeshesowned[k].boundarynodes.begin(); it!=submeshesowned[k].boundarynodes.end(); it++){
      idum = submeshesowned[k].nodes_gtl[*it];
      nodestosend[4*k*maxnodestosend+4*iloc] = (idx_t) *it;
      nodestosend[4*k*maxnodestosend+4*iloc+1] = submeshesowned[k].nodes[idum][0];
      nodestosend[4*k*maxnodestosend+4*iloc+2] = submeshesowned[k].nodes[idum][1];
      nodestosend[4*k*maxnodestosend+4*iloc+3] = submeshesowned[k].nodes[idum][2];
      iloc++;
    }
  }


  //Clear msgsend msgs recv
  for(idx_t pr=0;pr<nprocs;pr++){
    msgs_send[pr].clear();
    msgs_recv[pr].clear();
  }

  //Send Nodes vector
  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(msgs_send[ownerofsubmesh[*it]].count(submeshesowned[k].submeshid)==0){
        if(ownerofsubmesh[*it]!=me){
          MPI_Isend(&nodestosend[4*maxnodestosend*k], submeshesowned[k].boundarynodes.size()*4, REAL_T, ownerofsubmesh[*it], (*it)*ownerofsubmesh.size()+submeshesowned[k].submeshid, comm, &requestSendPtr[nreq]);
          nreq++;
        }
        else{
          potential_neighbors_boundary pnbound;
          pnbound.nodes.resize(submeshesowned[k].boundarynodes.size(),std::vector<real_t>(3,0.));
          for(idx_t i=0;i<submeshesowned[k].boundarynodes.size();i++){
            pnbound.nodes_ltg.insert(std::make_pair(i,nodestosend[4*maxnodestosend*k+i*4]));
            pnbound.nodes_gtl.insert(std::make_pair(nodestosend[4*maxnodestosend*k+i*4],i));
            pnbound.nodes[i][0] = nodestosend[4*maxnodestosend*k+i*4+1];
            pnbound.nodes[i][1] = nodestosend[4*maxnodestosend*k+i*4+2];
            pnbound.nodes[i][2] = nodestosend[4*maxnodestosend*k+i*4+3];
          }
          g_potentialneighbors.insert(std::make_pair(submeshesowned[k].submeshid,pnbound));
        }
        msgs_send[ownerofsubmesh[*it]].insert(submeshesowned[k].submeshid);
      }
    }
  }

  //Probe to get sizes
  idx_t *nnodestorecv = new idx_t[ntotreq_recv];
  idx_t maxnodestorecv=0;

  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        MPI_Probe(ownerofsubmesh[*it], (*it)+ownerofsubmesh.size()*submeshesowned[k].submeshid, comm, &stat);
        MPI_Get_count(&stat,REAL_T,&nnodestorecv[nreq]);
        nnodestorecv[nreq] /= 4;
        maxnodestorecv = std::max(maxnodestorecv,nnodestorecv[nreq]);
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        nreq++;
      }
    }
  }

  real_t *nodestorecv = new real_t[maxnodestorecv*4*ntotreq_recv];

  //Clear msgsend msgs recv
  for(idx_t pr=0;pr<nprocs;pr++){
    msgs_send[pr].clear();
    msgs_recv[pr].clear();
  }

  //IRecv nodes vector
  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        MPI_Recv(&nodestorecv[nreq*maxnodestorecv*4], nnodestorecv[nreq]*4, REAL_T, ownerofsubmesh[*it], (*it)+ownerofsubmesh.size()*submeshesowned[k].submeshid, comm, &stat);

        //Add nodes
        potential_neighbors_boundary pnbound;
        pnbound.nodes.resize(nnodestorecv[nreq],std::vector<real_t>(3,0.));
        for(idx_t i=0;i<nnodestorecv[nreq];i++){
          pnbound.nodes_ltg.insert(std::make_pair(i,(idx_t) nodestorecv[4*maxnodestorecv*nreq+i*4]));
          pnbound.nodes_gtl.insert(std::make_pair((idx_t) nodestorecv[4*maxnodestorecv*nreq+i*4],i));
          pnbound.nodes[i][0] = nodestorecv[4*maxnodestorecv*nreq+i*4+1];
          pnbound.nodes[i][1] = nodestorecv[4*maxnodestorecv*nreq+i*4+2];
          pnbound.nodes[i][2] = nodestorecv[4*maxnodestorecv*nreq+i*4+3];
        }
        g_potentialneighbors.insert(std::make_pair(*it,pnbound));
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        nreq++;
      }
    }
  }

  MPI_Waitall(ntotreq_send, requestSendPtr, statSend);

  //Construc contiguous vector of elems to send
  for(idx_t k=0; k<submeshesowned.size(); k++){
    iloc=0;
    for(it=submeshesowned[k].boundaryelems.begin(); it!=submeshesowned[k].boundaryelems.end(); it++){
      idum = submeshesowned[k].elems_ltg[*it];
      elemstosend[(1+esize)*(maxelemstosend*k+iloc)] = idum;
      for(idx_t j=0;j<esize;j++){
        elemstosend[(1+esize)*(maxelemstosend*k+iloc)+j+1] = submeshesowned[k].elems[*it][j];
      }
      iloc++;
    }
  }

  //Clear msgsend msgs recv
  for(idx_t pr=0;pr<nprocs;pr++){
    msgs_send[pr].clear();
    msgs_recv[pr].clear();
  }

  //Send elems vector
  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(msgs_send[ownerofsubmesh[*it]].count(submeshesowned[k].submeshid)==0){
        if(ownerofsubmesh[*it]!=me){
          MPI_Isend(&elemstosend[(1+esize)*maxelemstosend*k], submeshesowned[k].boundaryelems.size()*(1+esize), IDX_T, ownerofsubmesh[*it], (*it)*ownerofsubmesh.size()+submeshesowned[k].submeshid, comm, &requestSendPtr[nreq]);
          nreq++;
        }
        else{
          g_potentialneighbors[submeshesowned[k].submeshid].elems.resize(submeshesowned[k].boundaryelems.size(),std::vector<idx_t>(esize,0));
          for(idx_t i=0;i<submeshesowned[k].boundaryelems.size();i++){
            g_potentialneighbors[submeshesowned[k].submeshid].elems_ltg.insert(std::make_pair(i,elemstosend[(1+esize)*(maxelemstosend*k+i)]));
            g_potentialneighbors[submeshesowned[k].submeshid].elems_gtl.insert(std::make_pair(elemstosend[(1+esize)*(maxelemstosend*k+i)],i));
            for(idx_t j=0;j<esize;j++){
              g_potentialneighbors[submeshesowned[k].submeshid].elems[i][j] = elemstosend[(1+esize)*(maxelemstosend*k+i)+j+1];
            }
          }
        }
        msgs_send[ownerofsubmesh[*it]].insert(submeshesowned[k].submeshid);
      }
    }
  }

  //Probe to get sizes
  idx_t *nelemstorecv = new idx_t[ntotreq_recv];
  idx_t maxelemstorecv=0;

  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        MPI_Probe(ownerofsubmesh[*it], (*it)+ownerofsubmesh.size()*submeshesowned[k].submeshid, comm, &stat);
        MPI_Get_count(&stat,IDX_T,&nelemstorecv[nreq]);
        nelemstorecv[nreq] /= (1+esize);
        maxelemstorecv = std::max(maxelemstorecv,nelemstorecv[nreq]);
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        nreq++;
      }
    }
  }

  //Clear msgsend msgs recv
  for(idx_t pr=0;pr<nprocs;pr++){
    msgs_send[pr].clear();
    msgs_recv[pr].clear();
  }

  idx_t *elemstorecv = new idx_t[(1+esize)*maxelemstorecv*ntotreq_recv];

  //IRecv elems vector
  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors.begin(); it!= submeshesowned[k].potentialneighbors.end(); it++){
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        MPI_Recv(&elemstorecv[nreq*(1+esize)*maxelemstorecv], nelemstorecv[nreq]*(1+esize), IDX_T, ownerofsubmesh[*it], (*it)+ownerofsubmesh.size()*submeshesowned[k].submeshid, comm, &stat);

        //Add elems
        g_potentialneighbors[*it].elems.resize(nelemstorecv[nreq],std::vector<idx_t>(esize,0));
        for(idx_t i=0;i<nelemstorecv[nreq];i++){
          g_potentialneighbors[*it].elems_ltg.insert(std::make_pair(i,elemstorecv[(1+esize)*(maxelemstorecv*nreq+i)]));
          g_potentialneighbors[*it].elems_gtl.insert(std::make_pair(elemstorecv[(1+esize)*(maxelemstorecv*nreq+i)],i));
          for(idx_t j=0;j<esize;j++){
            g_potentialneighbors[*it].elems[i][j] = elemstorecv[(1+esize)*(maxelemstorecv*nreq+i)+j+1];
          }
        }
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        nreq++;
      }
    }
  }

  MPI_Waitall(ntotreq_send, requestSendPtr, statSend);

  delete [] nodestosend;
  delete [] nodestorecv;
  delete [] nnodestorecv;
  delete [] elemstosend;
  delete [] elemstorecv;
  delete [] nelemstorecv;
}

void Computepotentialneighbors(idx_t nsubmeshes, std::vector<submesh> &submeshesowned, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;

  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  /* std::vector<std::vector<std::vector<real_t> > > allextents(nsubmeshes,std::vector<std::vector<real_t> >(3, std::vector<real_t>(2, 0.))); */
  real_t *allextents = new real_t[nsubmeshes*3*2];
  for(idx_t k=0;k<nsubmeshes*2*3;k++){
    allextents[k] = 0.;
  }

  //Compute each submesh extents
  for(idx_t k=0; k<submeshesowned.size(); k++){
    submeshesowned[k].computemyextents(nsubmeshes);
    allextents[submeshesowned[k].submeshid*3*2] = submeshesowned[k].extents[0][0];
    allextents[submeshesowned[k].submeshid*3*2 + 1] = submeshesowned[k].extents[0][1];
    allextents[submeshesowned[k].submeshid*3*2 + 2] = submeshesowned[k].extents[1][0];
    allextents[submeshesowned[k].submeshid*3*2 + 3] = submeshesowned[k].extents[1][1];
    allextents[submeshesowned[k].submeshid*3*2 + 4] = submeshesowned[k].extents[2][0];
    allextents[submeshesowned[k].submeshid*3*2 + 5] = submeshesowned[k].extents[2][1];
  }

  //Share extents with everyone
  //MPI_Allreduce(MPI_IN_PLACE, &allextents[0][0][0], nsubmeshes*3*2, REAL_T, MPI_SUM, comm);
  //MPI_Allreduce(MPI_IN_PLACE, allextents.data(), 3*2, REAL_T, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, allextents, nsubmeshes*3*2, REAL_T, MPI_SUM, comm);

  real_t eps=1e-10;

  bool interx, intery, interz;
  //Compute potential neighbors
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(idx_t p=0; p<nsubmeshes; p++){
      if(submeshesowned[k].submeshid!=p){
        interx = ((allextents[submeshesowned[k].submeshid*6+1]>=allextents[p*6]-eps) and
            (allextents[submeshesowned[k].submeshid*6+1]<=allextents[p*6+1]+eps));
        interx = interx or ((allextents[submeshesowned[k].submeshid*6]>=allextents[p*6]-eps) and
            (allextents[submeshesowned[k].submeshid*6]<=allextents[p*6+1]+eps));
        interx = interx or ((allextents[submeshesowned[k].submeshid*6]<=allextents[p*6]+eps) and
            (allextents[submeshesowned[k].submeshid*6+1]>=allextents[p*6+1]-eps));

        intery = ((allextents[submeshesowned[k].submeshid*6+3]>=allextents[p*6+2]-eps) and
            (allextents[submeshesowned[k].submeshid*6+3]<=allextents[p*6+3]+eps));
        intery = intery or ((allextents[submeshesowned[k].submeshid*6+2]>=allextents[p*6+2]-eps) and (allextents[submeshesowned[k].submeshid*6+2]<=allextents[p*6+3]+eps));
        intery = intery or ((allextents[submeshesowned[k].submeshid*6+2]<=allextents[p*6+2]+eps) and (allextents[submeshesowned[k].submeshid*6+3]>=allextents[p*6+3]-eps));

        interz = ((allextents[submeshesowned[k].submeshid*6+5]>=allextents[p*6+4]-eps) and
            (allextents[submeshesowned[k].submeshid*6+5]<=allextents[p*6+5]+eps));
        interz = interz or ((allextents[submeshesowned[k].submeshid*6+4]>=allextents[p*6+4]-eps) and (allextents[submeshesowned[k].submeshid*6+4]<=allextents[p*6+5]+eps));
        interz = interz or ((allextents[submeshesowned[k].submeshid*6+4]<=allextents[p*6+4]+eps) and (allextents[submeshesowned[k].submeshid*6+5]>=allextents[p*6+5]-eps));

        if((interx and intery) and interz){
          submeshesowned[k].potentialneighbors.insert(p);
        }
      }
    }
  }
}

void AddElemsAndRenumber(std::vector<submesh> &submeshesowned){
  //Renumber send elems
  for(idx_t k=0; k<submeshesowned.size(); k++){
    //Add Recv elements to local partition
    for(std::set<idx_t>::iterator iter=submeshesowned[k].potentialneighbors.begin();
        iter!=submeshesowned[k].potentialneighbors.end();
        iter++){
      for(std::set<idx_t>::iterator it=submeshesowned[k].nodestorecv[*iter].begin();
          it!=submeshesowned[k].nodestorecv[*iter].end();
          it++){
        if(submeshesowned[k].nodes_gtl.find(*it)==submeshesowned[k].nodes_gtl.end()){
          submeshesowned[k].nodes_ltg.insert(std::make_pair(submeshesowned[k].nodes.size(), *it));
          submeshesowned[k].nodes_gtl.insert(std::make_pair(*it, submeshesowned[k].nodes.size()));
          submeshesowned[k].nodes.push_back(g_potentialneighbors[*iter].nodes[g_potentialneighbors[*iter].nodes_gtl[*it]]);
        }
      }
      for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
          it!=submeshesowned[k].elemstorecv[*iter].end();
          it++){
        idx_t itglob = g_potentialneighbors[*iter].elems_ltg[*it];
        submeshesowned[k].elems_ltg.insert(std::make_pair(submeshesowned[k].elems.size(), itglob));
        submeshesowned[k].elems_gtl.insert(std::make_pair(itglob, submeshesowned[k].elems.size()));
        submeshesowned[k].elems.push_back(g_potentialneighbors[*iter].elems[*it]);
      }
    }
    submeshesowned[k].nnodes = submeshesowned[k].nodes.size();
    submeshesowned[k].nelems = submeshesowned[k].elems.size();

    //Renumbering
    submeshesowned[k].erenumber.resize(submeshesowned[k].nelems, -1);
    submeshesowned[k].erenumberr.resize(submeshesowned[k].nelems, -1);
    idx_t itloc;
    idx_t ind=0;
    bool norsendorrecv;
    for(idx_t i=0; i<submeshesowned[k].nelems; i++){
      norsendorrecv = true;
      for(std::set<idx_t>::iterator iter=submeshesowned[k].potentialneighbors.begin();
          iter!=submeshesowned[k].potentialneighbors.end();
          iter++){

        if(submeshesowned[k].elemstosend[*iter].find(i)!=submeshesowned[k].elemstosend[*iter].end()){
          norsendorrecv = false;
        }

        if(g_potentialneighbors[*iter].elems_gtl.find(submeshesowned[k].elems_ltg[i])!=g_potentialneighbors[*iter].elems_gtl.end()){
          itloc = g_potentialneighbors[*iter].elems_gtl[submeshesowned[k].elems_ltg[i]];
          if(submeshesowned[k].elemstorecv[*iter].find(itloc)!=submeshesowned[k].elemstorecv[*iter].end()){
            norsendorrecv = false;
          }
        }
      }

      if(norsendorrecv){
        submeshesowned[k].erenumber[i] = ind;
        submeshesowned[k].erenumberr[ind] = i;
        ind++;
      }
    }

    for(std::set<idx_t>::iterator iter=submeshesowned[k].potentialneighbors.begin();
        iter!=submeshesowned[k].potentialneighbors.end();
        iter++){
      for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
          it!=submeshesowned[k].elemstosend[*iter].end();
          it++){
        if(submeshesowned[k].erenumber[*it] == -1){
          submeshesowned[k].erenumber[*it] = ind;
          submeshesowned[k].erenumberr[ind] = *it;
          ind++;
        }
      }
    }

    for(std::set<idx_t>::iterator iter=submeshesowned[k].potentialneighbors.begin();
        iter!=submeshesowned[k].potentialneighbors.end();
        iter++){
      for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
          it!=submeshesowned[k].elemstorecv[*iter].end();
          it++){
        idx_t itglobloc = submeshesowned[k].elems_gtl[g_potentialneighbors[*iter].elems_ltg[*it]];
        submeshesowned[k].erenumber[itglobloc] = ind;
        submeshesowned[k].erenumberr[ind] = itglobloc;
        ind++;
      }
    }


  }
}
