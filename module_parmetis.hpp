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

struct submesh{
  idx_t submeshid;
  idx_t nelems, nnodes, nboundaryelems;
  idx_t esize;

  std::vector<std::vector<real_t> > extents;
  std::set<idx_t> potentialneighbors;

  std::vector<std::vector<idx_t> > elems;
  std::vector<std::vector<idx_t> > neighbors;
  std::vector<std::vector<idx_t> > boundaryelems;
  std::vector<std::vector<real_t> > nodes;

  //Numbering maps for nodes and elements
  std::map<idx_t,idx_t> elems_gtl;
  std::map<idx_t,idx_t> elems_ltg;
  std::map<idx_t,idx_t> nodes_gtl;
  std::map<idx_t,idx_t> nodes_ltg;
  std::unordered_set<idx_t> nodesindex; 
  /* std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> edges; */
  /* std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> faces; */

  void Addelems(idx_t *elems, idx_t offset, idx_t nelems, idx_t esize);
  void buildEdges(std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges);
  void buildFaces(std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces);
  void buildElemConnectivity2D(std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges);
  void buildElemConnectivity3D(std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces);
  /* void buildEdges(); */
  /* void buildFaces(); */
  /* void buildElemConnectivity2D(); */
  /* void buildElemConnectivity3D(); */
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

  //Compute the list of partition to ba handled by me as parts
  idx_t firstsubmeshid;
  idx_t q=nsubmeshes/nprocs;
  idx_t r=nsubmeshes-q*nprocs;

  ownerofsubmesh.resize(nsubmeshes);

  idx_t nsubmeshesofp;
  for(idx_t p=0; p<nprocs; p++){
    if(me<r){
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

//ParallelReadMesh refactored from ParallelReadMesh in parmetis/programs/io.c
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
      if(ownerofsubmesh[i]!=me){//Send to idPartitions[i]
        nSendRequest+=1;
      }
    }
  } 

  idx_t nRecvRequest(0), maxRecv(0);
  for(idx_t i=0;i<ntotalsubmeshes;i++){
      if(ownerofsubmesh[i]==me){
        for(idx_t p=0;p<nprocs;p++){
          if(gsubmeshessizes[p*ntotalsubmeshes+i]>0 and p!=me){//Recv from p and add to my partition
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
        MPI_Isend((void *)submesheselems[i], (esize+1)*submeshessizes[i], IDX_T, ownerofsubmesh[i], 0, comm, &requestSendPtr[k]);
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
            MPI_Irecv((void *)(tmpElems+k*maxRecv*(esize+1)), (esize+1)*gsubmeshessizes[p*ntotalsubmeshes+i], IDX_T, p, MPI_ANY_TAG, comm, &requestRecvPtr[k]);
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

  delete [] tmpElems;
  delete [] requestSendPtr;
  delete [] requestRecvPtr;
  delete [] statSend;
  delete [] statRecv;
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

void writeVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
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
      for(idx_t p=0; p<esize; p++){
        outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].elems[i][p]] << " ";
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

    outFile << std::endl;
    outFile << "CELLS " << nboundaryelems << " " << nboundaryelems*(esize+1) << std::endl;
    for(idx_t i=0; i<nboundaryelems; i++){
      outFile << esize << " ";
      for(idx_t p=0; p<esize; p++){
        outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].boundaryelems[i][p]] << " ";
      }
      outFile << std::endl;
    }

    outFile << std::endl;
    outFile << "CELL_TYPES " << nboundaryelems << std::endl;
    idx_t cellType;
    if(dim==3) cellType=10; //tetrahedrons
    if(dim==2 and esize==3) cellType=5; //triangles
    if(dim==2 and esize==4) cellType=9; //quadrangles
    for(idx_t i=0; i<nboundaryelems; i++){
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

void Buildconnectivity(std::vector<submesh> &submeshesowned, idx_t dimension){
  if(dimension==2){
    for(idx_t k=0;k<submeshesowned.size();k++){
      std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> edges;
      submeshesowned[k].neighbors.resize(submeshesowned[k].nelems, std::vector<idx_t>(submeshesowned[k].esize,-1));
      submeshesowned[k].buildEdges(edges);
      submeshesowned[k].buildElemConnectivity2D(edges);
      /* submeshesowned[k].buildEdges(); */
      /* submeshesowned[k].buildElemConnectivity2D(); */
    }
  }
  if(dimension==3){
    for(idx_t k=0;k<submeshesowned.size();k++){
      std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> faces;
      submeshesowned[k].neighbors.resize(submeshesowned[k].nelems, std::vector<idx_t>(submeshesowned[k].esize,-1));

      auto t1 = std::chrono::high_resolution_clock::now();
      submeshesowned[k].buildFaces(faces);
      /* submeshesowned[k].buildFaces(); */
      auto t2 = std::chrono::high_resolution_clock::now();

      auto t3 = std::chrono::high_resolution_clock::now();
      submeshesowned[k].buildElemConnectivity3D(faces);
      /* submeshesowned[k].buildElemConnectivity3D(); */
      auto t4 = std::chrono::high_resolution_clock::now();
      auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
      auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count();
      /* std::cout << duration1/100000 << " " << duration2/100000 << "\n"; */
    }
  }
}

/* void submesh::buildEdges(){ */
void submesh::buildEdges(std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges){
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


void submesh::buildFaces(std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces){
/* void submesh::buildFaces(){ */
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

void submesh::buildElemConnectivity2D(std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges){
/* void submesh::buildElemConnectivity2D(){ */
  std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash>:: iterator it = edges.begin();

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
          neighbors[elem2][s] = elem2;
        }
      }
    }
    it++;
  }
}

void submesh::buildElemConnectivity3D(std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces){
/* void submesh::buildElemConnectivity3D(){ */
  std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash>:: iterator it = faces.begin();

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

void Findboundaryfromconnectivity(std::vector<submesh> &submeshesowned){

  for(idx_t k=0;k<submeshesowned.size();k++){
    //Count how many boundary elems there are, then allocate memory and fill
    idx_t nboundaryelems=0;
    for(idx_t i=0;i<submeshesowned[k].nelems;i++){
      if(submeshesowned[k].isBoundary(i)){
        nboundaryelems++;
      }
    }  

    submeshesowned[k].nboundaryelems = nboundaryelems;
    submeshesowned[k].boundaryelems.resize(nboundaryelems, std::vector<idx_t>(submeshesowned[k].esize,-1));

    idx_t iloc=0;
    for(idx_t i=0;i<submeshesowned[k].nelems;i++){
      if(submeshesowned[k].isBoundary(i)){
        for(idx_t j=0;j<submeshesowned[k].esize;j++){
          submeshesowned[k].boundaryelems[iloc][j] = submeshesowned[k].elems[i][j];
        }
        iloc++;
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
  MPI_Allreduce(MPI_IN_PLACE, (void *)allextents, nsubmeshes*3*2, REAL_T, MPI_SUM, comm);

  bool interx, intery, interz;
  //Compute potential neighbours
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(idx_t p=0; p<nsubmeshes; p++){
      interx = ((allextents[submeshesowned[k].submeshid*6+1]>allextents[p*6]) and
        (allextents[submeshesowned[k].submeshid*6+1]<allextents[p*6+1]));
      interx = interx or ((allextents[submeshesowned[k].submeshid*6]>allextents[p*6]) and
        (allextents[submeshesowned[k].submeshid*6]<allextents[p*6+1]));
      intery = ((allextents[submeshesowned[k].submeshid*6+3]>allextents[p*6+2]) and
        (allextents[submeshesowned[k].submeshid*6+3]<allextents[p*6+3]));
      intery = intery or ((allextents[submeshesowned[k].submeshid*6+2]>allextents[p*6+2]) and
        (allextents[submeshesowned[k].submeshid*6+2]<allextents[p*6+3]));
      interz = ((allextents[submeshesowned[k].submeshid*6+5]>allextents[p*6+4]) and
        (allextents[submeshesowned[k].submeshid*6+5]<allextents[p*6+5]));
      interz = interx or ((allextents[submeshesowned[k].submeshid*6+4]>allextents[p*6+4]) and
        (allextents[submeshesowned[k].submeshid*6+4]<allextents[p*6+5]));
      if((interx and intery) and interz){
        submeshesowned[k].potentialneighbors.insert(p);
      }
    }
  }

  std::set<idx_t>::iterator it;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it = submeshesowned[k].potentialneighbors.begin(); it!=submeshesowned[k].potentialneighbors.end(); it++){
      std::cout << submeshesowned[k].submeshid << " " << *it << std::endl;
    }
  }
}
