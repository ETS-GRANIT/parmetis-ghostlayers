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

  /* for(idx_t k=0; k<nsubmeshesowned; k++){ */
  /*   std::cout << me << " owns " << submeshesowned[k].submeshid << std::endl; */
  /* } */
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

    outFile << std::endl;
    outFile << "POINT_DATA " << nNodes << std::endl;
    outFile << "SCALARS b float" << std::endl;
    outFile << "LOOKUP_TABLE default" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << submeshesowned[k].nodes[i][2] << std::endl;
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
    /* std::cout << submeshesowned[k].submeshid << " " << submeshesowned[k].nelems << " " << ind << std::endl; */

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
    /* std::cout << submeshesowned[k].submeshid << " " << submeshesowned[k].nelems << " " << ind << std::endl; */

    for(std::set<idx_t>::iterator iter=submeshesowned[k].potentialneighbors.begin();
        iter!=submeshesowned[k].potentialneighbors.end();
        iter++){
      for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
          it!=submeshesowned[k].elemstorecv[*iter].end();
          it++){
        idx_t itglobloc = submeshesowned[k].elems_gtl[g_potentialneighbors[*iter].elems_ltg[*it]];
        submeshesowned[k].erenumber[itglobloc] = ind;
        submeshesowned[k].erenumber[ind] = itglobloc;
        ind++;
      }
    }
    /* std::cout << submeshesowned[k].submeshid << " " << submeshesowned[k].nelems << " " << ind <<  " " << ind-submeshesowned[k].nelems << std::endl; */


  }
}
