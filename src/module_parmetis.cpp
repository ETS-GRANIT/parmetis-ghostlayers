#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <algorithm>
#include <assert.h>
#include <chrono>
#include <thread>

#include "module_parmetis.hpp"

//Global variable of potential neighbors inner boundary to reduce memory usage
std::map<idx_t, potential_neighbors_boundary> g_potentialneighbors;

void submesh::computemyextents(){
  extents.resize(3,std::vector<real_t>(2, 0.));

  extents[0][0] = nodes[0*3+0];
  extents[0][1] = nodes[0*3+0];

  extents[1][0] = nodes[0*3+1];
  extents[1][1] = nodes[0*3+1];

  extents[2][0] = nodes[0*3+2];
  extents[2][1] = nodes[0*3+2];

  for(int i=1;i<get_nnodes();i++){
    if(nodes[i*3+0] < extents[0][0]) extents[0][0] = nodes[i*3+0];
    if(nodes[i*3+0] > extents[0][1]) extents[0][1] = nodes[i*3+0];

    if(nodes[i*3+1] < extents[1][0]) extents[1][0] = nodes[i*3+1];
    if(nodes[i*3+1] > extents[1][1]) extents[1][1] = nodes[i*3+1];

    if(nodes[i*3+2] < extents[2][0]) extents[2][0] = nodes[i*3+2];
    if(nodes[i*3+2] > extents[2][1]) extents[2][1] = nodes[i*3+2];
  }

}

bool submesh::isBoundary(idx_t ielem){
  for(idx_t i=0; i<esize; i++){
    if(neighbors[ielem*esize+i] == -1) return(1);
  }
  return(0);
}

void submesh::Addelems(idx_t *newelems, idx_t offset, idx_t nnewelems, idx_t newesize){
  esize = newesize; //should'nt change
  int nelems = get_nelems() + nnewelems;
  std::vector<idx_t> tmpelems(esize*nnewelems, 0);
  for(idx_t i=0;i<nnewelems;i++){
    elems_gtl.insert(std::make_pair(newelems[offset+i*(esize+1)],i+nelems-nnewelems));
    elems_ltg.insert(std::make_pair(i+nelems-nnewelems,newelems[offset+i*(esize+1)]));

    for(idx_t k=0; k<esize; k++){
      tmpelems[i*esize+k] = newelems[offset+i*(esize+1)+k+1];
      nodesindex.insert(tmpelems[i*esize+k]);
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

void buildBoundaryEdges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::set<idx_t> &elems_set){

  idx_t nelems=elems.size()/esize;

  std::set<idx_t>::iterator it=elems_set.begin();

  std::pair<idx_t, idx_t> dummyPair;
  idx_t sp;
  for(it=elems_set.begin();it!=elems_set.end();it++){
    for(idx_t s=0; s<esize; s++){
      sp = (s+1)%esize;
      dummyPair = std::make_pair(std::min(elems[*it*esize+s], elems[*it*esize+sp]), std::max(elems[*it*esize+s],elems[*it*esize+sp]));
      if(edges.find(dummyPair) != edges.end()){
        edges[dummyPair] = std::make_pair(edges[dummyPair].first, *it);
      }
      else{
        edges.insert({dummyPair, std::make_pair(*it, -1)});
      }
    }
  }
}

void buildEdges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges){

  idx_t nelems=elems.size()/esize;

  std::pair<idx_t, idx_t> dummyPair;
  idx_t sp;
  for(idx_t i=0;i<nelems;i++){
    for(idx_t s=0; s<esize; s++){
      sp = (s+1)%esize;
      dummyPair = std::make_pair(std::min(elems[i*esize+s], elems[i*esize+sp]), std::max(elems[i*esize+s],elems[i*esize+sp]));
      if(edges.find(dummyPair) != edges.end()){
        edges[dummyPair] = std::make_pair(edges[dummyPair].first, i);
      }
      else{
        edges.insert({dummyPair, std::make_pair(i, -1)});
      }
    }
  }
}

void buildBoundaryFaces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::set<idx_t> &elems_set){
  std::set<idx_t>::iterator it;
  idx_t nelems=elems.size()/esize;
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
        dummy[k] = elems[*it*esize+sp[k]];
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

void buildFaces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces){
  idx_t nelems=elems.size()/esize;
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
        dummy[k] = elems[i*esize+sp[k]];
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

void buildElemConnectivity2D(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::vector<idx_t> &neighbors){
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
        if(((elems[elem1*esize+s] == s1) and (elems[elem1*esize+sp] == s2)) or 
            ((elems[elem1*esize+s] == s2) and (elems[elem1*esize+sp] == s1))){
          neighbors[elem1*esize+s] = elem2;
        }
        if(((elems[elem2*esize+s] == s1) and (elems[elem2*esize+sp] == s2)) or 
            ((elems[elem2*esize+s] == s2) and (elems[elem2*esize+sp] == s1))){
          neighbors[elem2*esize+s] = elem1;
        }
      }
    }
    it++;
  }
}

void buildElemConnectivity3D(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::vector<idx_t> &neighbors){
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

        faceVert[0] = elems[elem1*esize+s];
        for(idx_t k=1; k<faceSize; k++){
          faceVert[k] = elems[elem1*esize+(s+k)%esize];
        }

        cond1 = true;
        for(idx_t k=0;k<faceSize;k++){
          cond1 = ((std::find(faceVert.begin(), faceVert.end(), it->first[k])!=faceVert.end()) and cond1);
        }
        if(cond1) neighbors[elem1*esize+s] = elem2;

        faceVert[0] = elems[elem2*esize+s];
        for(idx_t k=1; k<faceSize; k++){
          faceVert[k] = elems[elem2*esize+(s+k)%esize];
        }

        cond1 = true;
        for(idx_t k=0;k<faceSize;k++){
          cond1 = ((std::find(faceVert.begin(), faceVert.end(), it->first[k])!=faceVert.end()) and cond1);
        }
        if(cond1) neighbors[elem2*esize+s] = elem1;
      }
    }
    it++;
  }
}

void Buildconnectivity(std::vector<submesh> &submeshesowned, idx_t dimension){
  if(dimension==2){
    for(idx_t k=0;k<submeshesowned.size();k++){
      std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> edges;
      submeshesowned[k].neighbors.resize(submeshesowned[k].get_nelems()*submeshesowned[k].esize, -1);
      buildEdges(submeshesowned[k].esize, submeshesowned[k].elems, edges);
      buildElemConnectivity2D(submeshesowned[k].esize,submeshesowned[k].elems, edges, submeshesowned[k].neighbors);
    }
  }
  if(dimension==3){
    for(idx_t k=0;k<submeshesowned.size();k++){
      std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> faces;
      submeshesowned[k].neighbors.resize(submeshesowned[k].get_nelems()*submeshesowned[k].esize, -1);
      buildFaces(submeshesowned[k].esize,submeshesowned[k].elems, faces);
      buildElemConnectivity3D(submeshesowned[k].esize,submeshesowned[k].elems, faces, submeshesowned[k].neighbors);
    }
  }
}




void Findboundaryfromconnectivity(std::vector<submesh> &submeshesowned, idx_t method, idx_t numlayers){

  for(idx_t k=0;k<submeshesowned.size();k++){
    for(idx_t i=0;i<submeshesowned[k].get_nelems();i++){
      if(submeshesowned[k].isBoundary(i)){
        submeshesowned[k].boundaryelems.insert(i);

        /* if(method==1){ */
        //Add boundarynodes
        for(idx_t j=0; j<submeshesowned[k].esize; j++){
          if(submeshesowned[k].get_neighbors(i,j) == -1){
            submeshesowned[k].boundarynodes.insert(submeshesowned[k].get_elems(i,j));
            submeshesowned[k].boundarynodes.insert(submeshesowned[k].get_elems(i,(j+1)%submeshesowned[k].esize));
            submeshesowned[k].l1_boundarynodes.insert(submeshesowned[k].get_elems(i,j));
            submeshesowned[k].l1_boundarynodes.insert(submeshesowned[k].get_elems(i,(j+1)%submeshesowned[k].esize));
          }
        }
        /* } */
      }
    }  

    //If method=1 add cells that have only a point in the boundary points
    if(method==1){
      for(idx_t i=0;i<submeshesowned[k].get_nelems();i++){
        for(idx_t j=0;j<submeshesowned[k].esize;j++){
          if(submeshesowned[k].boundarynodes.find(submeshesowned[k].get_elems(i,j))!=submeshesowned[k].boundarynodes.end()){
            submeshesowned[k].boundaryelems.insert(i);
          }
        }
      }

      for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
          it!=submeshesowned[k].boundaryelems.end();it++){
        for(idx_t p=0; p<submeshesowned[k].esize; p++){
          submeshesowned[k].boundarynodes.insert(submeshesowned[k].get_elems(*it,p));
        }
      }
    }

    for(idx_t l=1; l<numlayers; l++){
      std::set<idx_t> tmpnewboundaryelems = submeshesowned[k].boundaryelems;
      std::set<idx_t>::iterator it;
      for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
        for(idx_t j=0;j<submeshesowned[k].esize;j++){
          if(submeshesowned[k].get_neighbors(*it,j) != -1){
            submeshesowned[k].boundaryelems.insert(submeshesowned[k].get_neighbors(*it,j));

            /* if(method==1){ */
            //Add boundarynodes
            for(idx_t j=0; j<submeshesowned[k].esize; j++){
              submeshesowned[k].boundarynodes.insert(submeshesowned[k].get_elems(*it,j));
            }
            /* } */

          }
        } 
      }

      if(method==1){
        for(idx_t i=0;i<submeshesowned[k].get_nelems();i++){
          for(idx_t j=0;j<submeshesowned[k].esize;j++){
            if(submeshesowned[k].boundarynodes.find(submeshesowned[k].get_elems(i,j))!=submeshesowned[k].boundarynodes.end()){
              submeshesowned[k].boundaryelems.insert(i);
            }
          }
        }
        for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
            it!=submeshesowned[k].boundaryelems.end();it++){
          for(idx_t p=0; p<submeshesowned[k].esize; p++){
            submeshesowned[k].boundarynodes.insert(submeshesowned[k].get_elems(*it,p));
          }
        }
      }
    }

    for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
        it!=submeshesowned[k].boundaryelems.end();it++){
      for(idx_t p=0; p<submeshesowned[k].esize; p++){
        submeshesowned[k].boundarynodes.insert(submeshesowned[k].get_elems(*it,p));
      }
    }

  }
}

void FindNodesElemsSendRecv(std::vector<submesh> &submeshesowned, idx_t dimension, idx_t method, idx_t numlayers){

  int esize = submeshesowned[0].esize;

  for(idx_t k; k<submeshesowned.size(); k++){
    //Build my boundary connectivity
    if(dimension==2){
      buildBoundaryEdges(esize,submeshesowned[k].elems, submeshesowned[k].boundary_edges, submeshesowned[k].boundaryelems);

      std::set<idx_t>::iterator iter;
      for(iter=submeshesowned[k].potentialneighbors.begin(); 
          iter!= submeshesowned[k].potentialneighbors.end(); iter++){

        if( not g_potentialneighbors[*iter].already_computed() ){

          g_potentialneighbors[*iter].neighbors.resize(g_potentialneighbors[*iter].get_nelems()*esize,-1);
          buildEdges(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].edges);
          buildElemConnectivity2D(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].edges, g_potentialneighbors[*iter].neighbors);
        }


        std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash>::iterator it_edges;
        for(it_edges=g_potentialneighbors[*iter].edges.begin();it_edges!=g_potentialneighbors[*iter].edges.end();it_edges++){
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
          for(idx_t i=0;i<g_potentialneighbors[*iter].get_nelems();i++){
            for(idx_t j=0;j<esize;j++){
              if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=submeshesowned[k].nodestorecv[*iter].end()){
                submeshesowned[k].elemstorecv[*iter].insert(i);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
              it!=submeshesowned[k].elemstorecv[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
            }
          }
        }

        //SEND If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
              it!=submeshesowned[k].boundaryelems.end();
              it++){
            for(idx_t j=0;j<esize;j++){
              if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].get_elems(*it,j))!=submeshesowned[k].nodestosend[*iter].end()){
                submeshesowned[k].elemstosend[*iter].insert(*it);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
              it!=submeshesowned[k].elemstosend[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].get_elems(*it,p));
            }
          }
        }

        //Add all the layers
        for(idx_t l=1; l<numlayers; l++){
          //Add to send list
          std::set<idx_t> tmpnewboundaryelems = submeshesowned[k].elemstosend[*iter];
          std::set<idx_t>::iterator it;
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(submeshesowned[k].get_neighbors(*it,j) != -1){
                submeshesowned[k].elemstosend[*iter].insert(submeshesowned[k].get_neighbors(*it,j));
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
                it!=submeshesowned[k].boundaryelems.end();
                it++){
              for(idx_t j=0;j<esize;j++){
                if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].get_elems(*it,j))!=submeshesowned[k].nodestosend[*iter].end()){
                  submeshesowned[k].elemstosend[*iter].insert(*it);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
                it!=submeshesowned[k].elemstosend[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].get_elems(*it,p));
              }
            }
          }

          //Add to recv list
          tmpnewboundaryelems.clear();
          tmpnewboundaryelems = submeshesowned[k].elemstorecv[*iter];
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(g_potentialneighbors[*iter].neighbors[*it*esize+j] != -1){
                submeshesowned[k].elemstorecv[*iter].insert(g_potentialneighbors[*iter].neighbors[*it*esize+j]);
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(idx_t i=0; i<g_potentialneighbors[*iter].get_nelems(); i++){
              for(idx_t j=0;j<esize;j++){
                if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=submeshesowned[k].nodestorecv[*iter].end()){
                  submeshesowned[k].elemstorecv[*iter].insert(i);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
                it!=submeshesowned[k].elemstorecv[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
              }
            }
          }


        }

        for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
            it!=submeshesowned[k].elemstorecv[*iter].end();it++){
          for(idx_t p=0; p<esize; p++){
            submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
          }
        }


      }
    }
    if(dimension==3){
      buildBoundaryFaces(esize,submeshesowned[k].elems, submeshesowned[k].boundary_faces, submeshesowned[k].boundaryelems);
      idx_t faceSize=3;//Tetrahedron

      std::set<idx_t>::iterator iter;
      for(iter=submeshesowned[k].potentialneighbors.begin(); 
          iter!= submeshesowned[k].potentialneighbors.end(); iter++){

        if( not g_potentialneighbors[*iter].already_computed() ){
          g_potentialneighbors[*iter].neighbors.resize(g_potentialneighbors[*iter].get_nelems()*esize,-1);

          buildFaces(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].faces);
          buildElemConnectivity3D(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].faces, g_potentialneighbors[*iter].neighbors);
        }

        std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash>::iterator it_faces;

        for(it_faces=g_potentialneighbors[*iter].faces.begin();it_faces!=g_potentialneighbors[*iter].faces.end();it_faces++){
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
          for(idx_t i=0;i<g_potentialneighbors[*iter].get_nelems();i++){
            for(idx_t j=0;j<esize;j++){
              if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=submeshesowned[k].nodestorecv[*iter].end()){
                submeshesowned[k].elemstorecv[*iter].insert(i);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
              it!=submeshesowned[k].elemstorecv[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
            }
          }
        }

        //SEND If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
              it!=submeshesowned[k].boundaryelems.end();
              it++){
            for(idx_t j=0;j<esize;j++){
              if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].get_elems(*it,j))!=submeshesowned[k].nodestosend[*iter].end()){
                submeshesowned[k].elemstosend[*iter].insert(*it);
              }
            }
          }
          for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
              it!=submeshesowned[k].elemstosend[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].get_elems(*it,p));
            }
          }
        }

        //Add all the layers
        for(idx_t l=1; l<numlayers; l++){
          //Add to send list
          std::set<idx_t> tmpnewboundaryelems = submeshesowned[k].elemstosend[*iter];
          std::set<idx_t>::iterator it;
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(submeshesowned[k].get_neighbors(*it,j) != -1){
                submeshesowned[k].elemstosend[*iter].insert(submeshesowned[k].get_neighbors(*it,j));
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(std::set<idx_t>::iterator it=submeshesowned[k].boundaryelems.begin();
                it!=submeshesowned[k].boundaryelems.end();
                it++){
              for(idx_t j=0;j<esize;j++){
                if(submeshesowned[k].nodestosend[*iter].find(submeshesowned[k].get_elems(*it,j))!=submeshesowned[k].nodestosend[*iter].end()){
                  submeshesowned[k].elemstosend[*iter].insert(*it);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
                it!=submeshesowned[k].elemstosend[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                submeshesowned[k].nodestosend[*iter].insert(submeshesowned[k].get_elems(*it,p));
              }
            }
          }

          //Add to recv list
          tmpnewboundaryelems.clear();
          tmpnewboundaryelems = submeshesowned[k].elemstorecv[*iter];
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(g_potentialneighbors[*iter].neighbors[*it*esize+j] != -1){
                submeshesowned[k].elemstorecv[*iter].insert(g_potentialneighbors[*iter].neighbors[*it*esize+j]);
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(idx_t i=0; i<g_potentialneighbors[*iter].get_nelems(); i++){
              for(idx_t j=0;j<esize;j++){
                if(submeshesowned[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=submeshesowned[k].nodestorecv[*iter].end()){
                  submeshesowned[k].elemstorecv[*iter].insert(i);
                }
              }
            }
            for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
                it!=submeshesowned[k].elemstorecv[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
              }
            }
          }


        }

        for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
            it!=submeshesowned[k].elemstorecv[*iter].end();it++){
          for(idx_t p=0; p<esize; p++){
            submeshesowned[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
          }
        }
      }

    }
  }

}

void ShareboundaryFixPotentialNeighbors(std::vector<submesh> &submeshesowned, std::vector<idx_t> &ownerofsubmesh, MPI_Comm comm){

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
    maxnodestosend = std::max(maxnodestosend, (idx_t) submeshesowned[k].l1_boundarynodes.size());
  }

  idx_t *nodestosend = new idx_t[maxnodestosend*submeshesowned.size()];

  //Compute total number of messages/request
  idx_t ntotreq_send=0;
  idx_t ntotreq_recv=0;
  idx_t nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors_extents.begin(); it!= submeshesowned[k].potentialneighbors_extents.end(); it++){
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
    for(it=submeshesowned[k].l1_boundarynodes.begin(); it!=submeshesowned[k].l1_boundarynodes.end(); it++){
      nodestosend[k*maxnodestosend+iloc] = (idx_t) *it;
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
    for(it=submeshesowned[k].potentialneighbors_extents.begin(); it!= submeshesowned[k].potentialneighbors_extents.end(); it++){
      if(msgs_send[ownerofsubmesh[*it]].count(submeshesowned[k].submeshid)==0){
        if(ownerofsubmesh[*it]!=me){
          MPI_Isend(&nodestosend[maxnodestosend*k], submeshesowned[k].l1_boundarynodes.size(), IDX_T, ownerofsubmesh[*it], (*it)*ownerofsubmesh.size()+submeshesowned[k].submeshid, comm, &requestSendPtr[nreq]);
          nreq++;
          msgs_send[ownerofsubmesh[*it]].insert(submeshesowned[k].submeshid);
        }
        else{
          //Check for node in common, if found add to potentialnieghbors
          idx_t k2 = submeshlocid[*it];
          for(idx_t i=0; i<submeshesowned[k].l1_boundarynodes.size(); i++){
            if(submeshesowned[k2].l1_boundarynodes.count(nodestosend[maxnodestosend*k+i]) != 0){
              /* submeshesowned[k].potentialneighbors.insert(*it); */
              submeshesowned[k2].potentialneighbors.insert(submeshesowned[k].submeshid);
              break;
            }
          }

        }
      }
    }
  }

  //Probe to get sizes
  idx_t *nnodestorecv = new idx_t[ntotreq_recv];
  idx_t maxnodestorecv=0;

  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors_extents.begin(); it!= submeshesowned[k].potentialneighbors_extents.end(); it++){
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        MPI_Probe(ownerofsubmesh[*it], (*it)+ownerofsubmesh.size()*submeshesowned[k].submeshid, comm, &stat);
        MPI_Get_count(&stat,IDX_T,&nnodestorecv[nreq]);
        maxnodestorecv = std::max(maxnodestorecv,nnodestorecv[nreq]);
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        nreq++;
      }
    }
  }

  idx_t *nodestorecv = new idx_t[maxnodestorecv*ntotreq_recv];

  //Clear msgsend msgs recv
  for(idx_t pr=0;pr<nprocs;pr++){
    msgs_send[pr].clear();
    msgs_recv[pr].clear();
  }

  std::vector<std::map<idx_t, idx_t> > msgs_req_recv(nprocs);

  //IRecv nodes vector
  nreq=0;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(it=submeshesowned[k].potentialneighbors_extents.begin(); it!= submeshesowned[k].potentialneighbors_extents.end(); it++){
      if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==0){
        MPI_Recv(&nodestorecv[nreq*maxnodestorecv], nnodestorecv[nreq], IDX_T, ownerofsubmesh[*it], (*it)+ownerofsubmesh.size()*submeshesowned[k].submeshid, comm, &stat);

        //Check for node in common, if found add to potentialnieghbors
        for(idx_t i=0; i<nnodestorecv[nreq]; i++){
          if(submeshesowned[k].l1_boundarynodes.count(nodestorecv[nreq*maxnodestorecv+i]) != 0){
            submeshesowned[k].potentialneighbors.insert(*it);
            break;
          }
        }
        msgs_recv[ownerofsubmesh[*it]].insert(*it);
        msgs_req_recv[ownerofsubmesh[*it]].insert({*it,nreq});
        nreq++;
      }
      else if(ownerofsubmesh[*it]!=me && msgs_recv[ownerofsubmesh[*it]].count(*it)==1){
        //Find what req the msg came with
        idx_t loc_req=msgs_req_recv[ownerofsubmesh[*it]][*it];
        //Check for node in common, if found add to potentialnieghbors
        for(idx_t i=0; i<nnodestorecv[loc_req]; i++){
          if(submeshesowned[k].l1_boundarynodes.count(nodestorecv[loc_req*maxnodestorecv+i]) != 0){
            submeshesowned[k].potentialneighbors.insert(*it);
            break;
          }
        }
      }
    }
  }

  MPI_Waitall(ntotreq_send, requestSendPtr, statSend);

  delete [] nodestosend;
  delete [] nodestorecv;
  delete [] nnodestorecv;
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
      nodestosend[4*k*maxnodestosend+4*iloc+1] = submeshesowned[k].get_nodes(idum,0);
      nodestosend[4*k*maxnodestosend+4*iloc+2] = submeshesowned[k].get_nodes(idum,1);
      nodestosend[4*k*maxnodestosend+4*iloc+3] = submeshesowned[k].get_nodes(idum,2);
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
          pnbound.nodes.resize(submeshesowned[k].boundarynodes.size()*3,0.);
          for(idx_t i=0;i<submeshesowned[k].boundarynodes.size();i++){
            pnbound.nodes_ltg.insert(std::make_pair(i,nodestosend[4*maxnodestosend*k+i*4]));
            pnbound.nodes_gtl.insert(std::make_pair(nodestosend[4*maxnodestosend*k+i*4],i));
            pnbound.get_nodes(i,0) = nodestosend[4*maxnodestosend*k+i*4+1];
            pnbound.get_nodes(i,1) = nodestosend[4*maxnodestosend*k+i*4+2];
            pnbound.get_nodes(i,2) = nodestosend[4*maxnodestosend*k+i*4+3];
          }
          g_potentialneighbors.insert(std::make_pair(submeshesowned[k].submeshid,pnbound));
          g_potentialneighbors[submeshesowned[k].submeshid].esize = submeshesowned[k].esize;
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
        pnbound.nodes.resize(nnodestorecv[nreq]*3,0.);
        for(idx_t i=0;i<nnodestorecv[nreq];i++){
          pnbound.nodes_ltg.insert(std::make_pair(i,(idx_t) nodestorecv[4*maxnodestorecv*nreq+i*4]));
          pnbound.nodes_gtl.insert(std::make_pair((idx_t) nodestorecv[4*maxnodestorecv*nreq+i*4],i));
          pnbound.get_nodes(i,0) = nodestorecv[4*maxnodestorecv*nreq+i*4+1];
          pnbound.get_nodes(i,1) = nodestorecv[4*maxnodestorecv*nreq+i*4+2];
          pnbound.get_nodes(i,2) = nodestorecv[4*maxnodestorecv*nreq+i*4+3];
        }
        g_potentialneighbors.insert(std::make_pair(*it,pnbound));
        g_potentialneighbors[*it].esize = submeshesowned[k].esize;
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
        elemstosend[(1+esize)*(maxelemstosend*k+iloc)+j+1] = submeshesowned[k].get_elems(*it,j);
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
          g_potentialneighbors[submeshesowned[k].submeshid].esize = esize;
          g_potentialneighbors[submeshesowned[k].submeshid].elems.resize(submeshesowned[k].boundaryelems.size()*esize,0);
          for(idx_t i=0;i<submeshesowned[k].boundaryelems.size();i++){
            g_potentialneighbors[submeshesowned[k].submeshid].elems_ltg.insert(std::make_pair(i,elemstosend[(1+esize)*(maxelemstosend*k+i)]));
            g_potentialneighbors[submeshesowned[k].submeshid].elems_gtl.insert(std::make_pair(elemstosend[(1+esize)*(maxelemstosend*k+i)],i));
            for(idx_t j=0;j<esize;j++){
              g_potentialneighbors[submeshesowned[k].submeshid].get_elems(i,j) = elemstosend[(1+esize)*(maxelemstosend*k+i)+j+1];
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
        g_potentialneighbors[*it].esize = esize;
        g_potentialneighbors[*it].elems.resize(nelemstorecv[nreq]*esize,0);
        for(idx_t i=0;i<nelemstorecv[nreq];i++){
          g_potentialneighbors[*it].elems_ltg.insert(std::make_pair(i,elemstorecv[(1+esize)*(maxelemstorecv*nreq+i)]));
          g_potentialneighbors[*it].elems_gtl.insert(std::make_pair(elemstorecv[(1+esize)*(maxelemstorecv*nreq+i)],i));
          for(idx_t j=0;j<esize;j++){
            g_potentialneighbors[*it].get_elems(i,j) = elemstorecv[(1+esize)*(maxelemstorecv*nreq+i)+j+1];
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
    submeshesowned[k].computemyextents();
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
          submeshesowned[k].potentialneighbors_extents.insert(p);
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
          submeshesowned[k].nodes_ltg.insert(std::make_pair(submeshesowned[k].get_nnodes(), *it));
          submeshesowned[k].nodes_gtl.insert(std::make_pair(*it, submeshesowned[k].get_nnodes()));
          for(int kk=0; kk<3; kk++){
            submeshesowned[k].nodes.push_back(g_potentialneighbors[*iter].get_nodes(g_potentialneighbors[*iter].nodes_gtl[*it],kk));
          }
        }
      }
      for(std::set<idx_t>::iterator it=submeshesowned[k].elemstorecv[*iter].begin();
          it!=submeshesowned[k].elemstorecv[*iter].end();
          it++){
        idx_t itglob = g_potentialneighbors[*iter].elems_ltg[*it];
        submeshesowned[k].elems_ltg.insert(std::make_pair(submeshesowned[k].get_nelems(), itglob));
        submeshesowned[k].elems_gtl.insert(std::make_pair(itglob, submeshesowned[k].get_nelems()));
        for(int kk=0; kk<submeshesowned[k].esize; kk++){
          submeshesowned[k].elems.push_back(g_potentialneighbors[*iter].get_elems(*it,kk));
        }
      }
    }
    /* submeshesowned[k].nnodes = submeshesowned[k].get_nnodes(); */
    /* submeshesowned[k].nelems = submeshesowned[k].get_nelems(); */

    //Renumbering
    submeshesowned[k].renumber_otn.resize(submeshesowned[k].get_nelems(), -1);
    submeshesowned[k].renumber_nto.resize(submeshesowned[k].get_nelems(), -1);
    idx_t itloc;
    idx_t ind=0;
    bool norsendorrecv;
    for(idx_t i=0; i<submeshesowned[k].get_nelems(); i++){
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
        submeshesowned[k].renumber_otn[i] = ind;
        submeshesowned[k].renumber_nto[ind] = i;
        ind++;
      }
    }

    for(std::set<idx_t>::iterator iter=submeshesowned[k].potentialneighbors.begin();
        iter!=submeshesowned[k].potentialneighbors.end();
        iter++){
      for(std::set<idx_t>::iterator it=submeshesowned[k].elemstosend[*iter].begin();
          it!=submeshesowned[k].elemstosend[*iter].end();
          it++){
        if(submeshesowned[k].renumber_otn[*it] == -1){
          submeshesowned[k].renumber_otn[*it] = ind;
          submeshesowned[k].renumber_nto[ind] = *it;
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
        submeshesowned[k].renumber_otn[itglobloc] = ind;
        submeshesowned[k].renumber_nto[ind] = itglobloc;
        ind++;
      }
    }
  }
}
