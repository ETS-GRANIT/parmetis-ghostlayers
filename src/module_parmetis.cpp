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

void partition::compute_extents(){
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

bool partition::is_boundary(idx_t ielem){
  for(idx_t i=0; i<esize; i++){
    if(neighbors[ielem*esize+i] == -1) return(1);
  }
  return(0);
}

void partition::add_elems(idx_t *newelems, idx_t offset, idx_t nnewelems, idx_t newesize){
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

void compute_partition_ownership(idx_t ntotparts, idx_t &ntotpartsowned, std::vector<partition> &parts, std::vector<idx_t> &ownerofpartition, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  //Compute the list of epartition to ba handled by me
  idx_t firstpartitionid;
  idx_t q=ntotparts/nprocs;
  idx_t r=ntotparts-q*nprocs;

  ownerofpartition.resize(ntotparts);

  idx_t ntotpartsofp;
  for(idx_t p=0; p<nprocs; p++){
    if(p<r){
      ntotpartsofp = q+1;
      firstpartitionid = p*(q+1);
    }
    else{
      ntotpartsofp = q;
      firstpartitionid = p*q + r;
    }

    for(idx_t k=0;k<ntotpartsofp;k++){
      ownerofpartition[k+firstpartitionid] = p;
    }
    if(me==p){
      ntotpartsowned = ntotpartsofp;
      parts.resize(ntotpartsowned);
      for(idx_t p=0;p<ntotpartsowned;p++){
        parts[p].partitionid = p+firstpartitionid; 
      }
    }
  }
}

void gather_partitions(idx_t*& elmdist, idx_t*& eind, idx_t*& epart, const idx_t esize, std::vector<partition> &parts, std::vector<idx_t> const ownerofpartition, MPI_Comm comm){

  std::stringstream line;
  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;
  MPI_Request *requestSendPtr, *requestRecvPtr;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  idx_t nelems = elmdist[me+1]-elmdist[me];

  //Generate parts
  std::map<idx_t, idx_t> iter;
  std::map<idx_t, idx_t*> partselems; //Needs to be properly freed
  std::map<idx_t, idx_t> partssizes;

  for(idx_t i=0;i<nelems;i++){
    if(partssizes.find(epart[i]) == partssizes.end()){
      partssizes.insert(std::make_pair(epart[i], 0));
    }
    partssizes[epart[i]] += 1;
  }

  for(idx_t i=0;i<nelems;i++){
    if(partselems.find(epart[i]) == partselems.end()){
      partselems.insert(std::make_pair(epart[i], new idx_t[((esize+1)*partssizes[epart[i]])]));
      iter.insert(std::make_pair(epart[i],0));
    }
    partselems[epart[i]][iter[epart[i]]*(esize+1)] = elmdist[me]+i;
    for(idx_t j=0;j<esize;j++){
      partselems[epart[i]][iter[epart[i]]*(esize+1)+j+1] = eind[esize*i+j];
    }
    iter[epart[i]] += 1;
  }

  //Communicate the size of each partition owned
  idx_t ntotalparts = ownerofpartition.size();
  idx_t *gpartssizes = new idx_t[(ntotalparts)*nprocs];

  for(idx_t i=0;i<ntotalparts*nprocs;i++){
    gpartssizes[i] = 0;
  }

  for(idx_t i=0;i<ntotalparts;i++){
    if(partselems.find(i) == partselems.end()){
      gpartssizes[me*ntotalparts+i] = 0;
    }
    else{
      gpartssizes[me*ntotalparts+i] = partssizes[i];
    }
  } 

  MPI_Allgather((void *)(gpartssizes+me*ntotalparts), ntotalparts, IDX_T, (void *)gpartssizes, ntotalparts, IDX_T, comm);

  std::map<idx_t, idx_t> partitionlocid;
  for(idx_t i=0;i<parts.size();i++){
    partitionlocid.insert(std::make_pair(parts[i].partitionid, i));
  }

  idx_t nSendRequest(0);
  for(idx_t i=0;i<ntotalparts;i++){
    if(gpartssizes[me*ntotalparts+i] > 0){
      if(ownerofpartition[i]!=me){
        nSendRequest+=1;
      }
    }
  } 

  idx_t nRecvRequest(0), maxRecv(0);
  for(idx_t i=0;i<ntotalparts;i++){
    if(ownerofpartition[i]==me){
      for(idx_t p=0;p<nprocs;p++){
        if(gpartssizes[p*ntotalparts+i]>0 and p!=me){
          nRecvRequest+=1;
          maxRecv=std::max((esize+1)*gpartssizes[p*ntotalparts+i],maxRecv);
        }
      }
    }
  }

  requestSendPtr = new MPI_Request[nSendRequest];
  requestRecvPtr = new MPI_Request[nRecvRequest];
  statSend = new MPI_Status[nSendRequest];
  statRecv = new MPI_Status[nRecvRequest];

  idx_t k=0;
  for(idx_t i=0;i<ntotalparts;i++){
    if(gpartssizes[me*ntotalparts+i] > 0){
      if(ownerofpartition[i]!=me){//Send to ownerofpartition[i]
        MPI_Isend((void *)partselems[i], (esize+1)*partssizes[i], IDX_T, ownerofpartition[i], me*ntotalparts+i, comm, &requestSendPtr[k]);
        k+=1;
      }
      else{//Add to my epartition
        parts[partitionlocid[i]].add_elems(partselems[i], 0, partssizes[i], esize);
      }
    }
  } 

  k=0;
  idx_t *tmpElems = new idx_t[(esize+1)*nRecvRequest*maxRecv];
  for(idx_t i=0;i<ntotalparts;i++){
    if(ownerofpartition[i]==me){
      for(idx_t p=0;p<nprocs;p++){
        if(gpartssizes[p*ntotalparts+i]>0 and p!=me){//Recv from p and add to my epartition
          MPI_Irecv((void *)(tmpElems+k*maxRecv*(esize+1)), (esize+1)*gpartssizes[p*ntotalparts+i], IDX_T, p, p*ntotalparts+i, comm, &requestRecvPtr[k]);
          k+=1;
        }
      }
    }
  }

  MPI_Waitall(nSendRequest, requestSendPtr, statSend);
  MPI_Waitall(nRecvRequest, requestRecvPtr, statRecv);

  k=0;
  for(idx_t i=0;i<ntotalparts;i++){
    if(ownerofpartition[i]==me){
      for(idx_t p=0;p<nprocs;p++){
        if(gpartssizes[p*ntotalparts+i]>0 and p!=me){//Recv from p and add to my epartition
          parts[partitionlocid[i]].add_elems(tmpElems, k*maxRecv*(esize+1), gpartssizes[p*ntotalparts+i], esize);
          k+=1;
        }
      }
    }
  }

  //Free of partselems
  std::map<idx_t, idx_t*>::iterator it;
  for(it=partselems.begin();it!=partselems.end();it++){
    delete [] it->second;
  }
}

void build_boundary_edges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::set<idx_t> &elems_set){

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

void build_edges(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges){

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

void build_faces(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces){
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

void build_elems_connectivity_2d(int esize, std::vector<idx_t> &elems, std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> &edges, std::vector<idx_t> &neighbors){
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

void build_elems_connectivity_3d(int esize, std::vector<idx_t> &elems, std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> &faces, std::vector<idx_t> &neighbors){
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

void build_connectivity(std::vector<partition> &parts, idx_t dimension){
  if(dimension==2){
    for(idx_t k=0;k<parts.size();k++){
      std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash> edges;
      parts[k].neighbors.resize(parts[k].get_nelems()*parts[k].esize, -1);
      build_edges(parts[k].esize, parts[k].elems, edges);
      build_elems_connectivity_2d(parts[k].esize,parts[k].elems, edges, parts[k].neighbors);
    }
  }
  if(dimension==3){
    for(idx_t k=0;k<parts.size();k++){
      std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash> faces;
      parts[k].neighbors.resize(parts[k].get_nelems()*parts[k].esize, -1);
      build_faces(parts[k].esize,parts[k].elems, faces);
      build_elems_connectivity_3d(parts[k].esize,parts[k].elems, faces, parts[k].neighbors);
    }
  }
}




void find_boundary_from_connectivity(std::vector<partition> &parts, idx_t method, idx_t numlayers){

  for(idx_t k=0;k<parts.size();k++){
    for(idx_t i=0;i<parts[k].get_nelems();i++){
      if(parts[k].is_boundary(i)){
        parts[k].boundaryelems.insert(i);

        //Add boundarynodes
        for(idx_t j=0; j<parts[k].esize; j++){
          if(parts[k].get_neighbors(i,j) == -1){
            parts[k].boundarynodes.insert(parts[k].get_elems(i,j));
            parts[k].boundarynodes.insert(parts[k].get_elems(i,(j+1)%parts[k].esize));
            parts[k].l1_boundarynodes.insert(parts[k].get_elems(i,j));
            parts[k].l1_boundarynodes.insert(parts[k].get_elems(i,(j+1)%parts[k].esize));
          }
        }
      }
    }  

    //If method=1 add cells that have only a point in the boundary points
    if(method==1){
      for(idx_t i=0;i<parts[k].get_nelems();i++){
        for(idx_t j=0;j<parts[k].esize;j++){
          if(parts[k].boundarynodes.find(parts[k].get_elems(i,j))!=parts[k].boundarynodes.end()){
            parts[k].boundaryelems.insert(i);
          }
        }
      }

      for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
          it!=parts[k].boundaryelems.end();it++){
        for(idx_t p=0; p<parts[k].esize; p++){
          parts[k].boundarynodes.insert(parts[k].get_elems(*it,p));
        }
      }
    }

    for(idx_t l=1; l<numlayers; l++){
      std::set<idx_t> tmpnewboundaryelems = parts[k].boundaryelems;
      std::set<idx_t>::iterator it;
      for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
        for(idx_t j=0;j<parts[k].esize;j++){
          if(parts[k].get_neighbors(*it,j) != -1){
            parts[k].boundaryelems.insert(parts[k].get_neighbors(*it,j));

            //Add boundarynodes
            for(idx_t j=0; j<parts[k].esize; j++){
              parts[k].boundarynodes.insert(parts[k].get_elems(*it,j));
            }

          }
        } 
      }

      if(method==1){
        for(idx_t i=0;i<parts[k].get_nelems();i++){
          for(idx_t j=0;j<parts[k].esize;j++){
            if(parts[k].boundarynodes.find(parts[k].get_elems(i,j))!=parts[k].boundarynodes.end()){
              parts[k].boundaryelems.insert(i);
            }
          }
        }
        for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
            it!=parts[k].boundaryelems.end();it++){
          for(idx_t p=0; p<parts[k].esize; p++){
            parts[k].boundarynodes.insert(parts[k].get_elems(*it,p));
          }
        }
      }
    }

    for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
        it!=parts[k].boundaryelems.end();it++){
      for(idx_t p=0; p<parts[k].esize; p++){
        parts[k].boundarynodes.insert(parts[k].get_elems(*it,p));
      }
    }

  }
}

void find_nodes_send_recv(std::vector<partition> &parts, idx_t dimension, idx_t method, idx_t numlayers){

  int esize = parts[0].esize;

  for(idx_t k; k<parts.size(); k++){
    //Build my boundary connectivity
    if(dimension==2){
      build_boundary_edges(esize,parts[k].elems, parts[k].boundary_edges, parts[k].boundaryelems);

      std::set<idx_t>::iterator iter;
      for(iter=parts[k].potentialneighbors.begin(); 
          iter!= parts[k].potentialneighbors.end(); iter++){

        if( not g_potentialneighbors[*iter].already_computed() ){

          g_potentialneighbors[*iter].neighbors.resize(g_potentialneighbors[*iter].get_nelems()*esize,-1);
          build_edges(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].edges);
          build_elems_connectivity_2d(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].edges, g_potentialneighbors[*iter].neighbors);
        }


        std::unordered_map<std::pair<idx_t,idx_t>, std::pair<idx_t,idx_t>, pair_idx_t_hash>::iterator it_edges;
        for(it_edges=g_potentialneighbors[*iter].edges.begin();it_edges!=g_potentialneighbors[*iter].edges.end();it_edges++){
          if(parts[k].boundary_edges.find(it_edges->first)!=parts[k].boundary_edges.end()){

            //Add elem to recv
            idx_t elem=it_edges->second.first;
            parts[k].elemstorecv[*iter].insert(elem);
            //Add boundarynodes to recv
            parts[k].nodestorecv[*iter].insert(it_edges->first.first);
            parts[k].nodestorecv[*iter].insert(it_edges->first.second);

            //Add elem to send
            elem=parts[k].boundary_edges[it_edges->first].first;
            parts[k].elemstosend[*iter].insert(elem);
            //Add boundarynodes to send
            parts[k].nodestosend[*iter].insert(it_edges->first.first);
            parts[k].nodestosend[*iter].insert(it_edges->first.second);


          }
        }

        //RECV If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(idx_t i=0;i<g_potentialneighbors[*iter].get_nelems();i++){
            for(idx_t j=0;j<esize;j++){
              if(parts[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=parts[k].nodestorecv[*iter].end()){
                parts[k].elemstorecv[*iter].insert(i);
              }
            }
          }
          for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
              it!=parts[k].elemstorecv[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
            }
          }
        }

        //SEND If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
              it!=parts[k].boundaryelems.end();
              it++){
            for(idx_t j=0;j<esize;j++){
              if(parts[k].nodestosend[*iter].find(parts[k].get_elems(*it,j))!=parts[k].nodestosend[*iter].end()){
                parts[k].elemstosend[*iter].insert(*it);
              }
            }
          }
          for(std::set<idx_t>::iterator it=parts[k].elemstosend[*iter].begin();
              it!=parts[k].elemstosend[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              parts[k].nodestosend[*iter].insert(parts[k].get_elems(*it,p));
            }
          }
        }

        //Add all the layers
        for(idx_t l=1; l<numlayers; l++){
          //Add to send list
          std::set<idx_t> tmpnewboundaryelems = parts[k].elemstosend[*iter];
          std::set<idx_t>::iterator it;
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(parts[k].get_neighbors(*it,j) != -1){
                parts[k].elemstosend[*iter].insert(parts[k].get_neighbors(*it,j));
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  parts[k].nodestosend[*iter].insert(parts[k].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
                it!=parts[k].boundaryelems.end();
                it++){
              for(idx_t j=0;j<esize;j++){
                if(parts[k].nodestosend[*iter].find(parts[k].get_elems(*it,j))!=parts[k].nodestosend[*iter].end()){
                  parts[k].elemstosend[*iter].insert(*it);
                }
              }
            }
            for(std::set<idx_t>::iterator it=parts[k].elemstosend[*iter].begin();
                it!=parts[k].elemstosend[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                parts[k].nodestosend[*iter].insert(parts[k].get_elems(*it,p));
              }
            }
          }

          //Add to recv list
          tmpnewboundaryelems.clear();
          tmpnewboundaryelems = parts[k].elemstorecv[*iter];
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(g_potentialneighbors[*iter].neighbors[*it*esize+j] != -1){
                parts[k].elemstorecv[*iter].insert(g_potentialneighbors[*iter].neighbors[*it*esize+j]);
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(idx_t i=0; i<g_potentialneighbors[*iter].get_nelems(); i++){
              for(idx_t j=0;j<esize;j++){
                if(parts[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=parts[k].nodestorecv[*iter].end()){
                  parts[k].elemstorecv[*iter].insert(i);
                }
              }
            }
            for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
                it!=parts[k].elemstorecv[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
              }
            }
          }


        }

        for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
            it!=parts[k].elemstorecv[*iter].end();it++){
          for(idx_t p=0; p<esize; p++){
            parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
          }
        }


      }
    }
    if(dimension==3){
      buildBoundaryFaces(esize,parts[k].elems, parts[k].boundary_faces, parts[k].boundaryelems);
      idx_t faceSize=3;//Tetrahedron

      std::set<idx_t>::iterator iter;
      for(iter=parts[k].potentialneighbors.begin(); 
          iter!= parts[k].potentialneighbors.end(); iter++){

        if( not g_potentialneighbors[*iter].already_computed() ){
          g_potentialneighbors[*iter].neighbors.resize(g_potentialneighbors[*iter].get_nelems()*esize,-1);

          build_faces(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].faces);
          build_elems_connectivity_3d(esize,g_potentialneighbors[*iter].elems, g_potentialneighbors[*iter].faces, g_potentialneighbors[*iter].neighbors);
        }

        std::unordered_map<std::vector<idx_t>, std::pair<idx_t,idx_t>, vector_idx_t_hash>::iterator it_faces;

        for(it_faces=g_potentialneighbors[*iter].faces.begin();it_faces!=g_potentialneighbors[*iter].faces.end();it_faces++){
          if(parts[k].boundary_faces.find(it_faces->first)!=parts[k].boundary_faces.end()){

            //Add elem to recv
            idx_t elem=it_faces->second.first;
            parts[k].elemstorecv[*iter].insert(elem);
            //Add boundarynodes to recv
            for(idx_t p=0; p<faceSize; p++){
              parts[k].nodestorecv[*iter].insert(it_faces->first[p]);
            }

            //Add elem to send
            elem=parts[k].boundary_faces[it_faces->first].first;
            parts[k].elemstosend[*iter].insert(elem);
            //Add boundarynodes to send
            for(idx_t p=0; p<faceSize; p++){
              parts[k].nodestosend[*iter].insert(it_faces->first[p]);
            }


          }
        }

        //RECV If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(idx_t i=0;i<g_potentialneighbors[*iter].get_nelems();i++){
            for(idx_t j=0;j<esize;j++){
              if(parts[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=parts[k].nodestorecv[*iter].end()){
                parts[k].elemstorecv[*iter].insert(i);
              }
            }
          }
          for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
              it!=parts[k].elemstorecv[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
            }
          }
        }

        //SEND If method=1 add cells that have only a point in the boundary points
        if(method==1){
          for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
              it!=parts[k].boundaryelems.end();
              it++){
            for(idx_t j=0;j<esize;j++){
              if(parts[k].nodestosend[*iter].find(parts[k].get_elems(*it,j))!=parts[k].nodestosend[*iter].end()){
                parts[k].elemstosend[*iter].insert(*it);
              }
            }
          }
          for(std::set<idx_t>::iterator it=parts[k].elemstosend[*iter].begin();
              it!=parts[k].elemstosend[*iter].end();
              it++){
            for(idx_t p=0; p<esize; p++){
              parts[k].nodestosend[*iter].insert(parts[k].get_elems(*it,p));
            }
          }
        }

        //Add all the layers
        for(idx_t l=1; l<numlayers; l++){
          //Add to send list
          std::set<idx_t> tmpnewboundaryelems = parts[k].elemstosend[*iter];
          std::set<idx_t>::iterator it;
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(parts[k].get_neighbors(*it,j) != -1){
                parts[k].elemstosend[*iter].insert(parts[k].get_neighbors(*it,j));
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  parts[k].nodestosend[*iter].insert(parts[k].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(std::set<idx_t>::iterator it=parts[k].boundaryelems.begin();
                it!=parts[k].boundaryelems.end();
                it++){
              for(idx_t j=0;j<esize;j++){
                if(parts[k].nodestosend[*iter].find(parts[k].get_elems(*it,j))!=parts[k].nodestosend[*iter].end()){
                  parts[k].elemstosend[*iter].insert(*it);
                }
              }
            }
            for(std::set<idx_t>::iterator it=parts[k].elemstosend[*iter].begin();
                it!=parts[k].elemstosend[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                parts[k].nodestosend[*iter].insert(parts[k].get_elems(*it,p));
              }
            }
          }

          //Add to recv list
          tmpnewboundaryelems.clear();
          tmpnewboundaryelems = parts[k].elemstorecv[*iter];
          for(it=tmpnewboundaryelems.begin();it!=tmpnewboundaryelems.end();it++){
            for(idx_t j=0;j<esize;j++){
              if(g_potentialneighbors[*iter].neighbors[*it*esize+j] != -1){
                parts[k].elemstorecv[*iter].insert(g_potentialneighbors[*iter].neighbors[*it*esize+j]);
                //Add boundarynodes
                for(idx_t j=0; j<esize; j++){
                  parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,j));
                }
              }
            } 
          }
          //If method=1 add cells that have only a point in the boundary points
          if(method==1){
            for(idx_t i=0; i<g_potentialneighbors[*iter].get_nelems(); i++){
              for(idx_t j=0;j<esize;j++){
                if(parts[k].nodestorecv[*iter].find(g_potentialneighbors[*iter].get_elems(i,j))!=parts[k].nodestorecv[*iter].end()){
                  parts[k].elemstorecv[*iter].insert(i);
                }
              }
            }
            for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
                it!=parts[k].elemstorecv[*iter].end();
                it++){
              for(idx_t p=0; p<esize; p++){
                parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
              }
            }
          }


        }

        for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
            it!=parts[k].elemstorecv[*iter].end();it++){
          for(idx_t p=0; p<esize; p++){
            parts[k].nodestorecv[*iter].insert(g_potentialneighbors[*iter].get_elems(*it,p));
          }
        }
      }

    }
  }

}

void share_boundaryFixPotentialNeighbors(std::vector<partition> &parts, std::vector<idx_t> &ownerofpartition, MPI_Comm comm){

  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;
  MPI_Request *requestSendPtr, *requestRecvPtr;

  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  std::vector<std::set<idx_t> > msgs_send(nprocs);
  std::vector<std::set<idx_t> > msgs_recv(nprocs);

  idx_t esize=parts[0].esize;


  std::set<idx_t>::iterator it;
  idx_t iloc(0);
  idx_t idum;

  //Construct partitionlocid (reverse lookup of mesh id)
  std::map<idx_t, idx_t> partitionlocid;
  for(idx_t i=0;i<parts.size();i++){
    partitionlocid.insert(std::make_pair(parts[i].partitionid, i));
  }

  //Find max number of nodes and elems to send
  idx_t maxnodestosend(0), maxelemstosend(0);
  for(idx_t k=0; k<parts.size(); k++){
    maxnodestosend = std::max(maxnodestosend, (idx_t) parts[k].l1_boundarynodes.size());
  }

  idx_t *nodestosend = new idx_t[maxnodestosend*parts.size()];

  //Compute total number of messages/request
  idx_t ntotreq_send=0;
  idx_t ntotreq_recv=0;
  idx_t nreq=0;
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors_extents.begin(); it!= parts[k].potentialneighbors_extents.end(); it++){
      if(ownerofpartition[*it]!=me && (msgs_send[ownerofpartition[*it]].count(parts[k].partitionid)==0)){
        msgs_send[ownerofpartition[*it]].insert(parts[k].partitionid);
        ntotreq_send++;
      }
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        msgs_recv[ownerofpartition[*it]].insert(*it);
        ntotreq_recv++;
      }
    }
  }

  requestSendPtr = new MPI_Request[ntotreq_send];

  //Construc contiguous vector of nodes  to send
  for(idx_t k=0; k<parts.size(); k++){
    iloc=0;
    for(it=parts[k].l1_boundarynodes.begin(); it!=parts[k].l1_boundarynodes.end(); it++){
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors_extents.begin(); it!= parts[k].potentialneighbors_extents.end(); it++){
      if(msgs_send[ownerofpartition[*it]].count(parts[k].partitionid)==0){
        if(ownerofpartition[*it]!=me){
          MPI_Isend(&nodestosend[maxnodestosend*k], parts[k].l1_boundarynodes.size(), IDX_T, ownerofpartition[*it], (*it)*ownerofpartition.size()+parts[k].partitionid, comm, &requestSendPtr[nreq]);
          nreq++;
          msgs_send[ownerofpartition[*it]].insert(parts[k].partitionid);
        }
        else{
          //Check for node in common, if found add to potentialnieghbors
          idx_t k2 = partitionlocid[*it];
          for(idx_t i=0; i<parts[k].l1_boundarynodes.size(); i++){
            if(parts[k2].l1_boundarynodes.count(nodestosend[maxnodestosend*k+i]) != 0){
              parts[k2].potentialneighbors.insert(parts[k].partitionid);
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors_extents.begin(); it!= parts[k].potentialneighbors_extents.end(); it++){
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        MPI_Probe(ownerofpartition[*it], (*it)+ownerofpartition.size()*parts[k].partitionid, comm, &stat);
        MPI_Get_count(&stat,IDX_T,&nnodestorecv[nreq]);
        maxnodestorecv = std::max(maxnodestorecv,nnodestorecv[nreq]);
        msgs_recv[ownerofpartition[*it]].insert(*it);
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors_extents.begin(); it!= parts[k].potentialneighbors_extents.end(); it++){
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        MPI_Recv(&nodestorecv[nreq*maxnodestorecv], nnodestorecv[nreq], IDX_T, ownerofpartition[*it], (*it)+ownerofpartition.size()*parts[k].partitionid, comm, &stat);

        //Check for node in common, if found add to potentialnieghbors
        for(idx_t i=0; i<nnodestorecv[nreq]; i++){
          if(parts[k].l1_boundarynodes.count(nodestorecv[nreq*maxnodestorecv+i]) != 0){
            parts[k].potentialneighbors.insert(*it);
            break;
          }
        }
        msgs_recv[ownerofpartition[*it]].insert(*it);
        msgs_req_recv[ownerofpartition[*it]].insert({*it,nreq});
        nreq++;
      }
      else if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==1){
        //Find what req the msg came with
        idx_t loc_req=msgs_req_recv[ownerofpartition[*it]][*it];
        //Check for node in common, if found add to potentialnieghbors
        for(idx_t i=0; i<nnodestorecv[loc_req]; i++){
          if(parts[k].l1_boundarynodes.count(nodestorecv[loc_req*maxnodestorecv+i]) != 0){
            parts[k].potentialneighbors.insert(*it);
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

void share_boundary(std::vector<partition> &parts, std::vector<idx_t> &ownerofpartition, MPI_Comm comm){

  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;
  MPI_Request *requestSendPtr, *requestRecvPtr;

  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  std::vector<std::set<idx_t> > msgs_send(nprocs);
  std::vector<std::set<idx_t> > msgs_recv(nprocs);

  idx_t esize=parts[0].esize;


  std::set<idx_t>::iterator it;
  idx_t iloc(0);
  idx_t idum;

  //Construct partitionlocid (reverse lookup of mesh id)
  std::map<idx_t, idx_t> partitionlocid;
  for(idx_t i=0;i<parts.size();i++){
    partitionlocid.insert(std::make_pair(parts[i].partitionid, i));
  }

  //Find max number of nodes and elems to send
  idx_t maxnodestosend(0), maxelemstosend(0);
  for(idx_t k=0; k<parts.size(); k++){
    maxnodestosend = std::max(maxnodestosend, (idx_t) parts[k].boundarynodes.size());
    maxelemstosend = std::max(maxelemstosend, (idx_t) parts[k].boundaryelems.size());
  }

  real_t *nodestosend = new real_t[maxnodestosend*4*parts.size()];
  idx_t *elemstosend = new idx_t[maxelemstosend*(1+esize)*parts.size()];

  //Compute total number of messages/request
  idx_t ntotreq_send=0;
  idx_t ntotreq_recv=0;
  idx_t nreq=0;
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(ownerofpartition[*it]!=me && (msgs_send[ownerofpartition[*it]].count(parts[k].partitionid)==0)){
        msgs_send[ownerofpartition[*it]].insert(parts[k].partitionid);
        ntotreq_send++;
      }
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        msgs_recv[ownerofpartition[*it]].insert(*it);
        ntotreq_recv++;
      }
    }
  }

  requestSendPtr = new MPI_Request[ntotreq_send];

  //Construc contiguous vector of nodes  to send
  for(idx_t k=0; k<parts.size(); k++){
    iloc=0;
    for(it=parts[k].boundarynodes.begin(); it!=parts[k].boundarynodes.end(); it++){
      idum = parts[k].nodes_gtl[*it];
      nodestosend[4*k*maxnodestosend+4*iloc] = (idx_t) *it;
      nodestosend[4*k*maxnodestosend+4*iloc+1] = parts[k].get_nodes(idum,0);
      nodestosend[4*k*maxnodestosend+4*iloc+2] = parts[k].get_nodes(idum,1);
      nodestosend[4*k*maxnodestosend+4*iloc+3] = parts[k].get_nodes(idum,2);
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(msgs_send[ownerofpartition[*it]].count(parts[k].partitionid)==0){
        if(ownerofpartition[*it]!=me){
          MPI_Isend(&nodestosend[4*maxnodestosend*k], parts[k].boundarynodes.size()*4, REAL_T, ownerofpartition[*it], (*it)*ownerofpartition.size()+parts[k].partitionid, comm, &requestSendPtr[nreq]);
          nreq++;
        }
        else{
          potential_neighbors_boundary pnbound;
          pnbound.nodes.resize(parts[k].boundarynodes.size()*3,0.);
          for(idx_t i=0;i<parts[k].boundarynodes.size();i++){
            pnbound.nodes_ltg.insert(std::make_pair(i,nodestosend[4*maxnodestosend*k+i*4]));
            pnbound.nodes_gtl.insert(std::make_pair(nodestosend[4*maxnodestosend*k+i*4],i));
            pnbound.get_nodes(i,0) = nodestosend[4*maxnodestosend*k+i*4+1];
            pnbound.get_nodes(i,1) = nodestosend[4*maxnodestosend*k+i*4+2];
            pnbound.get_nodes(i,2) = nodestosend[4*maxnodestosend*k+i*4+3];
          }
          g_potentialneighbors.insert(std::make_pair(parts[k].partitionid,pnbound));
          g_potentialneighbors[parts[k].partitionid].esize = parts[k].esize;
        }
        msgs_send[ownerofpartition[*it]].insert(parts[k].partitionid);
      }
    }
  }

  //Probe to get sizes
  idx_t *nnodestorecv = new idx_t[ntotreq_recv];
  idx_t maxnodestorecv=0;

  nreq=0;
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        MPI_Probe(ownerofpartition[*it], (*it)+ownerofpartition.size()*parts[k].partitionid, comm, &stat);
        MPI_Get_count(&stat,REAL_T,&nnodestorecv[nreq]);
        nnodestorecv[nreq] /= 4;
        maxnodestorecv = std::max(maxnodestorecv,nnodestorecv[nreq]);
        msgs_recv[ownerofpartition[*it]].insert(*it);
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        MPI_Recv(&nodestorecv[nreq*maxnodestorecv*4], nnodestorecv[nreq]*4, REAL_T, ownerofpartition[*it], (*it)+ownerofpartition.size()*parts[k].partitionid, comm, &stat);

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
        g_potentialneighbors[*it].esize = parts[k].esize;
        msgs_recv[ownerofpartition[*it]].insert(*it);
        nreq++;
      }
    }
  }

  MPI_Waitall(ntotreq_send, requestSendPtr, statSend);

  //Construc contiguous vector of elems to send
  for(idx_t k=0; k<parts.size(); k++){
    iloc=0;
    for(it=parts[k].boundaryelems.begin(); it!=parts[k].boundaryelems.end(); it++){
      idum = parts[k].elems_ltg[*it];
      elemstosend[(1+esize)*(maxelemstosend*k+iloc)] = idum;
      for(idx_t j=0;j<esize;j++){
        elemstosend[(1+esize)*(maxelemstosend*k+iloc)+j+1] = parts[k].get_elems(*it,j);
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(msgs_send[ownerofpartition[*it]].count(parts[k].partitionid)==0){
        if(ownerofpartition[*it]!=me){
          MPI_Isend(&elemstosend[(1+esize)*maxelemstosend*k], parts[k].boundaryelems.size()*(1+esize), IDX_T, ownerofpartition[*it], (*it)*ownerofpartition.size()+parts[k].partitionid, comm, &requestSendPtr[nreq]);
          nreq++;
        }
        else{
          g_potentialneighbors[parts[k].partitionid].esize = esize;
          g_potentialneighbors[parts[k].partitionid].elems.resize(parts[k].boundaryelems.size()*esize,0);
          for(idx_t i=0;i<parts[k].boundaryelems.size();i++){
            g_potentialneighbors[parts[k].partitionid].elems_ltg.insert(std::make_pair(i,elemstosend[(1+esize)*(maxelemstosend*k+i)]));
            g_potentialneighbors[parts[k].partitionid].elems_gtl.insert(std::make_pair(elemstosend[(1+esize)*(maxelemstosend*k+i)],i));
            for(idx_t j=0;j<esize;j++){
              g_potentialneighbors[parts[k].partitionid].get_elems(i,j) = elemstosend[(1+esize)*(maxelemstosend*k+i)+j+1];
            }
          }
        }
        msgs_send[ownerofpartition[*it]].insert(parts[k].partitionid);
      }
    }
  }

  //Probe to get sizes
  idx_t *nelemstorecv = new idx_t[ntotreq_recv];
  idx_t maxelemstorecv=0;

  nreq=0;
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        MPI_Probe(ownerofpartition[*it], (*it)+ownerofpartition.size()*parts[k].partitionid, comm, &stat);
        MPI_Get_count(&stat,IDX_T,&nelemstorecv[nreq]);
        nelemstorecv[nreq] /= (1+esize);
        maxelemstorecv = std::max(maxelemstorecv,nelemstorecv[nreq]);
        msgs_recv[ownerofpartition[*it]].insert(*it);
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
  for(idx_t k=0; k<parts.size(); k++){
    for(it=parts[k].potentialneighbors.begin(); it!= parts[k].potentialneighbors.end(); it++){
      if(ownerofpartition[*it]!=me && msgs_recv[ownerofpartition[*it]].count(*it)==0){
        MPI_Recv(&elemstorecv[nreq*(1+esize)*maxelemstorecv], nelemstorecv[nreq]*(1+esize), IDX_T, ownerofpartition[*it], (*it)+ownerofpartition.size()*parts[k].partitionid, comm, &stat);

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
        msgs_recv[ownerofpartition[*it]].insert(*it);
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

void compute_potential_neighbors_from_extents(idx_t ntotparts, std::vector<partition> &parts, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Status stat;
  MPI_Status *statSend, *statRecv;

  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&me);

  /* std::vector<std::vector<std::vector<real_t> > > allextents(ntotparts,std::vector<std::vector<real_t> >(3, std::vector<real_t>(2, 0.))); */
  real_t *allextents = new real_t[ntotparts*3*2];
  for(idx_t k=0;k<ntotparts*2*3;k++){
    allextents[k] = 0.;
  }

  //Compute each partition extents
  for(idx_t k=0; k<parts.size(); k++){
    parts[k].compute_extents();
    allextents[parts[k].partitionid*3*2] = parts[k].extents[0][0];
    allextents[parts[k].partitionid*3*2 + 1] = parts[k].extents[0][1];
    allextents[parts[k].partitionid*3*2 + 2] = parts[k].extents[1][0];
    allextents[parts[k].partitionid*3*2 + 3] = parts[k].extents[1][1];
    allextents[parts[k].partitionid*3*2 + 4] = parts[k].extents[2][0];
    allextents[parts[k].partitionid*3*2 + 5] = parts[k].extents[2][1];
  }

  //Share extents with everyone
  MPI_Allreduce(MPI_IN_PLACE, allextents, ntotparts*3*2, REAL_T, MPI_SUM, comm);

  real_t eps=1e-10;

  bool interx, intery, interz;
  //Compute potential neighbors
  for(idx_t k=0; k<parts.size(); k++){
    for(idx_t p=0; p<ntotparts; p++){
      if(parts[k].partitionid!=p){
        interx = ((allextents[parts[k].partitionid*6+1]>=allextents[p*6]-eps) and
            (allextents[parts[k].partitionid*6+1]<=allextents[p*6+1]+eps));
        interx = interx or ((allextents[parts[k].partitionid*6]>=allextents[p*6]-eps) and
            (allextents[parts[k].partitionid*6]<=allextents[p*6+1]+eps));
        interx = interx or ((allextents[parts[k].partitionid*6]<=allextents[p*6]+eps) and
            (allextents[parts[k].partitionid*6+1]>=allextents[p*6+1]-eps));

        intery = ((allextents[parts[k].partitionid*6+3]>=allextents[p*6+2]-eps) and
            (allextents[parts[k].partitionid*6+3]<=allextents[p*6+3]+eps));
        intery = intery or ((allextents[parts[k].partitionid*6+2]>=allextents[p*6+2]-eps) and (allextents[parts[k].partitionid*6+2]<=allextents[p*6+3]+eps));
        intery = intery or ((allextents[parts[k].partitionid*6+2]<=allextents[p*6+2]+eps) and (allextents[parts[k].partitionid*6+3]>=allextents[p*6+3]-eps));

        interz = ((allextents[parts[k].partitionid*6+5]>=allextents[p*6+4]-eps) and
            (allextents[parts[k].partitionid*6+5]<=allextents[p*6+5]+eps));
        interz = interz or ((allextents[parts[k].partitionid*6+4]>=allextents[p*6+4]-eps) and (allextents[parts[k].partitionid*6+4]<=allextents[p*6+5]+eps));
        interz = interz or ((allextents[parts[k].partitionid*6+4]<=allextents[p*6+4]+eps) and (allextents[parts[k].partitionid*6+5]>=allextents[p*6+5]-eps));

        if((interx and intery) and interz){
          parts[k].potentialneighbors_extents.insert(p);
        }
      }
    }
  }
}

void add_elemsAndRenumber(std::vector<partition> &parts){
  //Renumber send elems
  for(idx_t k=0; k<parts.size(); k++){
    //Add Recv elements to local epartition
    for(std::set<idx_t>::iterator iter=parts[k].potentialneighbors.begin();
        iter!=parts[k].potentialneighbors.end();
        iter++){
      for(std::set<idx_t>::iterator it=parts[k].nodestorecv[*iter].begin();
          it!=parts[k].nodestorecv[*iter].end();
          it++){
        if(parts[k].nodes_gtl.find(*it)==parts[k].nodes_gtl.end()){
          parts[k].nodes_ltg.insert(std::make_pair(parts[k].get_nnodes(), *it));
          parts[k].nodes_gtl.insert(std::make_pair(*it, parts[k].get_nnodes()));
          for(int kk=0; kk<3; kk++){
            parts[k].nodes.push_back(g_potentialneighbors[*iter].get_nodes(g_potentialneighbors[*iter].nodes_gtl[*it],kk));
          }
        }
      }
      for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin();
          it!=parts[k].elemstorecv[*iter].end();
          it++){
        idx_t itglob = g_potentialneighbors[*iter].elems_ltg[*it];
        parts[k].elems_ltg.insert(std::make_pair(parts[k].get_nelems(), itglob));
        parts[k].elems_gtl.insert(std::make_pair(itglob, parts[k].get_nelems()));
        for(int kk=0; kk<parts[k].esize; kk++){
          parts[k].elems.push_back(g_potentialneighbors[*iter].get_elems(*it,kk));
        }
      }
    }

    //Renumbering
    parts[k].renumber_otn.resize(parts[k].get_nelems(), -1);
    parts[k].renumber_nto.resize(parts[k].get_nelems(), -1);
    idx_t itloc;
    idx_t ind=0;
    bool norsendorrecv;
    for(idx_t i=0; i<parts[k].get_nelems(); i++){
      norsendorrecv = true;
      for(std::set<idx_t>::iterator iter=parts[k].potentialneighbors.begin();
          iter!=parts[k].potentialneighbors.end();
          iter++){

        if(parts[k].elemstosend[*iter].find(i)!=parts[k].elemstosend[*iter].end()){
          norsendorrecv = false;
        }

        if(g_potentialneighbors[*iter].elems_gtl.find(parts[k].elems_ltg[i])!=g_potentialneighbors[*iter].elems_gtl.end()){
          itloc = g_potentialneighbors[*iter].elems_gtl[parts[k].elems_ltg[i]];
          if(parts[k].elemstorecv[*iter].find(itloc)!=parts[k].elemstorecv[*iter].end()){
            norsendorrecv = false;
          }
        }
      }

      if(norsendorrecv){
        parts[k].renumber_otn[i] = ind;
        parts[k].renumber_nto[ind] = i;
        ind++;
      }
    }

    for(std::set<idx_t>::iterator iter=parts[k].potentialneighbors.begin();
        iter!=parts[k].potentialneighbors.end();
        iter++){

      std::vector<idx_t> sort_elemstosend;
      sort_elemstosend.assign(parts[k].elemstosend[*iter].begin(), parts[k].elemstosend[*iter].end());
      for(idx_t ie=0; ie<sort_elemstosend.size(); ie++){
        for(idx_t je=0; je<sort_elemstosend.size()-ie-1; je++){
          /* if(sort_elemstosend[je] > sort_elemstosend[je+1]){ */
          if(parts[k].elems_ltg[sort_elemstosend[je]] > parts[k].elems_ltg[sort_elemstosend[je+1]]){
            idx_t dummy = sort_elemstosend[je];
            sort_elemstosend[je] = sort_elemstosend[je+1];
            sort_elemstosend[je+1] = dummy;
          }
        }
        }

        for(idx_t ie=0; ie<sort_elemstosend.size(); ie++){
          if(parts[k].renumber_otn[sort_elemstosend[ie]] == -1){
            parts[k].renumber_otn[sort_elemstosend[ie]] = ind;
            parts[k].renumber_nto[ind] = sort_elemstosend[ie];
            ind++;
          }
        }


        /* for(std::set<idx_t>::iterator it=parts[k].elemstosend[*iter].begin(); */
        /*     it!=parts[k].elemstosend[*iter].end(); */
        /*     it++){ */
        /* if(parts[k].renumber_otn[*it] == -1){ */
        /*   parts[k].renumber_otn[*it] = ind; */
        /*   parts[k].renumber_nto[ind] = *it; */
        /*   ind++; */
        /* } */
        /* } */
      }

      for(std::set<idx_t>::iterator iter=parts[k].potentialneighbors.begin();
          iter!=parts[k].potentialneighbors.end();
          iter++){

        std::vector<idx_t> sort_elemstorecv;
        sort_elemstorecv.assign(parts[k].elemstorecv[*iter].begin(), parts[k].elemstorecv[*iter].end());
        for(idx_t ie=0; ie<sort_elemstorecv.size(); ie++){
          for(idx_t je=0; je<sort_elemstorecv.size()-ie-1; je++){
            /* if(sort_elemstorecv[je] > sort_elemstorecv[je+1]){ */
            if(g_potentialneighbors[*iter].elems_ltg[sort_elemstorecv[je]] > g_potentialneighbors[*iter].elems_ltg[sort_elemstorecv[je+1]]){
              idx_t dummy = sort_elemstorecv[je];
              sort_elemstorecv[je] = sort_elemstorecv[je+1];
              sort_elemstorecv[je+1] = dummy;
            }
          }
          }

          for(idx_t ie=0; ie < sort_elemstorecv.size(); ie++){
            idx_t itglobloc = parts[k].elems_gtl[g_potentialneighbors[*iter].elems_ltg[sort_elemstorecv[ie]]];
            parts[k].renumber_otn[itglobloc] = ind;
            parts[k].renumber_nto[ind] = itglobloc;
            ind++;
          }

          /* for(std::set<idx_t>::iterator it=parts[k].elemstorecv[*iter].begin(); */
          /*     it!=parts[k].elemstorecv[*iter].end(); */
          /*     it++){ */
          /*   idx_t itglobloc = parts[k].elems_gtl[g_potentialneighbors[*iter].elems_ltg[*it]]; */
          /*   parts[k].renumber_otn[itglobloc] = ind; */
          /*   parts[k].renumber_nto[ind] = itglobloc; */
          /*   ind++; */
          /* } */
        }
      }
    }
