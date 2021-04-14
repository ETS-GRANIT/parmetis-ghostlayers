#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include <algorithm>
#include <assert.h>
#include <chrono>
#include <thread>

#include "read_write.hpp"
#include "module_parmetis.hpp"

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
  zone = 1; //Reading zone 1 only
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
  cgsize_t TeBeg, TeEnd, *conn;
  int parentFlag;
  if(cg_section_read(index_file, base, zone, sec, secname, &type,
        &TeBeg, &TeEnd, &nBdry, &parentFlag) != CG_OK) cg_get_error();
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
  cgsize_t eBeg=elmdist[me]+TeBeg;
  cgsize_t eEnd=elmdist[me+1]+TeBeg-1;
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

void writeneighbors(std::vector<submesh> &submeshesowned, idx_t esize){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_neighbors.dat";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    for(idx_t i=0; i<submeshesowned[k].get_nelems(); i++){
      outFile << submeshesowned[k].elems_ltg[i] << " ";
      for(idx_t j=0; j<submeshesowned[k].esize; j++){
        idx_t value = submeshesowned[k].get_neighbors(i,j);
        if(value!=-1){
          outFile << submeshesowned[k].elems_ltg[submeshesowned[k].get_neighbors(i,j)] << " ";
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
    idx_t nNodes = submeshesowned[k].get_nnodes();
    idx_t nElems = submeshesowned[k].get_nelems();
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << submeshesowned[k].submeshid << "_output_cuteflow.txt";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "Table de coordonnees" << std::endl;
    outFile << nNodes << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << i+1 << " " << submeshesowned[k].get_nodes(i,0) << " " << submeshesowned[k].get_nodes(i,1) << " " << submeshesowned[k].get_nodes(i,2) << std::endl;
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems << " " << nElems*(esize+1) << std::endl;
    outFile << "Table de connectivité" << std::endl;
    outFile << nElems << std::endl;
    for(idx_t i=0; i<nElems; i++){
      outFile << i+1 << " ";
      indrenum = submeshesowned[k].renumber_otn[i];
      for(idx_t p=0; p<esize; p++){
        outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(indrenum,p)] << " ";
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
    fname << submeshesowned[k].submeshid << "_Mesh_Output.cgns";
    if(cg_open(fname.str().c_str(),CG_MODE_WRITE,&index_file)) cg_error_exit();
    if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base)) cg_error_exit();
    /* cgsize_t estart, eend; */
    cgsize_t isize[1][3];
    isize[0][0]=submeshesowned[k].get_nnodes();
    isize[0][1]=submeshesowned[k].get_nelems();
    isize[0][2]=0;
    std::stringstream zss;
    zss << "River_" << submeshesowned[k].submeshid;
    if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&index_zone) ) cg_error_exit();
    coord = new double[submeshesowned[k].get_nnodes()];
    for(idx_t i=0; i<submeshesowned[k].get_nnodes(); i++){
      coord[i] = submeshesowned[k].get_nodes(i,0);
    }
    if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateX", coord, &Cx)) cg_error_exit();
    for(idx_t i=0; i<submeshesowned[k].get_nnodes(); i++){
      coord[i] = submeshesowned[k].get_nodes(i,1);
    }
    if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateY", coord, &Cy)) cg_error_exit();
    for(idx_t i=0; i<submeshesowned[k].get_nnodes(); i++){
      coord[i] = submeshesowned[k].get_nodes(i,2);
    }
    if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateZ", coord, &Cz)) cg_error_exit();
    delete [] coord;
    int nstart=1, nend=submeshesowned[k].get_nelems();
    elems = new cgsize_t[esize*submeshesowned[k].get_nelems()];
    for(idx_t i=0; i<submeshesowned[k].get_nelems(); i++){
      for(idx_t j=0; j<esize; j++){
        elems[esize*i + j] = submeshesowned[k].renumber_otn[submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(i,j)]]+1;
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
          elemstosend[i] = submeshesowned[k].renumber_otn[*iter]+1;
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
        elemstorecv[0] = submeshesowned[k].renumber_otn[itglobloc]+1;
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
        if(dimensions==submeshesowned[k].get_nnodes()){
          location = CGNS_ENUMV(Vertex);
        }
        if(dimensions==submeshesowned[k].get_nelems()){
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
        isize[0][0]=submeshesowned[k].get_nnodes();
        isize[0][1]=submeshesowned[k].get_nelems();
        isize[0][2]=0;
        std::stringstream zss;
        zss << "River_" << submeshesowned[k].submeshid;
        if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&index_zone) ) cg_error_exit();
        coord = new double[submeshesowned[k].get_nnodes()];
        for(idx_t i=0; i<submeshesowned[k].get_nnodes(); i++){
          coord[i] = submeshesowned[k].get_nodes(i,0);
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateX", coord, &Cx)) cg_error_exit();
        for(idx_t i=0; i<submeshesowned[k].get_nnodes(); i++){
          coord[i] = submeshesowned[k].get_nodes(i,1);
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateY", coord, &Cy)) cg_error_exit();
        for(idx_t i=0; i<submeshesowned[k].get_nnodes(); i++){
          coord[i] = submeshesowned[k].get_nodes(i,2);
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateZ", coord, &Cz)) cg_error_exit();
        delete [] coord;
        int nstart=1, nend=submeshesowned[k].get_nelems();
        elems = new cgsize_t[esize*submeshesowned[k].get_nelems()];
        for(idx_t i=0; i<submeshesowned[k].get_nelems(); i++){
          for(idx_t j=0; j<esize; j++){
            elems[esize*i + j] = submeshesowned[k].renumber_otn[submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(i,j)]]+1;
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
              elemstosend[i] = submeshesowned[k].renumber_otn[*iter]+1;
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
            elemstorecv[0] = submeshesowned[k].renumber_otn[itglobloc]+1;
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
            if(dimensions==submeshesowned[k].get_nnodes()){
              location = CGNS_ENUMV(Vertex);
            }
            if(dimensions==submeshesowned[k].get_nelems()){
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
                  submeshesowned[k].arrays[nud-1][na-1].resize(submeshesowned[k].get_nelems());
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

void writeworecvVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t indrenum;
    idx_t nNodes = submeshesowned[k].get_nnodes();
    idx_t nElems = submeshesowned[k].get_nelems();
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
      outFile << submeshesowned[k].get_nodes(i,0) << " " << submeshesowned[k].get_nodes(i,1) << " " << submeshesowned[k].get_nodes(i,2) << std::endl;
    }

    idx_t totelemsrecv=0;
    for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[k].elemstorecv.begin(); it!= submeshesowned[k].elemstorecv.end(); it++){
      totelemsrecv += it->second.size();     
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems-totelemsrecv << " " << (nElems-totelemsrecv)*(esize+1) << std::endl;
    for(idx_t i=0; i<nElems; i++){
      /* indrenum = i; */
      indrenum = submeshesowned[k].renumber_nto[i];
      if(indrenum<nElems-totelemsrecv){
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(indrenum,p)] << " ";
        }
        outFile << std::endl;
      }
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

void writeVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t indrenum;
    idx_t nNodes = submeshesowned[k].get_nnodes();
    idx_t nElems = submeshesowned[k].get_nelems();
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
      outFile << submeshesowned[k].get_nodes(i,0) << " " << submeshesowned[k].get_nodes(i,1) << " " << submeshesowned[k].get_nodes(i,2) << std::endl;
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems << " " << nElems*(esize+1) << std::endl;
    for(idx_t i=0; i<nElems; i++){
      outFile << esize << " ";
      indrenum = submeshesowned[k].renumber_nto[i];
      /* indrenum = i; */
      for(idx_t p=0; p<esize; p++){
        outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(indrenum,p)] << " ";
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
    idx_t nNodes = submeshesowned[k].get_nnodes();

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
        outFile << submeshesowned[k].renumber_otn[*it1] << std::endl;
      }
    }

    for(std::map<idx_t,idx_t>::iterator it=trueneighborsrecv.begin();
        it!=trueneighborsrecv.end();it++){
      for(std::set<idx_t>::iterator it1=submeshesowned[k].elemstorecv[it->first].begin();
          it1!=submeshesowned[k].elemstorecv[it->first].end();it1++){
        idx_t itglobloc = submeshesowned[k].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*it1]];
        outFile << it->first << " " << it->second << " " << submeshesowned[k].renumber_otn[itglobloc] << std::endl;
        break;
      }
    }


  }

}
void writerecvVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){

  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].get_nnodes();
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
      outFile << submeshesowned[k].get_nodes(i,0) << " " << submeshesowned[k].get_nodes(i,1) << " " << submeshesowned[k].get_nodes(i,2) << std::endl;
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
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(itloc,p)] << " ";
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
      outFile << submeshesowned[k].get_nodes(i,2) << std::endl;
    }
  }
}

void writesendVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].get_nnodes();
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
      outFile << submeshesowned[k].get_nodes(i,0) << " " << submeshesowned[k].get_nodes(i,1) << " " << submeshesowned[k].get_nodes(i,2) << std::endl;
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
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(*it,p)] << " ";
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
      outFile << submeshesowned[k].get_nodes(i,2) << std::endl;
    }
  }
}
void writeboundaryVTK(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim){
  for(idx_t k=0; k<submeshesowned.size(); k++){
    idx_t nNodes = submeshesowned[k].get_nnodes();
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
      outFile << submeshesowned[k].get_nodes(i,0) << " " << submeshesowned[k].get_nodes(i,1) << " " << submeshesowned[k].get_nodes(i,2) << std::endl;
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
          outFile << submeshesowned[k].nodes_gtl[submeshesowned[k].get_elems(*it,p)] << " ";
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
      outFile << submeshesowned[k].get_nodes(i,2) << std::endl;
    }
  }
}

void updateNodes(std::vector<submesh> &submeshesowned, std::string nodesFilename, MPI_Comm comm){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  for(idx_t k=0; k<submeshesowned.size(); k++){
    /* submeshesowned[k].nnodes = submeshesowned[k].nodesindex.size(); */
    submeshesowned[k].nodes.resize(3*submeshesowned[k].nodesindex.size(), 0.);
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
        submeshesowned[k].get_nodes(nloc[k],0) = sx;
        submeshesowned[k].get_nodes(nloc[k],1) = sy;
        submeshesowned[k].get_nodes(nloc[k],2) = sz;
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

  /* std::cout << me << " IUPN " << gnnodes << " " << gnelems << std::endl; */

  if(cg_ncoords(index_file, base, zone, &nCoords) != CG_OK) cg_get_error();
  double *x, *y, *z;


  //Resize each subdomain nodes vector
  for(idx_t k=0; k<submeshesowned.size(); k++){
    /* submeshesowned[k].nnodes = submeshesowned[k].nodesindex.size(); */
    submeshesowned[k].nodes.resize(3*submeshesowned[k].nodesindex.size(), 0.);
  }

  //Construct a global_nodes vector so that we only check once if we own a node
  std::map<idx_t, std::set<idx_t>> global_nodes;
  for(idx_t k=0; k<submeshesowned.size(); k++){
    for(std::unordered_set<idx_t>::iterator it=submeshesowned[k].nodesindex.begin();
        it!=submeshesowned[k].nodesindex.end(); it++){
      global_nodes[*it].insert(k);
    }
  }

  //Read nodes in batches
  std::vector<idx_t> nloc(submeshesowned.size(), 0);
  //int nbatches=1;
  int nbatches=nprocs;
  /* cgsize_t maxnnodesbatch = gnnodes - ((nbatches-1)*gnnodes/nbatches + 1) + 1; */
  cgsize_t nnodesperbatch = gnnodes/nbatches;
  cgsize_t maxnnodesbatch = nnodesperbatch + nbatches + 1;
  x = new double[maxnnodesbatch];
  y = new double[maxnnodesbatch];
  z = new double[maxnnodesbatch];
  for(int b=0; b<nbatches; b++){
    cgsize_t firstnode = b*nnodesperbatch + 1;
    cgsize_t lastnode = firstnode + nnodesperbatch-1;
    if(b==(nbatches-1)){//last batch
      lastnode = gnnodes;
    }
    cgsize_t nnodesbatch = (lastnode-firstnode+1);

    if(cg_coord_info(index_file, base, zone, 1, &dataType, zonename) != CG_OK) cg_get_error();
    if(cg_coord_read(index_file, base, zone, zonename, CGNS_ENUMV(RealDouble),
          &firstnode, &lastnode, x) != CG_OK) cg_get_error();

    if(cg_coord_info(index_file, base, zone, 2, &dataType, zonename) != CG_OK) cg_get_error();
    if(cg_coord_read(index_file, base, zone, zonename, CGNS_ENUMV(RealDouble),
          &firstnode, &lastnode, y) != CG_OK) cg_get_error();

    if(cg_coord_info(index_file, base, zone, 3, &dataType, zonename) != CG_OK) cg_get_error();
    if(cg_coord_read(index_file, base, zone, zonename, CGNS_ENUMV(RealDouble),
          &firstnode, &lastnode, z) != CG_OK) cg_get_error();

    real_t sx, sy, sz;
    int j=0;
    for(idx_t i=firstnode-1; i<lastnode; i++){
      if(global_nodes.count(i)!=0){
        for(std::set<idx_t>::iterator it=global_nodes[i].begin();
            it!=global_nodes[i].end(); it++){
          int k=*it;
          submeshesowned[k].get_nodes(nloc[k],0) = x[j];
          submeshesowned[k].get_nodes(nloc[k],1) = y[j];
          submeshesowned[k].get_nodes(nloc[k],2) = z[j];
          submeshesowned[k].nodes_ltg.insert({nloc[k], i});
          submeshesowned[k].nodes_gtl.insert({i, nloc[k]});
          nloc[k] += 1;
        }
      }
      j++;
    }
  }
  if (cg_close(index_file)) cg_error_exit();
}

void writeMeshPCGNS_wos(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm(MPI_COMM_WORLD);

  int index_file, index_base, index_zone, n_bases, base, physDim, cellDim, nZones, zone, index_coord;
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
  cgsize_t *elems, *elemstosend, *elemstorecv;
  double *coord;
  int Cx, Cy, Cz;

  /* if(cgp_mpi_comm(MPI_COMM_WORLD)) cgp_error_exit(); */
  if(cgp_pio_mode(CGP_INDEPENDENT)) cgp_error_exit();
  if(cgp_open("Mesh_Output_pcgns_wos.cgns",CG_MODE_WRITE,&index_file)) cgp_error_exit();


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

  if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base) || cg_goto(index_file, index_base,"end")) cg_error_exit();


  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  std::vector<std::string> zms(ownerofsubmesh.size());
  std::vector<std::string> ems(ownerofsubmesh.size());
  for(int k=0; k<ownerofsubmesh.size(); k++){
    std::stringstream zss;
    zss << "River_" << k;
    zms[k] = zss.str();
  }

  int ns[ownerofsubmesh.size()];
  int ne[ownerofsubmesh.size()];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k] == me){
      ns[k] = submeshesowned[submeshlocid[k]].get_nnodes();
      ne[k] = submeshesowned[submeshlocid[k]].get_nelems();
    }
    else{
      ns[k] = 0;
      ne[k] = 0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, ns, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, ne, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);

  cgsize_t nsendtot[ownerofsubmesh.size()];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k]==me){
      int kk=submeshlocid[k];
      std::map<idx_t, std::set<idx_t> >:: iterator it;
      nsendtot[k] = 0;
      for(it=submeshesowned[kk].elemstosend.begin(); it!=submeshesowned[kk].elemstosend.end(); it++){
        nsendtot[k] += it->second.size();
      }
    }
    else{
      nsendtot[k] = 0; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, nsendtot, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);

  int cgzones[ownerofsubmesh.size()][5];
  int index_array[ownerofsubmesh.size()];

  cgsize_t estart, eend;
  cgsize_t isize[1][3];


  //Phase 1
  for(int k=0; k<ownerofsubmesh.size(); k++){
    int index_bc;
    isize[0][0]=ns[k];
    isize[0][1]=ne[k];
    isize[0][2]=0;
    if(cg_zone_write(index_file,index_base,zms[k].c_str(),*isize,CGNS_ENUMV(Unstructured),&cgzones[k][0]) ) cg_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateX",&cgzones[k][1])) cgp_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateY",&cgzones[k][2])) cgp_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateZ",&cgzones[k][3])) cgp_error_exit();
    estart = 1;
    eend = ne[k];
    if (cgp_section_write(index_file, index_base, cgzones[k][0], "Elements", elementType, estart, eend, 0, &cgzones[k][4])) cgp_error_exit();

    /* if(cgp_boco_write(index_file,index_base,cgzones[k][0],"Test BC Parallel",CGNS_ENUMV(BCInflow),CGNS_ENUMV(PointList),42,&index_bc)) cg_error_exit(); */
    /* if(cg_boco_write(index_file,index_base,cgzones[k][0],"Test BC Parallel",CGNS_ENUMV(BCInflow),CGNS_ENUMV(PointList),42,NULL,&index_bc)) cg_error_exit(); */

    /* if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit(); */
    /* std::string name="ElemsToSend"; */
    /* if(cg_user_data_write(name.c_str())) cg_error_exit(); */
    /* if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit(); */
    /* if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit(); */
  }


  //Phase 2
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k] == me){
      int kk = submeshlocid[k];

      cgsize_t nnodes=submeshesowned[kk].get_nnodes();
      estart=1; eend=nnodes;
      coord = new double[nnodes];
      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,0);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][1],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,1);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][2],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,2);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][3],&estart,&eend,coord)) cgp_error_exit();

      delete [] coord;
      cgsize_t nelems=submeshesowned[kk].get_nelems();
      estart=1;eend=ne[k];
      elems = new cgsize_t[esize*nelems];
      for(idx_t i=0; i<submeshesowned[kk].get_nelems(); i++){
        for(idx_t j=0; j<esize; j++){
          elems[esize*i + j] = submeshesowned[kk].renumber_nto[submeshesowned[kk].nodes_gtl[submeshesowned[kk].get_elems(i,j)]]+1;
        }
      }
      if(cgp_elements_write_data(index_file,index_base,cgzones[k][0],cgzones[k][4],estart,eend,elems)) cgp_error_exit();
      delete [] elems;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(cgp_close(index_file)) cgp_error_exit();
}

void writeMeshPCGNS(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm(MPI_COMM_WORLD);

  int index_file, index_base, index_zone, n_bases, base, physDim, cellDim, nZones, zone, index_coord, index_bc;
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
  cgsize_t *elems, *elemstosend, *elemstorecv;
  double *coord;
  int Cx, Cy, Cz;

  /* if(cgp_mpi_comm(MPI_COMM_WORLD)) cgp_error_exit(); */



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

  std::vector<std::string> zms(ownerofsubmesh.size());
  for(int k=0; k<ownerofsubmesh.size(); k++){
    std::stringstream zss;
    zss << "River_" << k;
    zms[k] = zss.str();
  }

  int cgzones[ownerofsubmesh.size()][5];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    cgzones[k][0] = 0;
    cgzones[k][1] = 0;
    cgzones[k][2] = 0;
    cgzones[k][3] = 0;
    cgzones[k][4] = 0;
  }


  MPI_Barrier(MPI_COMM_WORLD);

  //Sequential phase
  for(int p=0; p<nprocs; p++){
    if(me==p){
      if(me==0){
        if(cg_open("Mesh_Output_pcgns.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
        if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base)) cg_error_exit();
      }
      else{
        if(cg_open("Mesh_Output_pcgns.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();
      }
      index_base=1;
      for(int k=0; k<submeshesowned.size(); k++){
        /* cgsize_t estart, eend; */
        cgsize_t isize[1][3];
        isize[0][0]=submeshesowned[k].get_nnodes();
        isize[0][1]=submeshesowned[k].get_nelems();
        isize[0][2]=0;
        std::stringstream zss;
        zss << "River_" << submeshesowned[k].submeshid;
        std::string name="ElemsToSend";
        /* if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&index_zone) ) cg_error_exit(); */
        if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&cgzones[submeshesowned[k].submeshid][0]) ) cg_error_exit();
        /* std::cout << submeshesowned[k].submeshid << " " << cgzones[submeshesowned[k].submeshid][0] << std::endl; */
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_ordinal_write(submeshesowned[k].submeshid));
        if(cg_user_data_write(name.c_str())) cg_error_exit();
        for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[k].elemstosend.begin(); it!=submeshesowned[k].elemstosend.end(); it++){
          if(it->second.size()>0){
            elemstosend = new cgsize_t[it->second.size()];
            cgsize_t i=0;
            for(std::set<idx_t>::iterator iter=it->second.begin(); iter!=it->second.end(); iter++){
              elemstosend[i] = submeshesowned[k].renumber_otn[*iter]+1;
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
            elemstorecv[0] = submeshesowned[k].renumber_otn[itglobloc]+1;
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
            if(cg_boco_write(index_file,index_base,cgzones[submeshesowned[k].submeshid][0],submeshesowned[k].boundary_conditions_names[bc].c_str(),submeshesowned[k].boundary_conditions_types[bc],CGNS_ENUMV(PointList),submeshesowned[k].boundary_conditions[bc].size(),bcnodes,&index_bc)) cg_error_exit();
            cg_boco_gridlocation_write(index_file,index_base,cgzones[submeshesowned[k].submeshid][0],index_bc,CGNS_ENUMV(Vertex));
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
            if(dimensions==submeshesowned[k].get_nnodes()){
              location = CGNS_ENUMV(Vertex);
            }
            if(dimensions==submeshesowned[k].get_nelems()){
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

  MPI_Allreduce(MPI_IN_PLACE, cgzones, ownerofsubmesh.size()*5, MPI_INT, MPI_SUM, comm);

  if(cgp_pio_mode(CGP_INDEPENDENT)) cgp_error_exit();
  if(cgp_open("Mesh_Output_pcgns.cgns",CG_MODE_MODIFY,&index_file)) cgp_error_exit();
  index_base=1;

  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  std::vector<std::string> ems(ownerofsubmesh.size());

  int ns[ownerofsubmesh.size()];
  int ne[ownerofsubmesh.size()];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k] == me){
      ns[k] = submeshesowned[submeshlocid[k]].get_nnodes();
      ne[k] = submeshesowned[submeshlocid[k]].get_nelems();
    }
    else{
      ns[k] = 0;
      ne[k] = 0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, ns, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, ne, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);

  cgsize_t nsendtot[ownerofsubmesh.size()];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k]==me){
      int kk=submeshlocid[k];
      std::map<idx_t, std::set<idx_t> >:: iterator it;
      nsendtot[k] = 0;
      for(it=submeshesowned[kk].elemstosend.begin(); it!=submeshesowned[kk].elemstosend.end(); it++){
        nsendtot[k] += it->second.size();
      }
    }
    else{
      nsendtot[k] = 0; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, nsendtot, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);

  int index_array[ownerofsubmesh.size()];

  cgsize_t estart, eend;
  cgsize_t isize[1][3];

  //Corection of cgzones
  /* std::cout << "fou" << std::endl; */
  for(int ki=0; ki<ownerofsubmesh.size(); ki++){
    char dumn[50];
    int k;
    if(cg_zone_read(index_file,index_base,ki+1,dumn,*isize)) cg_error_exit();
    /* std::cout << ki << " fou " << dumn << std::endl; */
    if(cg_goto(index_file, index_base, dumn, 0, "end")) cg_error_exit();
    if(cg_ordinal_read(&k)) cg_error_exit();
    /* std::cout << "cor 1 " << k << std::endl; */
    cgzones[k][0] = ki+1;
    /* std::cout << "cor " << dumn << " " << k << " " << ki+1 << std::endl; */

  }

  //Phase 1
  for(int k=0; k<ownerofsubmesh.size(); k++){
    isize[0][0]=ns[k];
    isize[0][1]=ne[k];
    isize[0][2]=0;
    /* if(cg_zone_write(index_file,index_base,zms[k].c_str(),*isize,CGNS_ENUMV(Unstructured),&cgzones[k][0]) ) cg_error_exit(); */
    /* tmpl if(cg_zone_read(index_file, base, zone, zonename, sizes) != CG_OK) cg_get_error(); */
    /* if(cg_zone_read(index_file,index_base,&zn,dumn,*isize)) cg_error_exit(); */
    /* if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit(); */
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateX",&cgzones[k][1])) cgp_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateY",&cgzones[k][2])) cgp_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateZ",&cgzones[k][3])) cgp_error_exit();
    estart = 1;
    eend = ne[k];
    if (cgp_section_write(index_file, index_base, cgzones[k][0], "Elements", elementType, estart, eend, 0, &cgzones[k][4])) cgp_error_exit();

    /* if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit(); */
    /* std::string name="ElemsToSend"; */
    /* if(cg_user_data_write(name.c_str())) cg_error_exit(); */
    /* if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit(); */
    /* if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit(); */
  }

  MPI_Barrier(MPI_COMM_WORLD);
  /* std::cout << "phase 1 ok" << std::endl; */
  //Phase 2
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k] == me){
      int kk = submeshlocid[k];

      cgsize_t nnodes=submeshesowned[kk].get_nnodes();
      estart=1; eend=nnodes;
      coord = new double[nnodes];
      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,0);
      }
      /* std::cout << index_base << " " << cgzones[k][0] << " " << cgzones[k][1] << " " << estart << " " << eend << std::endl; */
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][1],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,1);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][2],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,2);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][3],&estart,&eend,coord)) cgp_error_exit();

      delete [] coord;
      cgsize_t nelems=submeshesowned[kk].get_nelems();
      estart=1;eend=ne[k];
      elems = new cgsize_t[esize*nelems];
      for(idx_t i=0; i<submeshesowned[kk].get_nelems(); i++){
        for(idx_t j=0; j<esize; j++){
          elems[esize*i + j] = submeshesowned[kk].renumber_nto[submeshesowned[kk].nodes_gtl[submeshesowned[kk].get_elems(i,j)]]+1;
        }
      }
      if(cgp_elements_write_data(index_file,index_base,cgzones[k][0],cgzones[k][4],estart,eend,elems)) cgp_error_exit();
      delete [] elems;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(cgp_close(index_file)) cgp_error_exit();
  MPI_Barrier(MPI_COMM_WORLD);



}

void writeMeshPCGNS_ch(std::vector<submesh> &submeshesowned, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofsubmesh){
  idx_t nprocs, me;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  MPI_Comm comm(MPI_COMM_WORLD);

  int index_file, index_base, index_zone, n_bases, base, physDim, cellDim, nZones, zone, index_coord;
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
  cgsize_t *elems, *elemstosend, *elemstorecv;
  double *coord;
  int Cx, Cy, Cz;

  /* if(cgp_mpi_comm(MPI_COMM_WORLD)) cgp_error_exit(); */
  if(cgp_pio_mode(CGP_INDEPENDENT)) cgp_error_exit();
  if(cgp_open("Mesh_Output_pcgns_ch.cgns",CG_MODE_WRITE,&index_file)) cgp_error_exit();


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

  if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base) || cg_goto(index_file, index_base,"end")) cg_error_exit();


  std::map<idx_t, idx_t> submeshlocid;
  for(idx_t i=0;i<submeshesowned.size();i++){
    submeshlocid.insert(std::make_pair(submeshesowned[i].submeshid, i));
  }

  std::vector<std::string> zms(ownerofsubmesh.size());
  std::vector<std::string> ems(ownerofsubmesh.size());
  for(int k=0; k<ownerofsubmesh.size(); k++){
    std::stringstream zss;
    zss << "River_" << k;
    zms[k] = zss.str();
  }

  int ns[ownerofsubmesh.size()];
  int ne[ownerofsubmesh.size()];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k] == me){
      ns[k] = submeshesowned[submeshlocid[k]].get_nnodes();
      ne[k] = submeshesowned[submeshlocid[k]].get_nelems();
    }
    else{
      ns[k] = 0;
      ne[k] = 0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, ns, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, ne, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);

  cgsize_t nsendtot[ownerofsubmesh.size()];
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k]==me){
      int kk=submeshlocid[k];
      std::map<idx_t, std::set<idx_t> >:: iterator it;
      nsendtot[k] = 0;
      for(it=submeshesowned[kk].elemstosend.begin(); it!=submeshesowned[kk].elemstosend.end(); it++){
        nsendtot[k] += it->second.size();
      }
    }
    else{
      nsendtot[k] = 0; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, nsendtot, ownerofsubmesh.size(), MPI_INT, MPI_SUM, comm);

  int cgzones[ownerofsubmesh.size()][5];
  int index_array[ownerofsubmesh.size()];

  cgsize_t estart, eend;
  cgsize_t isize[1][3];


  //Phase 1
  for(int k=0; k<ownerofsubmesh.size(); k++){
    isize[0][0]=ns[k];
    isize[0][1]=ne[k];
    isize[0][2]=0;
    if(cg_zone_write(index_file,index_base,zms[k].c_str(),*isize,CGNS_ENUMV(Unstructured),&cgzones[k][0]) ) cg_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateX",&cgzones[k][1])) cgp_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateY",&cgzones[k][2])) cgp_error_exit();
    if(cgp_coord_write(index_file,index_base,cgzones[k][0],CGNS_ENUMV(RealDouble),"CoordinateZ",&cgzones[k][3])) cgp_error_exit();
    estart = 1;
    eend = ne[k];
    if (cgp_section_write(index_file, index_base, cgzones[k][0], "Elements", elementType, estart, eend, 0, &cgzones[k][4])) cgp_error_exit();

  }

  int nnei, fnei, snei;

  //Phase 2
  for(int k=0; k<ownerofsubmesh.size(); k++){
    if(ownerofsubmesh[k] == me){
      int kk = submeshlocid[k];

      cgsize_t nnodes=submeshesowned[kk].get_nnodes();
      estart=1; eend=nnodes;
      coord = new double[nnodes];
      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,0);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][1],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,1);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][2],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<submeshesowned[kk].get_nnodes(); i++){
        coord[i] = submeshesowned[kk].get_nodes(i,2);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][3],&estart,&eend,coord)) cgp_error_exit();

      delete [] coord;
      cgsize_t nelems=submeshesowned[kk].get_nelems();
      estart=1;eend=ne[k];
      elems = new cgsize_t[esize*nelems];
      for(idx_t i=0; i<submeshesowned[kk].get_nelems(); i++){
        for(idx_t j=0; j<esize; j++){
          elems[esize*i + j] = submeshesowned[kk].renumber_nto[submeshesowned[kk].nodes_gtl[submeshesowned[kk].get_elems(i,j)]]+1;
        }
      }
      if(cgp_elements_write_data(index_file,index_base,cgzones[k][0],cgzones[k][4],estart,eend,elems)) cgp_error_exit();
      delete [] elems;

      //Ghost cells informations
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      std::string name="ElemsToSend";
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      int tmp_nnei=0;
      for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[kk].elemstosend.begin(); it!=submeshesowned[kk].elemstosend.end(); it++){
        if(it->second.size()>0){
          tmp_nnei++;
        }
      }
      nnei = tmp_nnei;
      MPI_Bcast(&nnei, 1, MPI_INT, me, comm);
      for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[kk].elemstosend.begin(); it!=submeshesowned[kk].elemstosend.end(); it++){
        if(it->second.size()>0){
          fnei = it->first;
          snei = it->second.size();
          MPI_Bcast(&fnei, 1, MPI_INT, me, comm);
          MPI_Bcast(&snei, 1, MPI_INT, me, comm);
          elemstosend = new cgsize_t[it->second.size()];
          cgsize_t i=0;
          for(std::set<idx_t>::iterator iter=it->second.begin(); iter!=it->second.end(); iter++){
            elemstosend[i] = submeshesowned[kk].renumber_otn[*iter]+1;
            i++;
          }
          MPI_Bcast(elemstosend, snei, MPI_INT, me, comm);
          std::stringstream ssname;
          ssname << it->first;
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit();
          if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
          if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
          /* std::cout << it->second.size() << std::endl; */
          if(cg_ptset_write(CGNS_ENUMV(PointList), it->second.size(), elemstosend)) cg_error_exit();
          if(cg_ordinal_write(it->first));
          delete [] elemstosend;
        }
      }
      name="ElemsToRecv";
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      for(std::map<idx_t, std::set<idx_t> >::iterator it=submeshesowned[kk].elemstorecv.begin(); it!=submeshesowned[kk].elemstorecv.end(); it++){
        if(it->second.size()>0){
          fnei = it->first;
          snei = it->second.size();
          MPI_Bcast(&fnei, 1, MPI_INT, me, comm);
          MPI_Bcast(&snei, 1, MPI_INT, me, comm);
          std::stringstream ssname;
          ssname << it->first;
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToRecv", 0, "end")) cg_error_exit();
          if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToRecv", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
          if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
          elemstorecv = new cgsize_t[2];
          std::set<idx_t>::iterator iter=it->second.begin();
          idx_t itglobloc = submeshesowned[kk].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*iter]];
          elemstorecv[0] = submeshesowned[kk].renumber_otn[itglobloc]+1;
          elemstorecv[1] = it->second.size();
          MPI_Bcast(elemstorecv, 2, MPI_INT, me, comm);
          if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
          if(cg_ordinal_write(it->first));
          delete [] elemstorecv;
        }
      }
    }
    else{
      //Ghost cells informations
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      std::string name="ElemsToSend";
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      MPI_Bcast(&nnei, 1, MPI_INT, ownerofsubmesh[k], comm);
      for(int knei=0; knei<nnei; knei++){
        MPI_Bcast(&fnei, 1, MPI_INT, ownerofsubmesh[k], comm);
        MPI_Bcast(&snei, 1, MPI_INT, ownerofsubmesh[k], comm);
        elemstosend = new cgsize_t[snei];
        MPI_Bcast(elemstosend, snei, MPI_INT, ownerofsubmesh[k], comm);
        std::stringstream ssname;
        ssname << fnei;
        if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit();
        if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
        if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
        if(cg_ptset_write(CGNS_ENUMV(PointList), snei, elemstosend)) cg_error_exit();
        if(cg_ordinal_write(fnei));
        delete [] elemstosend;
      }
      name="ElemsToRecv";
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      for(int knei=0; knei<nnei; knei++){
        MPI_Bcast(&fnei, 1, MPI_INT, ownerofsubmesh[k], comm);
        MPI_Bcast(&snei, 1, MPI_INT, ownerofsubmesh[k], comm);
        std::stringstream ssname;
        ssname << fnei;
        if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToRecv", 0, "end")) cg_error_exit();
        if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
        if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToRecv", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
        elemstorecv = new cgsize_t[2];
        MPI_Bcast(elemstorecv, 2, MPI_INT, ownerofsubmesh[k], comm);
        if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
        if(cg_ordinal_write(fnei));
        delete [] elemstorecv;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(cgp_close(index_file)) cgp_error_exit();
}
