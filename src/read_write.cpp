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

//Global variables to store name an type of all boundary conditions
//The issue with local variables is that it is diffucult to send BCType_t via MPI
std::vector<std::string> boundary_conditions_names;
std::vector<CGNS_ENUMV(BCType_t)> boundary_conditions_types;

std::vector<std::string>  ud_names;
std::vector<std::vector<std::string> > ar_names;

void parallel_read_mesh_cgns(idx_t*& elmdist, idx_t*& eptr, idx_t*& eind, idx_t*& epart, idx_t& esize, idx_t& dim, const idx_t numberingstart, std::string filename, MPI_Comm comm){

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
  epart = new idx_t[nelems];
  eptr = new idx_t[nelems+1];

  for(idx_t i=0;i<nelems;i++){
    for(idx_t j=0;j<esize;j++){
      eind[i*esize+j] = (idx_t) (conn[i*esize+j]-1);
    }
  }
  for(idx_t i=0;i<=nelems;i++){
    eptr[i] = esize*i;
  }

  if (cg_close(index_file)) cg_error_exit();

}

void read_boundary_conditions(std::vector<partition> &parts, std::string filename){
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
  boundary_conditions_names.resize(nBocos);
  boundary_conditions_types.resize(nBocos);
  for(int k=0; k<parts.size(); k++){
    parts[k].boundary_conditions.resize(nBocos);
  }
  for(int boco=1; boco<=nBocos; boco++)
  {
    if(cg_boco_info(index_file, base, zone, boco, boconame, &bocoType,
          &ptsetType, &nBCNodes, &normalIndex,
          &normListFlag, &normDataType, &nDataSet) != CG_OK) cg_get_error();
    bcnodes = new cgsize_t[nBCNodes];
    CGNS_ENUMV(GridLocation_t) location;
    if(cg_boco_gridlocation_read(index_file, base, zone, boco, &location) != CG_OK) cg_get_error();
    /* assert(location==CGNS_ENUMV(Vertex)); */
    if(location!=CGNS_ENUMV(Vertex)){
      std::cout << "Boundary condition not on vertex is ignored" << std::endl;
    }
    else{
      if(cg_boco_read(index_file, base, zone, boco, bcnodes,
            NULL) != CG_OK) cg_get_error();

      boundary_conditions_names[boco-1] = boconame;
      boundary_conditions_types[boco-1] = bocoType;

      for(cgsize_t i=0; i<nBCNodes; i++){
        bcnodes[i] --;
        for(int k=0; k<parts.size(); k++){
          if(parts[k].nodes_gtl.count(bcnodes[i]) != 0){
            parts[k].boundary_conditions[boco-1].insert(bcnodes[i]);
          }
        }
      }
    }

    delete [] bcnodes;
  }
}

void write_cgns_single(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition){
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

  std::map<idx_t, idx_t> partitionlocid;
  for(idx_t i=0;i<parts.size();i++){
    partitionlocid.insert(std::make_pair(parts[i].partitionid, i));
  }

  std::vector<std::string> zms(ownerofpartition.size());
  for(int k=0; k<ownerofpartition.size(); k++){
    std::stringstream zss;
    zss << "Zone_" << k;
    zms[k] = zss.str();
  }

  if(me==0){
    if(cg_open("Mesh_Output.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
    if(cg_base_write(index_file,"Base",icelldim,iphysdim,&index_base)) cg_error_exit();
    if(cg_close(index_file)) cg_error_exit();
  }
  MPI_Barrier(comm);

  cgsize_t cgzones[ownerofpartition.size()][4];

  for(int p=0; p<nprocs; p++){
    if(me==p){
      if(cg_open("Mesh_Output.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();
      index_base=1;
      for(int k=0; k<parts.size(); k++){
        cgsize_t isize[1][3];
        isize[0][0]=parts[k].get_nnodes();
        isize[0][1]=parts[k].get_nelems();
        isize[0][2]=0;
        std::stringstream zss;
        zss << "Zone_" << parts[k].partitionid;
        if(cg_zone_write(index_file,index_base,zss.str().c_str(),*isize,CGNS_ENUMV(Unstructured),&index_zone) ) cg_error_exit();
        coord = new double[parts[k].get_nnodes()];
        for(idx_t i=0; i<parts[k].get_nnodes(); i++){
          coord[i] = parts[k].get_nodes(i,0);
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateX", coord, &Cx)) cg_error_exit();
        for(idx_t i=0; i<parts[k].get_nnodes(); i++){
          coord[i] = parts[k].get_nodes(i,1);
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateY", coord, &Cy)) cg_error_exit();
        for(idx_t i=0; i<parts[k].get_nnodes(); i++){
          coord[i] = parts[k].get_nodes(i,2);
        }
        if(cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble), "CoordinateZ", coord, &Cz)) cg_error_exit();
        delete [] coord;
        int nstart=1, nend=parts[k].get_nelems();
        elems = new cgsize_t[esize*parts[k].get_nelems()];
        for(idx_t i=0; i<parts[k].get_nelems(); i++){
          for(idx_t j=0; j<esize; j++){
            int oldi = parts[k].renumber_nto[i];
            elems[esize*i + j] = parts[k].nodes_gtl[parts[k].get_elems(oldi,j)]+1;
            /* elems[esize*i + j] = parts[k].renumber_otn[parts[k].nodes_gtl[parts[k].get_elems(i,j)]]+1; */
          }
        }
        if(cg_section_write(index_file,index_base,index_zone,"Elements",elementType,nstart,nend,0,elems,&index_section)) cg_error_exit();
        delete [] elems;
        std::string name="ElemsToSend";
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_user_data_write(name.c_str())) cg_error_exit();
        for(std::map<idx_t, std::set<idx_t> >::iterator it=parts[k].elemstosend.begin(); it!=parts[k].elemstosend.end(); it++){
          if(it->second.size()>0){
            elemstosend = new cgsize_t[it->second.size()];
            cgsize_t i=0;
            for(std::set<idx_t>::iterator iter=it->second.begin(); iter!=it->second.end(); iter++){
              elemstosend[i] = parts[k].renumber_otn[*iter]+1;
              i++;
            }


            /* //BS as part of new renumbering correction */
            for(idx_t ie=0; ie < it->second.size(); ie++){
              for(idx_t je=0; je < it->second.size()-ie-1; je++){
                /* if(parts[k].elems_ltg[parts[k].renumber_nto[elemstosend[je]-1]] > parts[k].elems_ltg[parts[k].renumber_nto[elemstosend[je+1]-1]]) { */
                /* if(parts[k].elems_ltg[elemstosend[je]-1] > parts[k].elems_ltg[elemstosend[je+1]-1]) { */
                if(parts[k].elems_ltg[parts[k].renumber_nto[elemstosend[je]-1]] > parts[k].elems_ltg[parts[k].renumber_nto[elemstosend[je+1]-1]]) {
                  idx_t dummy = elemstosend[je];
                  elemstosend[je] = elemstosend[je+1];
                  elemstosend[je+1] = dummy;
                }
              }
              }


            std::stringstream ssname;
            ssname << it->first;
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit();
            if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToSend", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
            if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
            if(cg_ptset_write(CGNS_ENUMV(PointList), it->second.size(), elemstosend)) cg_error_exit();
            if(cg_ordinal_write(it->first));
            delete [] elemstosend;
          }
        }
        name="ElemsToRecv";
        if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_user_data_write(name.c_str())) cg_error_exit();
        for(std::map<idx_t, std::set<idx_t> >::iterator it=parts[k].elemstorecv.begin(); it!=parts[k].elemstorecv.end(); it++){
          if(it->second.size()>0){
            std::stringstream ssname;
            ssname << it->first;
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToRecv", 0, "end")) cg_error_exit();
            if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "ElemsToRecv", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
            if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
            elemstorecv = new cgsize_t[2];
            std::set<idx_t>::iterator iter=it->second.begin();
            idx_t itglobloc = parts[k].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*iter]];
            elemstorecv[0] = parts[k].renumber_otn[itglobloc]+1;
            elemstorecv[1] = it->second.size();
            if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
            if(cg_ordinal_write(it->first));
            delete [] elemstorecv;
          }
        }

        //Write boundary conditions
        for(int bc=0; bc<parts[k].boundary_conditions.size(); bc++){
          if(parts[k].boundary_conditions[bc].size()>0){
            cgsize_t bcnodes[parts[k].boundary_conditions[bc].size()];
            cgsize_t nbc=0;
            for(std::set<idx_t>::iterator it=parts[k].boundary_conditions[bc].begin();
                it!=parts[k].boundary_conditions[bc].end();it++){
              bcnodes[nbc] = parts[k].nodes_gtl[*it] ;
              nbc++;
            }
            if(cg_boco_write(index_file,index_base,index_zone,boundary_conditions_names[bc].c_str(),boundary_conditions_types[bc],CGNS_ENUMV(PointList),parts[k].boundary_conditions[bc].size(),bcnodes,&index_bc)) cg_error_exit();
            cg_boco_gridlocation_write(index_file,index_base,index_zone,index_bc,CGNS_ENUMV(Vertex));
          }
        }

        int nuserdata = ud_names.size();
        for(int nud=1; nud<=nuserdata; nud++){
          int narrays = ar_names[nud-1].size();
          if(narrays>0){
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, "end")) cg_error_exit();
            if(cg_user_data_write(ud_names[nud-1].c_str())) cg_error_exit();
            cgsize_t dimensions=parts[k].arrays[nud-1][0].size();
            CGNS_ENUMV(GridLocation_t) location;
            if(dimensions==parts[k].get_nnodes()){
              location = CGNS_ENUMV(Vertex);
            }
            if(dimensions==parts[k].get_nelems()){
              location = CGNS_ENUMV(CellCenter);
            }
            if(cg_goto(index_file, index_base, zss.str().c_str(), 0, ud_names[nud-1].c_str(), 0, "end")) cg_error_exit();
            if(cg_gridlocation_write(location)) cg_error_exit();
            for(int na=1; na<=narrays; na++){
              int rank=1;
              if(cg_array_write(ar_names[nud-1][na-1].c_str(), CGNS_ENUMV(RealDouble), 1, &dimensions, parts[k].arrays[nud-1][na-1].data())) cg_error_exit();

            }
          }
        }

      }
      if(cg_close(index_file)) cg_error_exit();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void read_arrays(std::vector<partition> &parts, std::string filename, MPI_Comm comm){
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
  for(int nud=1; nud<=nuserdata; nud++){
    char udname[40];
    if(cg_user_data_read(nud, udname)) cg_error_exit();
    if(cg_goto(index_file, index_base, zonename, 0, udname, 0, "end")) cg_error_exit();
    if(cg_narrays(&narrays)) cg_error_exit();
    for(int na=1; na<=narrays; na++){
      char aname[40];
      CGNS_ENUMV(DataType_t) dataType;
      int rank;
      cgsize_t dimensions;
      CGNS_ENUMV(GridLocation_t) location;
      if(cg_gridlocation_read(&location)) cg_error_exit();
      assert(location==CGNS_ENUMV(CellCenter) or location==CGNS_ENUMV(Vertex));
      if(cg_array_info(na,aname,&dataType,&rank,&dimensions)) cg_error_exit();
      if(dataType == CGNS_ENUMV(RealDouble)){
        double *array = new double[dimensions];
        if(cg_array_read(na, array)) cg_error_exit();
        std::vector<idx_t> nloc(parts.size(), 0);
        //Resize if first time
        if(nud==1 and na==1){
            ud_names.resize(nuserdata);
            ar_names.resize(nuserdata, std::vector<std::string>(narrays));
          for(idx_t k=0; k<parts.size(); k++){
            /* parts[k].ud_names.resize(nuserdata); */
            /* parts[k].ar_names.resize(nuserdata, std::vector<std::string>(narrays)); */
            parts[k].arrays.resize(nuserdata, std::vector<std::vector<double> >(narrays));
          }
        }

        ud_names[nud-1] = std::string(udname);
        ar_names[nud-1][na-1] = std::string(aname);


        for(idx_t i=0; i<dimensions; i++){
          for(idx_t k=0; k<parts.size(); k++){
            if(location=CGNS_ENUMV(CellCenter)){
              if(parts[k].elems_gtl.count(i)!=0){
                if(parts[k].arrays[nud-1][na-1].size() == 0) {
                  /* parts[k].ud_names[nud-1] = std::string(udname); */
                  /* parts[k].ar_names[nud-1][na-1] = std::string(aname); */
                  parts[k].arrays[nud-1][na-1].resize(parts[k].get_nelems());
                }
                int iloc = parts[k].elems_gtl[i];
                parts[k].arrays[nud-1][na-1][iloc] = array[i];
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

void write_wo_recv_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim){
  for(idx_t k=0; k<parts.size(); k++){
    idx_t indrenum;
    idx_t nNodes = parts[k].get_nnodes();
    idx_t nElems = parts[k].get_nelems();
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << parts[k].partitionid << "_output_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << parts[k].get_nodes(i,0) << " " << parts[k].get_nodes(i,1) << " " << parts[k].get_nodes(i,2) << std::endl;
    }

    idx_t totelemsrecv=0;
    for(std::map<idx_t, std::set<idx_t> >::iterator it=parts[k].elemstorecv.begin(); it!= parts[k].elemstorecv.end(); it++){
      totelemsrecv += it->second.size();     
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems-totelemsrecv << " " << (nElems-totelemsrecv)*(esize+1) << std::endl;
    for(idx_t i=0; i<nElems; i++){
      indrenum = parts[k].renumber_nto[i];
      if(indrenum<nElems-totelemsrecv){
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << parts[k].nodes_gtl[parts[k].get_elems(indrenum,p)] << " ";
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

void write_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim){
  for(idx_t k=0; k<parts.size(); k++){
    idx_t indrenum;
    idx_t nNodes = parts[k].get_nnodes();
    idx_t nElems = parts[k].get_nelems();
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << parts[k].partitionid << "_output_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << parts[k].get_nodes(i,0) << " " << parts[k].get_nodes(i,1) << " " << parts[k].get_nodes(i,2) << std::endl;
    }

    outFile << std::endl;
    outFile << "CELLS " << nElems << " " << nElems*(esize+1) << std::endl;
    for(idx_t i=0; i<nElems; i++){
      outFile << esize << " ";
      indrenum = parts[k].renumber_nto[i];
      for(idx_t p=0; p<esize; p++){
        outFile << parts[k].nodes_gtl[parts[k].get_elems(indrenum,p)] << " ";
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

void write_recv_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim){

  for(idx_t k=0; k<parts.size(); k++){
    idx_t nNodes = parts[k].get_nnodes();
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << parts[k].partitionid << "_output_recv_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << parts[k].get_nodes(i,0) << " " << parts[k].get_nodes(i,1) << " " << parts[k].get_nodes(i,2) << std::endl;
    }



    idx_t nbelemsrecv(0);
    std::set<idx_t>::iterator iter;
    for(iter=parts[k].potentialneighbors.begin();iter!=parts[k].potentialneighbors.end();iter++){
      nbelemsrecv +=  parts[k].elemstorecv[*iter].size();
    }

    outFile << std::endl;
    outFile << "CELLS " << nbelemsrecv << " " << nbelemsrecv*(esize+1) << std::endl;
    std::set<idx_t>::iterator it;
    for(iter=parts[k].potentialneighbors.begin();iter!=parts[k].potentialneighbors.end();iter++){
      for(it=parts[k].elemstorecv[*iter].begin(); it!=parts[k].elemstorecv[*iter].end(); it++){
        idx_t itloc = parts[k].elems_gtl[g_potentialneighbors[*iter].elems_ltg[*it]];
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << parts[k].nodes_gtl[parts[k].get_elems(itloc,p)] << " ";
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
      outFile << parts[k].get_nodes(i,2) << std::endl;
    }
  }
}

void write_send_vtk(std::vector<partition> &parts, idx_t esize, idx_t dim){
  for(idx_t k=0; k<parts.size(); k++){
    idx_t nNodes = parts[k].get_nnodes();
    std::stringstream prefixedOutputFilename;
    prefixedOutputFilename << parts[k].partitionid << "_output_send_paraview.vtk";
    std::ofstream outFile(prefixedOutputFilename.str());
    outFile << std::setprecision(16);
    outFile << "# vtk DataFile Version 1.0" << std::endl;
    outFile << "Genrated from ParMetis Ghost" << std::endl;
    outFile << "ASCII" << std::endl << std::endl;
    outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    outFile << "POINTS " << nNodes << " float" << std::endl;
    for(idx_t i=0; i<nNodes; i++){
      outFile << parts[k].get_nodes(i,0) << " " << parts[k].get_nodes(i,1) << " " << parts[k].get_nodes(i,2) << std::endl;
    }

    idx_t nbelemssend(0);
    std::set<idx_t>::iterator iter;
    for(iter=parts[k].potentialneighbors.begin();iter!=parts[k].potentialneighbors.end();iter++){
      nbelemssend +=  parts[k].elemstosend[*iter].size();
    }

    outFile << std::endl;
    outFile << "CELLS " << nbelemssend << " " << nbelemssend*(esize+1) << std::endl;
    std::set<idx_t>::iterator it;
    for(iter=parts[k].potentialneighbors.begin();iter!=parts[k].potentialneighbors.end();iter++){
      for(it=parts[k].elemstosend[*iter].begin(); it!=parts[k].elemstosend[*iter].end(); it++){
        outFile << esize << " ";
        for(idx_t p=0; p<esize; p++){
          outFile << parts[k].nodes_gtl[parts[k].get_elems(*it,p)] << " ";
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
      outFile << parts[k].get_nodes(i,2) << std::endl;
    }
  }
}

void read_nodes_cgns(std::vector<partition> &parts, std::string filename, MPI_Comm comm){
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
  double *x, *y, *z;

  //Resize each subdomain nodes vector
  for(idx_t k=0; k<parts.size(); k++){
    parts[k].nodes.resize(3*parts[k].nodesindex.size(), 0.);
  }

  //Construct a global_nodes vector so that we only check once if we own a node
  std::map<idx_t, std::set<idx_t>> global_nodes;
  for(idx_t k=0; k<parts.size(); k++){
    for(std::unordered_set<idx_t>::iterator it=parts[k].nodesindex.begin();
        it!=parts[k].nodesindex.end(); it++){
      global_nodes[*it].insert(k);
    }
  }

  //Read nodes in batches
  std::vector<idx_t> nloc(parts.size(), 0);
  int nbatches=nprocs;
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
          parts[k].get_nodes(nloc[k],0) = x[j];
          parts[k].get_nodes(nloc[k],1) = y[j];
          parts[k].get_nodes(nloc[k],2) = z[j];
          parts[k].nodes_ltg.insert({nloc[k], i});
          parts[k].nodes_gtl.insert({i, nloc[k]});
          nloc[k] += 1;
        }
      }
      j++;
    }
  }
  if (cg_close(index_file)) cg_error_exit();
}

void write_pcgns_without_send_recv_info(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition){
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


  std::map<idx_t, idx_t> partitionlocid;
  for(idx_t i=0;i<parts.size();i++){
    partitionlocid.insert(std::make_pair(parts[i].partitionid, i));
  }

  std::vector<std::string> zms(ownerofpartition.size());
  std::vector<std::string> ems(ownerofpartition.size());
  for(int k=0; k<ownerofpartition.size(); k++){
    std::stringstream zss;
    zss << "Zone_" << k;
    zms[k] = zss.str();
  }

  int ns[ownerofpartition.size()];
  int ne[ownerofpartition.size()];
  for(int k=0; k<ownerofpartition.size(); k++){
    if(ownerofpartition[k] == me){
      ns[k] = parts[partitionlocid[k]].get_nnodes();
      ne[k] = parts[partitionlocid[k]].get_nelems();
    }
    else{
      ns[k] = 0;
      ne[k] = 0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, ns, ownerofpartition.size(), MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, ne, ownerofpartition.size(), MPI_INT, MPI_SUM, comm);

  cgsize_t nsendtot[ownerofpartition.size()];
  for(int k=0; k<ownerofpartition.size(); k++){
    if(ownerofpartition[k]==me){
      int kk=partitionlocid[k];
      std::map<idx_t, std::set<idx_t> >:: iterator it;
      nsendtot[k] = 0;
      for(it=parts[kk].elemstosend.begin(); it!=parts[kk].elemstosend.end(); it++){
        nsendtot[k] += it->second.size();
      }
    }
    else{
      nsendtot[k] = 0; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, nsendtot, ownerofpartition.size(), MPI_INT, MPI_SUM, comm);

  int cgzones[ownerofpartition.size()][5];
  int index_array[ownerofpartition.size()];

  cgsize_t estart, eend;
  cgsize_t isize[1][3];


  //Phase 1
  for(int k=0; k<ownerofpartition.size(); k++){
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
  }


  //Phase 2
  for(int k=0; k<ownerofpartition.size(); k++){
    if(ownerofpartition[k] == me){
      int kk = partitionlocid[k];

      cgsize_t nnodes=parts[kk].get_nnodes();
      estart=1; eend=nnodes;
      coord = new double[nnodes];
      for(idx_t i=0; i<parts[kk].get_nnodes(); i++){
        coord[i] = parts[kk].get_nodes(i,0);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][1],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<parts[kk].get_nnodes(); i++){
        coord[i] = parts[kk].get_nodes(i,1);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][2],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<parts[kk].get_nnodes(); i++){
        coord[i] = parts[kk].get_nodes(i,2);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][3],&estart,&eend,coord)) cgp_error_exit();

      delete [] coord;
      cgsize_t nelems=parts[kk].get_nelems();
      estart=1;eend=ne[k];
      elems = new cgsize_t[esize*nelems];
      for(idx_t i=0; i<parts[kk].get_nelems(); i++){
        for(idx_t j=0; j<esize; j++){
          int oldi = parts[kk].renumber_nto[i];
          elems[esize*i + j] = parts[kk].nodes_gtl[parts[kk].get_elems(oldi,j)]+1;
          /* elems[esize*i + j] = parts[kk].renumber_nto[parts[kk].nodes_gtl[parts[kk].get_elems(i,j)]]+1; */
        }
      }
      if(cgp_elements_write_data(index_file,index_base,cgzones[k][0],cgzones[k][4],estart,eend,elems)) cgp_error_exit();
      delete [] elems;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(cgp_close(index_file)) cgp_error_exit();
}


void write_pcgns_hybird_with_send_recv_info(std::vector<partition> &parts, idx_t esize, idx_t dim, std::vector<idx_t> &ownerofpartition){
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


  std::map<idx_t, idx_t> partitionlocid;
  for(idx_t i=0;i<parts.size();i++){
    partitionlocid.insert(std::make_pair(parts[i].partitionid, i));
  }

  std::vector<std::string> zms(ownerofpartition.size());
  std::vector<std::string> ems(ownerofpartition.size());
  for(int k=0; k<ownerofpartition.size(); k++){
    std::stringstream zss;
    zss << "Zone_" << k;
    zms[k] = zss.str();
  }

  int ns[ownerofpartition.size()];
  int ne[ownerofpartition.size()];
  for(int k=0; k<ownerofpartition.size(); k++){
    if(ownerofpartition[k] == me){
      ns[k] = parts[partitionlocid[k]].get_nnodes();
      ne[k] = parts[partitionlocid[k]].get_nelems();
    }
    else{
      ns[k] = 0;
      ne[k] = 0;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, ns, ownerofpartition.size(), MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, ne, ownerofpartition.size(), MPI_INT, MPI_SUM, comm);

  cgsize_t nsendtot[ownerofpartition.size()];
  for(int k=0; k<ownerofpartition.size(); k++){
    if(ownerofpartition[k]==me){
      int kk=partitionlocid[k];
      std::map<idx_t, std::set<idx_t> >:: iterator it;
      nsendtot[k] = 0;
      for(it=parts[kk].elemstosend.begin(); it!=parts[kk].elemstosend.end(); it++){
        nsendtot[k] += it->second.size();
      }
    }
    else{
      nsendtot[k] = 0; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, nsendtot, ownerofpartition.size(), MPI_INT, MPI_SUM, comm);

  int cgzones[ownerofpartition.size()][5];
  int index_array[ownerofpartition.size()];

  cgsize_t estart, eend;
  cgsize_t isize[1][3];


  //Phase 1
  for(int k=0; k<ownerofpartition.size(); k++){
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
  for(int k=0; k<ownerofpartition.size(); k++){
    if(ownerofpartition[k] == me){
      int kk = partitionlocid[k];

      cgsize_t nnodes=parts[kk].get_nnodes();
      estart=1; eend=nnodes;
      coord = new double[nnodes];
      for(idx_t i=0; i<parts[kk].get_nnodes(); i++){
        coord[i] = parts[kk].get_nodes(i,0);
      }

      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][1],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<parts[kk].get_nnodes(); i++){
        coord[i] = parts[kk].get_nodes(i,1);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][2],&estart,&eend,coord)) cgp_error_exit();

      for(idx_t i=0; i<parts[kk].get_nnodes(); i++){
        coord[i] = parts[kk].get_nodes(i,2);
      }
      if(cgp_coord_write_data(index_file,index_base,cgzones[k][0],cgzones[k][3],&estart,&eend,coord)) cgp_error_exit();

      delete [] coord;
      cgsize_t nelems=parts[kk].get_nelems();
      estart=1;eend=ne[k];
      elems = new cgsize_t[esize*nelems];
      for(idx_t i=0; i<parts[kk].get_nelems(); i++){
        for(idx_t j=0; j<esize; j++){
          int oldi = parts[kk].renumber_nto[i];
          elems[esize*i + j] = parts[kk].nodes_gtl[parts[kk].get_elems(oldi,j)]+1;
          /* elems[esize*i + j] = parts[kk].renumber_nto[parts[kk].nodes_gtl[parts[kk].get_elems(i,j)]]+1; */
        }
      }
      if(cgp_elements_write_data(index_file,index_base,cgzones[k][0],cgzones[k][4],estart,eend,elems)) cgp_error_exit();
      delete [] elems;

      //Zone ordinal number
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      if(cg_ordinal_write(k)) cg_error_exit();

      //Ghost cells informations
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      std::string name="ElemsToSend";
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      int tmp_nnei=0;
      for(std::map<idx_t, std::set<idx_t> >::iterator it=parts[kk].elemstosend.begin(); it!=parts[kk].elemstosend.end(); it++){
        if(it->second.size()>0){
          tmp_nnei++;
        }
      }
      nnei = tmp_nnei;
      MPI_Bcast(&nnei, 1, MPI_INT, me, comm);
      for(std::map<idx_t, std::set<idx_t> >::iterator it=parts[kk].elemstosend.begin(); it!=parts[kk].elemstosend.end(); it++){
        if(it->second.size()>0){
          fnei = it->first;
          snei = it->second.size();
          MPI_Bcast(&fnei, 1, MPI_INT, me, comm);
          MPI_Bcast(&snei, 1, MPI_INT, me, comm);
          elemstosend = new cgsize_t[it->second.size()];
          cgsize_t i=0;
          for(std::set<idx_t>::iterator iter=it->second.begin(); iter!=it->second.end(); iter++){
            elemstosend[i] = parts[kk].renumber_otn[*iter]+1;
            i++;
          }

          /* //BS as part of new renumbering correction */
          for(idx_t ie=0; ie < it->second.size(); ie++){
            for(idx_t je=0; je < it->second.size()-ie-1; je++){
              if(parts[kk].elems_ltg[parts[kk].renumber_nto[elemstosend[je]-1]] > parts[kk].elems_ltg[parts[kk].renumber_nto[elemstosend[je+1]-1]]) {
                idx_t dummy = elemstosend[je];
                elemstosend[je] = elemstosend[je+1];
                elemstosend[je+1] = dummy;
              }
            }
            }


          MPI_Bcast(elemstosend, snei, MPI_INT, me, comm);
          std::stringstream ssname;
          ssname << it->first;
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, "end")) cg_error_exit();
          if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToSend", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
          if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
          if(cg_ptset_write(CGNS_ENUMV(PointList), it->second.size(), elemstosend)) cg_error_exit();
          if(cg_ordinal_write(it->first));
          delete [] elemstosend;
        }
      }
      name="ElemsToRecv";
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      for(std::map<idx_t, std::set<idx_t> >::iterator it=parts[kk].elemstorecv.begin(); it!=parts[kk].elemstorecv.end(); it++){
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
          idx_t itglobloc = parts[kk].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*iter]];
          elemstorecv[0] = parts[kk].renumber_otn[itglobloc]+1;

          for(iter=it->second.begin(); iter!=it->second.end(); iter++){
            itglobloc = parts[kk].elems_gtl[g_potentialneighbors[it->first].elems_ltg[*iter]];
            elemstorecv[0] = std::min(elemstorecv[0], parts[kk].renumber_otn[itglobloc]+1);
          }

          elemstorecv[1] = it->second.size();
          MPI_Bcast(elemstorecv, 2, MPI_INT, me, comm);
          if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
          if(cg_ordinal_write(it->first));
          delete [] elemstorecv;
        }
      }

      // Boundary conditions
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      //Write boundary conditions
      int index_bc;
      for(int bc=0; bc<parts[kk].boundary_conditions.size(); bc++){
        int bcsize = parts[kk].boundary_conditions[bc].size();
        MPI_Bcast(&bcsize, 1,  MPI_INT, me, comm);
        if(parts[kk].boundary_conditions[bc].size()>0){
          /* std::cout << me << " sender " << parts[kk].boundary_conditions[bc].size() << std::endl; */
          cgsize_t bcnodes[parts[kk].boundary_conditions[bc].size()];
          cgsize_t ibc=0;
          for(std::set<idx_t>::iterator it=parts[kk].boundary_conditions[bc].begin();
              it!=parts[kk].boundary_conditions[bc].end();it++){
            bcnodes[ibc] = parts[kk].nodes_gtl[*it] + 1;
            /* if(boundary_conditions_names[bc]=="Inflow nodes"){ */
            /*   std::cout << parts[kk].partitionid << " " << *it << " " << parts[kk].nodes_gtl[*it] << std::endl; */
            /* } */
            ibc++;
          }
          
          MPI_Bcast(&bcnodes, parts[kk].boundary_conditions[bc].size(), MPI_INT, me, comm);
          int index_zone = cgzones[k][0];
          if(cg_boco_write(index_file,index_base,index_zone,boundary_conditions_names[bc].c_str(),boundary_conditions_types[bc],CGNS_ENUMV(PointList),parts[kk].boundary_conditions[bc].size(),bcnodes,&index_bc)) cg_error_exit();

          if(cg_boco_gridlocation_write(index_file,index_base,index_zone,index_bc,CGNS_ENUMV(Vertex))) cg_error_exit();

          /* std::cout << "send " << index_zone << " " << boundary_conditions_names[bc] << std::endl; */
        }
      }


      //Write user data
      int nuserdata = ud_names.size();
      for(int nud=1; nud<=nuserdata; nud++){
        int narrays = ar_names[nud-1].size();
        if(narrays>0){
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
          if(cg_user_data_write(ud_names[nud-1].c_str())) cg_error_exit();
          cgsize_t dimensions=parts[kk].arrays[nud-1][0].size();
          CGNS_ENUMV(GridLocation_t) location;
          MPI_Bcast(&dimensions, 1, MPI_INT, me, comm);
          int id_location;
          if(dimensions==parts[kk].get_nnodes()){
            location = CGNS_ENUMV(Vertex);
            id_location = 0;
          }
          if(dimensions==parts[kk].get_nelems()){
            location = CGNS_ENUMV(CellCenter);
            id_location = 1;
          }
          MPI_Bcast(&id_location, 1, MPI_INT, ownerofpartition[k], comm);
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, ud_names[nud-1].c_str(), 0, "end")) cg_error_exit();
          if(cg_gridlocation_write(location)) cg_error_exit();
          for(int na=1; na<=narrays; na++){
            int rank=1;
            MPI_Bcast(parts[kk].arrays[nud-1][na-1].data(), dimensions, MPI_DOUBLE, me, comm);
            if(cg_array_write(ar_names[nud-1][na-1].c_str(), CGNS_ENUMV(RealDouble), 1, &dimensions, parts[kk].arrays[nud-1][na-1].data())) cg_error_exit();

          }
        }

      }

    }
    else{
      //Zone ordinal number
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      if(cg_ordinal_write(k)) cg_error_exit();
      //Ghost cells informations
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      std::string name="ElemsToSend";
      if(cg_user_data_write(name.c_str())) cg_error_exit();
      MPI_Bcast(&nnei, 1, MPI_INT, ownerofpartition[k], comm);
      for(int knei=0; knei<nnei; knei++){
        MPI_Bcast(&fnei, 1, MPI_INT, ownerofpartition[k], comm);
        MPI_Bcast(&snei, 1, MPI_INT, ownerofpartition[k], comm);
        elemstosend = new cgsize_t[snei];
        MPI_Bcast(elemstosend, snei, MPI_INT, ownerofpartition[k], comm);
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
        MPI_Bcast(&fnei, 1, MPI_INT, ownerofpartition[k], comm);
        MPI_Bcast(&snei, 1, MPI_INT, ownerofpartition[k], comm);
        std::stringstream ssname;
        ssname << fnei;
        if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToRecv", 0, "end")) cg_error_exit();
        if(cg_user_data_write(ssname.str().c_str())) cg_error_exit();
        if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "ElemsToRecv", 0, ssname.str().c_str(), 0, "end")) cg_error_exit();
        if(cg_gridlocation_write(CGNS_ENUMV(CellCenter))) cg_error_exit();
        elemstorecv = new cgsize_t[2];
        MPI_Bcast(elemstorecv, 2, MPI_INT, ownerofpartition[k], comm);
        if(cg_ptset_write(CGNS_ENUMV(PointRange), 2, elemstorecv)) cg_error_exit();
        if(cg_ordinal_write(fnei));
        delete [] elemstorecv;
      }

      int index_bc;
      // Boundary conditions
      if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
      //Write boundary conditions
      int nbc=parts[0].boundary_conditions.size();
      for(int bc=0; bc<nbc; bc++){
        int bcsize=0;
        MPI_Bcast(&bcsize, 1,  MPI_INT, ownerofpartition[k], comm);
        if(bcsize>0){
          int *bcnodes = new int[bcsize];
          MPI_Bcast(bcnodes, bcsize, MPI_INT, ownerofpartition[k], comm);

          int index_zone = cgzones[k][0];

          if(cg_boco_write(index_file,index_base,index_zone,boundary_conditions_names[bc].c_str(),boundary_conditions_types[bc],CGNS_ENUMV(PointList),bcsize,bcnodes,&index_bc)) cg_error_exit();
          if(cg_boco_gridlocation_write(index_file,index_base,index_zone,index_bc,CGNS_ENUMV(Vertex))) cg_error_exit();
          /* std::cout << "recv " << index_zone << " " << boundary_conditions_names[bc] << std::endl; */
        }
      }



      //Write user data
      int nuserdata = ud_names.size();
      for(int nud=1; nud<=nuserdata; nud++){
        int narrays = ar_names[nud-1].size();
        if(narrays>0){
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, "end")) cg_error_exit();
          if(cg_user_data_write(ud_names[nud-1].c_str())) cg_error_exit();
          int dimensions, id_location;
          CGNS_ENUMV(GridLocation_t) location;
          MPI_Bcast(&dimensions, 1, MPI_INT, ownerofpartition[k], comm);
          double *array = new double[dimensions];
          MPI_Bcast(&id_location, 1, MPI_INT, ownerofpartition[k], comm);
          if(id_location==0){
            location = CGNS_ENUMV(Vertex);
          }
          if(id_location==1){
            location = CGNS_ENUMV(CellCenter);
          }
          if(cg_goto(index_file, index_base, zms[k].c_str(), 0, ud_names[nud-1].c_str(), 0, "end")) cg_error_exit();
          if(cg_gridlocation_write(location)) cg_error_exit();
          for(int na=1; na<=narrays; na++){
            int rank=1;
            MPI_Bcast(array, dimensions, MPI_DOUBLE, ownerofpartition[k], comm);
            if(cg_array_write(ar_names[nud-1][na-1].c_str(), CGNS_ENUMV(RealDouble), 1, &dimensions, array)) cg_error_exit();

          }
        }

      }

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(cgp_close(index_file)) cgp_error_exit();
}
