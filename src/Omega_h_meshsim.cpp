#include "Omega_h_file.hpp"

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"

namespace Omega_h {

namespace meshsim {

namespace {

int classId(pEntity e) {
  pGEntity g = EN_whatIn(e);
  assert(g);
  return GEN_tag(g);
}

void read_internal(pParMesh sm, Mesh* mesh) {
  (void)mesh;
  pMesh m = PM_mesh(sm, 0);
  Int max_dim;
  if (M_numRegions(m)) {
    max_dim = 3;
  } else if (M_numFaces(m)) {
    max_dim = 2;
  } else if (M_numEdges(m)) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }

  //get the types of elements
  Omega_h_Family family = OMEGA_H_SIMPLEX;
  RIter regions = M_regionIter(m);
  pRegion rgn;
  LO i = 0;
  while ((rgn = (pRegion) RIter_next(regions))) {
    if(R_topoType(rgn) != Rtet)
      Omega_h_fail("Non-simplex element found!\n");
    ++i;
  }
  RIter_delete(regions);
  std::vector<int> ent_nodes[4];
  std::vector<int> ent_class_ids[4];
  std::vector<int> ent_matches[3];
  std::vector<int> ent_match_classId[3];
  //write vertex coords into node_coords and
  //  write vertex ids into ent_nodes
  const int numVtx = M_numVertices(m);
  ent_nodes[0].reserve(numVtx);
  ent_class_ids[0].reserve(numVtx);
  ent_matches[0].reserve(numVtx);
  ent_match_classId[0].reserve(numVtx);
  HostWrite<Real> host_coords(numVtx*max_dim);
  VIter vertices = M_vertexIter(m);
  pVertex vtx;
  i = 0;
  int count_matches = 0;
  int count_matched = 0;
  while ((vtx = (pVertex) VIter_next(vertices))) {
    double xyz[3];
    V_coord(vtx,xyz);
    if( max_dim < 3 && xyz[2] != 0 )
      Omega_h_fail("The z coordinate must be zero for a 2d mesh!\n");
    for(int j=0; j<max_dim; j++) {
      host_coords[i * max_dim + j] = xyz[j];
    }
    ent_nodes[0].push_back(EN_id(vtx));
    ent_class_ids[0].push_back(classId(vtx));
    ++i;

    pPList matches = EN_getMatchingEnts(vtx, 0, 0);
    void *iterM = 0;
    pVertex match;
    count_matches = 0;
    while((match = (pVertex)PList_next(matches, &iterM))) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(vtx))) {
      //if (EN_id(match) != EN_id(vtx)) {
        printf("vtx %d is a match to original vtx %d\n", EN_id(match), EN_id(vtx));
        ent_matches[0].push_back(EN_id(match));
        ent_match_classId[0].push_back(classId(match));
        ++count_matches;
        ++count_matched;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");
      }
      else if (PList_size(matches)==1) {
      //else if ((Plist_size(matches)==1) && (EN_id(match) == EN_id(vtx))) {
        printf("vtx %d has no match\n", EN_id(vtx));
        ent_matches[0].push_back(-1);
        ent_match_classId[0].push_back(-1);
        //this '-1' has been put in for tag visualization.might not need if
        //we use CSR
      }
    }
    PList_delete(matches);
    //incase of count_matches>1 will need to comeup with different data
    //structure to store matches of variable size. One idea is to use CSR
  }
  VIter_delete(vertices);
  printf("matched verts=%d \n", count_matched);
  //get the ids of vertices bounding each edge
  const int numEdges = M_numEdges(m);
  ent_nodes[1].reserve(numEdges*2);
  ent_class_ids[1].reserve(numEdges);
  EIter edges = M_edgeIter(m);
  pEdge edge;
  while ((edge = (pEdge) EIter_next(edges))) {
    for(int j=1; j>=0; --j) {
      vtx = E_vertex(edge,j);
      ent_nodes[1].push_back(EN_id(vtx));
    }
    ent_class_ids[1].push_back(classId(edge));

    pPList matches = EN_getMatchingEnts(edge, 0, 0);
    void *iterM = 0;
    pEdge match;
    count_matches = 0;
    while((match = (pEdge)PList_next(matches, &iterM))) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(edge))) {
      //if (EN_id(match) != EN_id(edge)) {
        ent_matches[1].push_back(EN_id(match));
        ent_match_classId[1].push_back(classId(match));
        ++count_matches;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");
      }
      else if (PList_size(matches)==1) {
      //else {
        ent_matches[1].push_back(-1);
        ent_match_classId[1].push_back(-1);
      }
    }
    PList_delete(matches);
  }
  EIter_delete(edges);
  //get the ids of vertices bounding each face
  const int numFaces = M_numFaces(m);
  ent_nodes[2].reserve(numFaces*3);
  ent_class_ids[2].reserve(numFaces);
  FIter faces = M_faceIter(m);
  pFace face;
  while ((face = (pFace) FIter_next(faces))) {
    pPList verts = F_vertices(face,1);
    assert(PList_size(verts) == 3);
    void *iter = 0; // must initialize to 0
    while((vtx = (pVertex)PList_next(verts, &iter)))
      ent_nodes[2].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[2].push_back(classId(face));

    pPList matches = EN_getMatchingEnts(face, 0, 0);
    void *iterM = 0;
    pFace match;
    count_matches = 0;
    while((match = (pFace)PList_next(matches, &iterM))) {
      if ((PList_size(matches)>1) && (EN_id(match) != EN_id(face))) {
      //if (EN_id(match) != EN_id(face)) {
        ent_matches[2].push_back(EN_id(match));
        ent_match_classId[2].push_back(classId(match));
        ++count_matches;
        if (count_matches > 1) Omega_h_fail("Error:matches per entity > 1\n");
      }
      else if (PList_size(matches)==1) {
      //else {
        ent_matches[2].push_back(-1);
        ent_match_classId[2].push_back(-1);
      }
    }
    PList_delete(matches);
  }
  FIter_delete(faces);
  //get the ids of vertices bounding each region
  const int numRegions = M_numRegions(m);
  ent_nodes[3].reserve(numRegions*4);
  ent_class_ids[3].reserve(numRegions);
  regions = M_regionIter(m);
  while ((rgn = (pRegion) RIter_next(regions))) {
    pPList verts = R_vertices(rgn,1);
    assert(PList_size(verts) == 4);
    void *iter = 0; // must initialize to 0
    while((vtx = (pVertex)PList_next(verts, &iter)))
      ent_nodes[3].push_back(EN_id(vtx));
    PList_delete(verts);
    ent_class_ids[3].push_back(classId(rgn));
  }
  RIter_delete(regions);
  //flatten the ent_nodes and ent_class_ids arrays
  for (Int ent_dim = max_dim; ent_dim >= 0; --ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size()) / neev;
    HostWrite<LO> host_ev2v(ndim_ents * neev);
    HostWrite<LO> host_class_id(ndim_ents);
    for (i = 0; i < ndim_ents; ++i) {
      for (Int j = 0; j < neev; ++j) {
        host_ev2v[i * neev + j] =
            ent_nodes[ent_dim][static_cast<std::size_t>(i * neev + j)];
      }
      host_class_id[i] = ent_class_ids[ent_dim][static_cast<std::size_t>(i)];
    }
    auto eqv2v = Read<LO>(host_ev2v.write());
    if (ent_dim == max_dim) {
      build_from_elems_and_coords(
          mesh, family, max_dim, eqv2v, host_coords.write());
    }
    classify_equal_order(mesh, ent_dim, eqv2v, host_class_id.write());
  }
  finalize_classification(mesh);

  //add matches info to mesh
  for (Int ent_dim = 0; ent_dim < 3; ++ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size())/neev;
    HostWrite<LO> host_matches(ndim_ents);
    HostWrite<LO> host_match_classId(ndim_ents);
    for (i = 0; i < ndim_ents; ++i) {
      host_matches[i] = ent_matches[ent_dim][static_cast<std::size_t>(i)];
      host_match_classId[i] = ent_match_classId[ent_dim][static_cast<std::size_t>(i)];
      //how to store when multiple matches on different parts in parallel?
    }
    auto matches = Read<LO>(host_matches.write());
    auto match_classId = Read<LO>(host_match_classId.write());
    mesh->add_tag(ent_dim, "matches", 1, matches);
    mesh->add_tag(ent_dim, "match_classId", 1, match_classId);
  }
  //
}

}  // end anonymous namespace

Mesh read(filesystem::path const& mesh_fname, filesystem::path const& mdl_fname,
    CommPtr comm) {
  SimPartitionedMesh_start(NULL,NULL);
  SimModel_start();
  Sim_readLicenseFile(NULL);
  pNativeModel nm = NULL;
  pProgress p = NULL;
    printf("ok0\n");
  pGModel g = GM_load(mdl_fname.c_str(), nm, p);
    printf("ok1\n");
  pParMesh sm = PM_load(mesh_fname.c_str(), g, p);
    printf("ok2\n");
  auto mesh = Mesh(comm->library());
  meshsim::read_internal(sm, &mesh);
  M_release(sm);
  GM_release(g);
  SimModel_stop();
  SimPartitionedMesh_stop();
  return mesh;
}

}  // namespace meshsim

}  // end namespace Omega_h
