#include "Omega_h_beziers.hpp"
#include "Omega_h_vector.hpp"
#include "Omega_h_scalar.hpp"

LO checkValidity_2d(Mesh *mesh, LOs new_tris) {

  auto fv2v = mesh->ask_down(2, 0).ab2b;
  auto fe2e = mesh->get_adj(2, 1).ab2b;
  auto vertCtrlPts = mesh->get_ctrlPts(0);
  auto edgeCtrlPts = mesh->get_ctrlPts(1);

  Write<I8> is_valid(new_tris.size(), 1);

  auto check_validity = OMEGA_H_LAMBDA (LO i) {
    auto tri = new_tris(i);

    
  };
  parallel_for(new_tris.size(), std::move(check_validity));

  apf::Element* elem = apf::createElement(mesh->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> elemNodes;
  apf::getVectorNodes(elem,elemNodes);

  apf::NewArray<double> nodes(order*(2*order-1));
  // have to use this function because its for x-y plane, and
  // the other method used in 3D does not work in those cases
  getTriJacDetNodes(order,elemNodes,nodes);

  // check vertices
  apf::Downward verts;
  mesh->getDownward(e,0,verts);
  for (int i = 0; i < 3; ++i){
    if(nodes[i] < minAcceptable){
      return i+2;
    }
  }

  apf::MeshEntity* edges[3];
  mesh->getDownward(e,1,edges);
  double minJ = 0, maxJ = 0;
  // Vertices will already be flagged in the first check
  for (int edge = 0; edge < 3; ++edge){
    for (int i = 0; i < 2*(order-1)-1; ++i){
      if (nodes[3+edge*(2*(order-1)-1)+i] < minAcceptable){
        minJ = -1e10;
        apf::NewArray<double> edgeNodes(2*(order-1)+1);
        if(algorithm < 2){
          edgeNodes[0] = nodes[apf::tri_edge_verts[edge][0]];
          edgeNodes[2*(order-1)] = nodes[apf::tri_edge_verts[edge][1]];
          for (int j = 0; j < 2*(order-1)-1; ++j)
            edgeNodes[j+1] = nodes[3+edge*(2*(order-1)-1)+j];
          /*
          if(algorithm == 1){
            getJacDetByElevation(apf::Mesh::EDGE,2*(order-1),edgeNodes,minJ,maxJ);
          } else {
            // allows recursion stop on first "conclusive" invalidity
            bool done = false;
            getJacDetBySubdivision(apf::Mesh::EDGE,2*(order-1),
                0,edgeNodes,minJ,maxJ,done);
          }
          */
        } else {
          edgeNodes[0] = nodes[apf::tri_edge_verts[edge][0]];
          edgeNodes[1] = nodes[apf::tri_edge_verts[edge][1]];
          for (int j = 0; j < 2*(order-1)-1; ++j)
            edgeNodes[j+2] = nodes[3+edge*(2*(order-1)-1)+j];
          bool done = false;
          bool quality = false;
          getJacDetBySubdivisionMatrices(apf::Mesh::EDGE,2*(order-1),
              0,subdivisionCoeffs[1],edgeNodes,minJ,maxJ,done,quality);
        }
        if(minJ < minAcceptable){
          return 8+edge;
        }
      }
    }
  }

  bool done = false;
  for (int i = 0; i < (2*order-3)*(2*order-4)/2; ++i){
    if (nodes[6*(order-1)+i] < minAcceptable){
      minJ = -1e10;
      /*
      if(algorithm == 1)
        getJacDetByElevation(apf::Mesh::TRIANGLE,2*(order-1),nodes,minJ,maxJ);
      else if(algorithm == 2){
        bool quality = false;
        getJacDetBySubdivisionMatrices(apf::Mesh::TRIANGLE,2*(order-1),
            0,subdivisionCoeffs[2],nodes,minJ,maxJ,done,quality);
      } else {
        getJacDetBySubdivision(apf::Mesh::TRIANGLE,2*(order-1),
            0,nodes,minJ,maxJ,done);
      }
      */
      if(minJ < minAcceptable){
        return 14;
      }
    }
  }
  return 1;
}
