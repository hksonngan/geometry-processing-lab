#include "prescribedmeancurvature.h"

PrescribedMeanCurvature::PrescribedMeanCurvature()
{
}



void ExplicitAnisotropicMeanCurvature::smooth(int _iterations)
{

    for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS) ; o_it != PluginFunctions::objectsEnd(); ++o_it) {

    bool selectionExists = false;
    double step = 0.00001;
    double lambda = 0.5;
    double r = 10;
    //there should be some singular vertices that make the smooth vector become 0
    //it's because of the mean curvature He

    if ( o_it->dataType( DATA_TRIANGLE_MESH ) ) {

        // Get the mesh to work on
      //TriMesh* mesh = PluginFunctions::triMesh(*o_it);
        TriMeshObject * meshObject = PluginFunctions::triMeshObject(o_it);
        TriMesh* mesh = meshObject->mesh();

      // Property for the active mesh to store original point positions
      OpenMesh::VPropHandleT< TriMesh::Normal > smoothVector;
      OpenMesh::VPropHandleT< double > areaStar;

      // Add a property to the mesh to store mean curvature and area
      mesh->add_property( smoothVector, "explicitAnisotropicMeanCurvature" );
      mesh->add_property( areaStar, "areaStar" );

      mesh->request_vertex_normals();
      mesh->request_vertex_colors();
      mesh->request_face_normals();


      mesh->update_normals();

      for ( int i = 0 ; i < _iterations ; ++i )
      {

          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              mesh->property(smoothVector,v_it).vectorize(0.0f);
              mesh->property(areaStar,v_it) = 0;
              selectionExists |= mesh->status(v_it).selected();
          }





          //last step update all vertices, so the geometry has not changed in the previous calculation
          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {

              //mesh->set_point(v_it, mesh->point(v_it) + updateVector);
              //TriMesh::Normal updateVector;
              //updateVector.vectorize(0);

              for (TriMesh::VertexEdgeIter ve_it=mesh->ve_iter(v_it.handle()); ve_it; ++ve_it)
              {

                  TriMesh::Normal edgeNormal;
                  TriMesh::Scalar meanCurvature = edgeMeanCurvature(mesh, ve_it.handle(), edgeNormal);

                  mesh->property(smoothVector, v_it) += 0.5*meanCurvature*anisotropicWeight(meanCurvature, lambda, r)*edgeNormal;



              }

              for (TriMesh::VertexFaceIter vf_it=mesh->vf_iter(v_it.handle()); vf_it; ++vf_it)
              {

                  mesh->property(areaStar, v_it) += faceArea(mesh, vf_it.handle());
              }


          }

          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              if(selectionExists && mesh->status(v_it).selected() == false) {
                continue;
              }

              TriMesh::Scalar area = mesh->property(areaStar, v_it);
              TriMesh::Normal updateVector = mesh->property(smoothVector, v_it);

              if(selectionExists)
              {
                  printf("area: %f update length: %f\n", area, updateVector.norm());

              }

              mesh->set_point(v_it, mesh->point(v_it) - (3*step/area)*updateVector);
          }

          mesh->update_normals();


      }// Iterations end


      mesh->update_normals();



      updateLineNode(meshObject, smoothVector, areaStar);


      emit updatedObject( o_it->id(), UPDATE_ALL );

      // Create backup
      emit createBackup(o_it->id(), "ExplicitAnisotropicMeanCurvature Smoothing", UPDATE_ALL );

   } else if ( o_it->dataType( DATA_POLY_MESH ) ) {

       // Get the mesh to work on
      PolyMesh* mesh = PluginFunctions::polyMesh(*o_it);

      // Property for the active mesh to store original point positions
      OpenMesh::VPropHandleT< TriMesh::Point > origPositions;

      // Add a property to the mesh to store original vertex positions
      mesh->add_property( origPositions, "SmootherPlugin_Original_Positions" );

      for ( int i = 0 ; i < _iterations ; ++i ) {

         // Copy original positions to backup ( in Vertex property )
         PolyMesh::VertexIter v_it, v_end=mesh->vertices_end();
         for (v_it=mesh->vertices_begin(); v_it!=v_end; ++v_it) {
            mesh->property( origPositions, v_it ) = mesh->point(v_it);
            // See if at least one vertex has been selected
            selectionExists |= mesh->status(v_it).selected();
         }

         // Do one smoothing step (For each point of the mesh ... )
         for (v_it=mesh->vertices_begin(); v_it!=v_end; ++v_it) {

            if(selectionExists && mesh->status(v_it).selected() == false) {
              continue;
            }

            PolyMesh::Point point = PolyMesh::Point(0.0,0.0,0.0);

            // Flag, to skip boundary vertices
            bool skip = false;

            // ( .. for each Outoing halfedge .. )
            PolyMesh::VertexOHalfedgeIter voh_it(*mesh,v_it);
            for ( ; voh_it; ++voh_it ) {
               // .. add the (original) position of the Neighbour ( end of the outgoing halfedge )
               point += mesh->property( origPositions, mesh->to_vertex_handle(voh_it) );

               // Check if the current Halfedge is a boundary halfedge
               // If it is, abort and keep the current vertex position
               if ( mesh->is_boundary( voh_it.handle() ) ) {
                  skip = true;
                  break;
               }

            }

            // Devide by the valence of the current vertex
            point /= mesh->valence( v_it );

            if ( ! skip ) {
               // Set new position for the mesh if its not on the boundary
               mesh->point(v_it) = point;
            }
         }

      }// Iterations end

      // Remove the property
      mesh->remove_property( origPositions );

      mesh->update_normals();

      emit updatedObject( o_it->id() , UPDATE_ALL);

      // Create backup
      emit createBackup(o_it->id(), "ExplicitAnisotropicMeanCurvature Smoothing", UPDATE_ALL);

    } else {

      emit log(LOGERR, "DataType not supported.");
    }
  }

  // Show script logging
  emit scriptInfo("simpleLaplace(" + QString::number(_iterations) + ")");

  emit updateView();
}


double PrescribedMeanCurvature::edgeMeanCurvature(TriMesh *_mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal)
{

    //dihedral = 0 on the boundary
    double dihedral = 0;

    TriMesh::HalfedgeHandle  hh1, hh2;
    TriMesh::FaceHandle    fh1, fh2;

    hh1 = _mesh->halfedge_handle(_eh, 0);
    hh2 = _mesh->halfedge_handle(_eh, 1);

    //normal is the edge vector rotated by 90 degree
    if ( _mesh->is_boundary( hh1 ) )
    {

        fh2 = _mesh->face_handle(hh2);
        TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);
        TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh2)) - _mesh->point(_mesh->from_vertex_handle(hh2));
        normal = (edgeVector%n2);

    }else if ( _mesh->is_boundary( hh2 ) )
    {

        fh1 = _mesh->face_handle(hh1);
        TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
        TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh1)) - _mesh->point(_mesh->from_vertex_handle(hh1));
        normal = (edgeVector%n1);

    }else
    {

        fh1 = _mesh->face_handle(hh1);
        fh2 = _mesh->face_handle(hh2);

        TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
        TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);

        normal = (n1+n2);
        normal.normalize();

        dihedral = _mesh->calc_dihedral_angle(_eh);

        if (dihedral < 0) normal *= -1;

    }

    double edgeLength = _mesh->calc_edge_length(_eh);
    return 2*edgeLength*cos(dihedral/2.0);

}

double PrescribedMeanCurvature::anisotropicWeight(double curvature, double lambda, double r)
{

    double weight;

    //printf("lambda curvature %f %f \n", lambda, curvature);

    if (fabs(curvature) <= lambda)
    {
        weight = 1;
        //printf("weight %f \n", weight);
    }
    else
    {
        weight = (lambda*lambda)/(r*pow(lambda-fabs(curvature), 2) + lambda*lambda);
    }

    //printf("weight %f \n", weight);

    return weight;
}

TriMesh::Scalar PrescribedMeanCurvature::faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh)
{
    // calaculate face area
    TriMesh::Point  p1, p2, p3;
    TriMesh::VertexHandle v1, v2, v3;
    TriMesh::Scalar area;

#define heh halfedge_handle
#define nheh next_halfedge_handle
#define tvh to_vertex_handle
#define fvh from_vertex_handle
    v1 = _mesh->tvh(_mesh->heh(fh));
    v2 = _mesh->fvh(_mesh->heh(fh));
    v3 = _mesh->tvh(_mesh->nheh(_mesh->heh(fh)));
#undef heh
#undef nheh
#undef tvh
#undef fvh

    p1 = _mesh->point(v1);
    p2 = _mesh->point(v2);
    p3 = _mesh->point(v3);

    area = ((p2 - p1) % (p3 - p1)).norm();
    area /= 2.0;

    return area;
}
