#include "prescribedmeancurvature.h"

PrescribedMeanCurvature::PrescribedMeanCurvature()
{
}



void PrescribedMeanCurvature::smooth(int _iterations, TriMeshObject * meshObject)
{


    bool selectionExists = false;

    /*
      bunny: lambda = 0.1 and number of feature vertices: 1389 in total 8810
      cylinder: lambda = 0.32 or 0.31
      cylinder does not manifest the problem with curved edges, but the bunny does
      espectially in the ears regions having curved edges (or high curvature features)
      even with smaller time step
    */
    double lambda = 0.1;


        TriMesh* mesh = meshObject->mesh();

      // Property for the active mesh to store original point positions
      OpenMesh::VPropHandleT< TriMesh::Normal > anisoMeanCurvature;
      //the smoothedAMC is for temporary use
      OpenMesh::VPropHandleT< TriMesh::Normal > smoothedAMC;
      OpenMesh::VPropHandleT< double > areaStar;
      OpenMesh::VPropHandleT< double > anisoPMC;
      OpenMesh::VPropHandleT< double > smoothedAPMC;
      OpenMesh::VPropHandleT< bool > isFeature;
      OpenMesh::VPropHandleT< TriMesh::Normal > volumeGradientProp;

      // Add a property to the mesh to store mean curvature and area
      mesh->add_property( anisoMeanCurvature, "explicitAnisotropicMeanCurvature" );
      mesh->add_property( smoothedAMC, "smoothedAnisotropicMeanCurvature" );
      mesh->add_property( areaStar, "areaStar" );
      mesh->add_property( anisoPMC, "anisotropicPrescribedMeanCurvature" );
      mesh->add_property( smoothedAPMC, "smoothedAPMC" );
      mesh->add_property( isFeature, "isFeature" );
      mesh->add_property( volumeGradientProp, "volumeGradientProp" );

      mesh->request_vertex_normals();
      mesh->request_vertex_colors();
      mesh->request_face_normals();

      unsigned int count(meshObject->mesh()->n_vertices());


      mesh->update_normals();

      for ( int i = 0 ; i < _iterations ; ++i )
      {

          unsigned int noFeatureVertices = 0;

          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              mesh->property(anisoMeanCurvature,v_it).vectorize(0.0f);
              mesh->property(smoothedAMC,v_it).vectorize(0.0f);
              mesh->property(volumeGradientProp,v_it).vectorize(0.0f);
              mesh->property(areaStar,v_it) = 0;
              mesh->property(anisoPMC,v_it) = 0;
              mesh->property(smoothedAPMC,v_it) = 0;
              mesh->property(isFeature,v_it) = false;
              selectionExists |= mesh->status(v_it).selected();
          }


          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              TriMesh::Normal isotropic;
              isotropic.vectorize(0);

              for (TriMesh::VertexEdgeIter ve_it=mesh->ve_iter(v_it.handle()); ve_it; ++ve_it)
              {

                  TriMesh::Normal edgeNormal;
                  TriMesh::Scalar meanCurvature = edgeMeanCurvature(mesh, ve_it.handle(), edgeNormal);
                  double weight = anisotropicWeight(meanCurvature, lambda, R);

                  mesh->property(anisoMeanCurvature, v_it) += 0.5*meanCurvature*weight*edgeNormal;

                  isotropic += 0.5*meanCurvature*edgeNormal;

              }

              if (mesh->property(anisoMeanCurvature, v_it) != isotropic)
              {
                  mesh->property(isFeature,v_it) = true;
                  noFeatureVertices++;
                  //TriMesh::Normal aniso = mesh->property(smoothVector, v_it);
                  //if (noFeatureVertices%10 == 0) printf("aniso x %f iso x %f y %f %f z %f %f \n",
                  //       aniso[0], isotropic[0], aniso[1], isotropic[1], aniso[2], isotropic[2]);

              }

              for (TriMesh::VertexFaceIter vf_it=mesh->vf_iter(v_it.handle()); vf_it; ++vf_it)
              {

                  mesh->property(areaStar, v_it) += faceArea(mesh, vf_it.handle());
                  //volume gradient
                  TriMesh::Normal volGrad;
                  volumeGradient(mesh, vf_it.handle(), v_it, volGrad);
                  mesh->property(volumeGradientProp, v_it) += volGrad;
              }

          }

          smoothAnisotropicMeanCurvature(mesh, anisoMeanCurvature, smoothedAMC, volumeGradientProp, isFeature);

          computePMCFunction(mesh, anisoMeanCurvature, volumeGradientProp, isFeature, anisoPMC, smoothedAPMC);

          printf("number of feature vertices: %d in total %d \n", noFeatureVertices, count);

          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              if(selectionExists && mesh->status(v_it).selected() == false) {
                continue;
              }

              TriMesh::Scalar area = mesh->property(areaStar, v_it);
              TriMesh::Normal updateVector = mesh->property(anisoMeanCurvature, v_it);
              updateVector -= mesh->property(smoothedAPMC, v_it)*mesh->property(volumeGradientProp, v_it);

              mesh->set_point(v_it, mesh->point(v_it) - (3*TIME_STEP/area)*updateVector);
          }

          mesh->update_normals();


      }// Iterations end


      mesh->update_normals();
      updateLineNode(meshObject, anisoMeanCurvature, areaStar);


      // Remove the property
      mesh->remove_property( anisoMeanCurvature );
      mesh->remove_property( areaStar );
      mesh->remove_property( isFeature );
      mesh->remove_property( volumeGradientProp );
      mesh->remove_property( smoothedAMC );
      mesh->remove_property( smoothedAPMC );
      mesh->remove_property( anisoPMC );

}


void PrescribedMeanCurvature::
computePMCFunction(TriMesh *_mesh
                   , const OpenMesh::VPropHandleT< TriMesh::Normal > & anisoMeanCurvature
                   , const OpenMesh::VPropHandleT< TriMesh::Normal > & volumeGrad
                   , const OpenMesh::VPropHandleT< bool > & isFeature
                   , const OpenMesh::VPropHandleT< double > & anisoPMC
                   , const OpenMesh::VPropHandleT< double > & smoothedAPMC)
{
    for (TriMesh::VertexIter v_it=_mesh->vertices_begin(); v_it!=_mesh->vertices_end(); ++v_it)
    {
        TriMesh::Normal volGrad = _mesh->property(volumeGrad, v_it);
        TriMesh::Normal anisoMCVec = _mesh->property(anisoMeanCurvature, v_it);
        TriMesh::Scalar pmc = anisoMCVec.norm()/volGrad.norm();

        _mesh->property(anisoPMC, v_it) = pmc;

    }

    //smooth by averaging over the neighborhood
    for (TriMesh::VertexIter v_it=_mesh->vertices_begin(); v_it!=_mesh->vertices_end(); ++v_it)
    {

        bool feature = _mesh->property(isFeature, v_it);
        TriMesh::Scalar avg = 0;
        int i = 0;

        TriMesh::VertexVertexIter vv_iter(*_mesh,v_it);
        for ( ; vv_iter; ++vv_iter )
        {
            //if the center vertex is feature vertex we only average over feature neighbors
            if (feature && !_mesh->property(isFeature, vv_iter)) continue;

            avg += _mesh->property(anisoPMC, vv_iter);
            i++;

        }

        avg /= i;
        _mesh->property(smoothedAPMC, v_it) = avg;

    }

}


void PrescribedMeanCurvature::
smoothAnisotropicMeanCurvature(TriMesh *_mesh
                               , const OpenMesh::VPropHandleT< TriMesh::Normal > & anisoMeanCurvature
                               , const OpenMesh::VPropHandleT< TriMesh::Normal > & smoothedAMC
                               , const OpenMesh::VPropHandleT< TriMesh::Normal > & volumeGrad
                               , const OpenMesh::VPropHandleT< bool > & isFeature)
{
    for (TriMesh::VertexIter v_it=_mesh->vertices_begin(); v_it!=_mesh->vertices_end(); ++v_it)
    {

        if (!_mesh->property(isFeature, v_it)) continue;


        TriMesh::Point p, q, r1, r2;
        p = _mesh->point(v_it);
        double totalAngle = 0;
        TriMesh::Normal unitVec;
        TriMesh::Normal volGradient;


        // ( .. for each Incoming halfedge .. )
        TriMesh::VertexIHalfedgeIter vih_it(*_mesh,v_it);
        for ( ; vih_it; ++vih_it )
        {

            double sumAngle = 0;
            TriMesh::VertexHandle qvh = _mesh->from_vertex_handle(vih_it);
            q = _mesh->point(qvh);

            if (!_mesh->is_boundary(vih_it))
            {

                r1 = _mesh->point(_mesh->to_vertex_handle(_mesh->next_halfedge_handle(vih_it)));
                sumAngle += calAngle(p, q, r1);

            }

            TriMesh::HalfedgeHandle opp = _mesh->opposite_halfedge_handle(vih_it);

            if (!_mesh->is_boundary(opp))
            {
                r2 = _mesh->point(_mesh->from_vertex_handle(_mesh->prev_halfedge_handle(opp)));
                sumAngle += calAngle(p, q, r2);

            }

            totalAngle += sumAngle;

            _mesh->property(smoothedAMC, v_it) += sumAngle*_mesh->property(anisoMeanCurvature, qvh);

        }

        unitVec = 0.5*(_mesh->property(anisoMeanCurvature, v_it) + _mesh->property(smoothedAMC, v_it)/totalAngle);
        unitVec.normalize();

        volGradient = _mesh->property(volumeGrad, v_it);
        double sign = volGradient|unitVec;
        if (sign < 0) unitVec *= -1;

        _mesh->property(volumeGrad, v_it) = unitVec;

    }
}





double PrescribedMeanCurvature::
calAngle(const TriMesh::Point & p, const TriMesh::Point & q, const TriMesh::Point & r)
{
    TriMesh::Normal pq = q - p;
    TriMesh::Normal pr = r - p;
    pq.normalize();
    pr.normalize();

    return acos((pq|pr));
}


double PrescribedMeanCurvature::
edgeMeanCurvature(TriMesh *_mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal)
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

    //printf("lambda %f curvature %f \n", lambda, curvature);

    if (fabs(curvature) <= lambda)
    {
        weight = 1;
    }
    else
    {
        weight = (lambda*lambda)/(r*pow(lambda-fabs(curvature), 2) + lambda*lambda);
    }

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


/*
 * circulate around the oriented face using FaceHalfEdgeIter to extract p, q, r in CCW
 *
 * take (qxr)/6
 */
void
PrescribedMeanCurvature::
volumeGradient(TriMesh *_mesh, TriMesh::FaceHandle fh, TriMesh::VertexHandle vh, TriMesh::Normal & gradient)
{

    TriMesh::Point  q, r;
    TriMesh::VertexHandle vq, vr;

    //from p to q to r
    for (TriMesh::FaceHalfedgeIter fhh_it=_mesh->fh_iter(fh); fhh_it; ++fhh_it)
    {
        if (_mesh->from_vertex_handle(fhh_it) == vh)
        {
            vq = _mesh->to_vertex_handle(fhh_it);
            vr = _mesh->to_vertex_handle(_mesh->next_halfedge_handle(fhh_it));
        }
    }


    gradient = q%r;
    gradient /= 6;

}


void
PrescribedMeanCurvature::
updateLineNode(TriMeshObject * _meshObject, OpenMesh::VPropHandleT< TriMesh::Normal > & anisoMeanCurvature, OpenMesh::VPropHandleT< double >& areaStar)
{
  ACG::SceneGraph::LineNode * node = getLineNode(_meshObject);
  //OpenMesh::VPropHandleT< TriMesh::Normal > smoothVector;

  node->clear();

  for (TriMesh::VertexIter vit = _meshObject->mesh()->vertices_begin();
                          vit != _meshObject->mesh()->vertices_end(); ++vit)
  {
    TriMesh::Point  p = _meshObject->mesh()->point(vit);
    TriMesh::Normal n = _meshObject->mesh()->property(anisoMeanCurvature, vit);
    TriMesh::Scalar coefficient = _meshObject->mesh()->property(areaStar, vit);
    coefficient = 3/coefficient;//no time step involved
    addLine(node, p, p+coefficient*n, Color(255,0,0) );
  }
}

ACG::SceneGraph::LineNode *
PrescribedMeanCurvature::
getLineNode(TriMeshObject * _meshObject)
{
  ACG::SceneGraph::LineNode * line_node = 0;

  // get or add line node
  if( !_meshObject->hasAdditionalNode( "NormalEstimationPlugin", "LineNode" ) )
  {
    line_node = new ACG::SceneGraph::LineNode( ACG::SceneGraph::LineNode::LineSegmentsMode, _meshObject->manipulatorNode() );

    if( !_meshObject->addAdditionalNode(line_node, QString("NormalEstimationPlugin"), QString("LineNode") ) )
    {
      std::cerr << "NormalEstimationPlugin::getLineNode(): could not add line node.\n";
      return 0;
    }
    line_node->clear();
    line_node->set_line_width(1.0);
    line_node->set_base_color( ACG::Vec4f(255,0,255,255) );
    line_node->show();
  }
  else
    _meshObject->getAdditionalNode ( line_node, "NormalEstimationPlugin", "LineNode" );

  return line_node;
}

void
PrescribedMeanCurvature::
addLine( ACG::SceneGraph::LineNode * _line_node, Vec3d _p0, Vec3d _p1, Color _color )
{
  _line_node->add_line( _p0, _p1 );
  _line_node->add_color(_color);
}
