/*===========================================================================*\
*                                                                            *
*                              OpenFlipper                                   *
*      Copyright (C) 2001-2011 by Computer Graphics Group, RWTH Aachen       *
*                           www.openflipper.org                              *
*                                                                            *
*--------------------------------------------------------------------------- *
*  This file is part of OpenFlipper.                                         *
*                                                                            *
*  OpenFlipper is free software: you can redistribute it and/or modify       *
*  it under the terms of the GNU Lesser General Public License as            *
*  published by the Free Software Foundation, either version 3 of            *
*  the License, or (at your option) any later version with the               *
*  following exceptions:                                                     *
*                                                                            *
*  If other files instantiate templates or use macros                        *
*  or inline functions from this file, or you compile this file and          *
*  link it with other files to produce an executable, this file does         *
*  not by itself cause the resulting executable to be covered by the         *
*  GNU Lesser General Public License. This exception does not however        *
*  invalidate any other reasons why the executable file might be             *
*  covered by the GNU Lesser General Public License.                         *
*                                                                            *
*  OpenFlipper is distributed in the hope that it will be useful,            *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*  GNU Lesser General Public License for more details.                       *
*                                                                            *
*  You should have received a copy of the GNU LesserGeneral Public           *
*  License along with OpenFlipper. If not,                                   *
*  see <http://www.gnu.org/licenses/>.                                       *
*                                                                            *
\*===========================================================================*/

/*===========================================================================*\
*                                                                            *
*   $Revision: 14436 $                                                       *
*   $LastChangedBy: phan-anh $                                               *
*   $Date: 2012-04-23 15:54:08 +0200 (Mon, 23 Apr 2012) $                    *
*                                                                            *
\*===========================================================================*/



#include "ExplicitAnisotropicMeanCurvature.hh"
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"


ExplicitAnisotropicMeanCurvature::ExplicitAnisotropicMeanCurvature() :
        iterationsSpinbox_(0)
{

}

ExplicitAnisotropicMeanCurvature::~ExplicitAnisotropicMeanCurvature()
{

}

void ExplicitAnisotropicMeanCurvature::initializePlugin()
{
   // Create the Toolbox Widget
   QWidget* toolBox = new QWidget();
   QGridLayout* layout = new QGridLayout(toolBox);

   QPushButton* smoothButton = new QPushButton("&Smooth",toolBox);
   smoothButton->setToolTip(tr("Smooths an Object using ExplicitAnisotropicMeanCurvature Smoothing."));
   smoothButton->setWhatsThis(tr("Smooths an Object using ExplicitAnisotropicMeanCurvature Smoothing. Use the Smooth Plugin for more options."));



   iterationsSpinbox_ =  new QSpinBox(toolBox) ;
   iterationsSpinbox_->setMinimum(1);
   iterationsSpinbox_->setMaximum(1000);
   iterationsSpinbox_->setSingleStep(1);
   iterationsSpinbox_->setToolTip(tr("The number of the smooting operations."));
   iterationsSpinbox_->setWhatsThis(tr("Give the number, how often the ExplicitAnisotropicMeanCurvature Smoothing should modify the object."));

   QLabel* label = new QLabel("Iterations:");

   layout->addWidget( label             , 0, 0);
   layout->addWidget( smoothButton      , 1, 1);
   layout->addWidget( iterationsSpinbox_, 0, 1);

   layout->addItem(new QSpacerItem(10,10,QSizePolicy::Expanding,QSizePolicy::Expanding),2,0,1,2);

   connect( smoothButton, SIGNAL(clicked()), this, SLOT(smooth()) );

   QIcon* toolIcon = new QIcon(OpenFlipper::Options::iconDirStr()+OpenFlipper::Options::dirSeparator()+"smoother1.png");
   emit addToolbox( tr("ExplicitAnisotropicMeanCurvature Smoother") , toolBox, toolIcon );
}

void ExplicitAnisotropicMeanCurvature::pluginsInitialized() {
    
    // Emit slot description
    emit setSlotDescription(tr("smooth(int)"),   tr("Smooth mesh using the ExplicitAnisotropicMeanCurvature."),
                            QStringList(tr("iterations")), QStringList(tr("Number of iterations")));
}

/** \brief simpleLaplace
 *
 *  Smooth mesh using the Laplace operator
 *  with uniform weights.
 */
void ExplicitAnisotropicMeanCurvature::smooth() {

    int iterations = 1;
    
    if(!OpenFlipper::Options::nogui()) {
        iterations = iterationsSpinbox_->value();
    }
    
    smooth(iterations);
}

/** \brief simpleLaplace
 *
 * Smooth mesh using the Laplace operator
 * with uniform weights.
 *
 * @param _iterations Number of iterations
 */
void ExplicitAnisotropicMeanCurvature::smooth(int _iterations) {
    
    for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS) ; o_it != PluginFunctions::objectsEnd(); ++o_it) {

    bool selectionExists = false;
    double step = 0.0001;
    double lambda = 0.3;
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


double ExplicitAnisotropicMeanCurvature::edgeMeanCurvature(TriMesh *_mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal)
{


    double dihedral;

    TriMesh::HalfedgeHandle  hh1, hh2;
    TriMesh::FaceHandle    fh1, fh2;

    hh1 = _mesh->halfedge_handle(_eh, 0);
    hh2 = _mesh->halfedge_handle(_eh, 1);

    fh1 = _mesh->face_handle(hh1);
    fh2 = _mesh->face_handle(hh2);

    TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
    TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);

    normal = (n1+n2);
    //normal /= normal.norm();//or normal.normalize();
    normal.normalize();

    double edgeLength = _mesh->calc_edge_length(_eh);

    #define PI 3.14159265

    //dihedral = PI - acos(n1|n2);
    //dihedral = acos(n1|n2);
    dihedral = _mesh->calc_dihedral_angle(_eh);

    //printf("dihedralMesh %f dihedralNorm %f \n", dihedral*180/PI, acos(n1|n2)*180/PI);

    if (dihedral < 0) normal *= -1;

    return 2*edgeLength*cos(dihedral/2.0);
    //return 2*edgeLength*fabs(sin(dihedral/2.0));

}

double ExplicitAnisotropicMeanCurvature::anisotropicWeight(double curvature, double lambda, double r)
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

TriMesh::Scalar ExplicitAnisotropicMeanCurvature::faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh)
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


void
ExplicitAnisotropicMeanCurvature::
updateLineNode(TriMeshObject * _meshObject, OpenMesh::VPropHandleT< TriMesh::Normal > & smoothVector, OpenMesh::VPropHandleT< double >& areaStar)
{
  ACG::SceneGraph::LineNode * node = getLineNode(_meshObject);
  //OpenMesh::VPropHandleT< TriMesh::Normal > smoothVector;

  node->clear();

  for (TriMesh::VertexIter vit = _meshObject->mesh()->vertices_begin();
                          vit != _meshObject->mesh()->vertices_end(); ++vit)
  {
    TriMesh::Point  p = _meshObject->mesh()->point(vit);
    TriMesh::Normal n = _meshObject->mesh()->property(smoothVector, vit);
    TriMesh::Scalar length = _meshObject->mesh()->property(areaStar, vit);
    addLine(node, p, p+length*50*n, Color(255,0,0) );
  }
}

ACG::SceneGraph::LineNode *
ExplicitAnisotropicMeanCurvature::
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
ExplicitAnisotropicMeanCurvature::
addLine( ACG::SceneGraph::LineNode * _line_node, Vec3d _p0, Vec3d _p1, Color _color )
{
  _line_node->add_line( _p0, _p1 );
  _line_node->add_color(_color);
}


Q_EXPORT_PLUGIN2( explicitAnisotropicMeanCurvature , ExplicitAnisotropicMeanCurvature );

