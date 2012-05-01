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
    double step = 0.001;
    double lambda = 0.1;
    double r = 10;

    if ( o_it->dataType( DATA_TRIANGLE_MESH ) ) {

        // Get the mesh to work on
      TriMesh* mesh = PluginFunctions::triMesh(*o_it);

      // Property for the active mesh to store original point positions
      OpenMesh::VPropHandleT< TriMesh::Normal > smoothVector;
      OpenMesh::VPropHandleT< double > areaStar;

      // Add a property to the mesh to store mean curvature and area
      mesh->add_property( smoothVector, "explicitAnisotropicMeanCurvature" );
      mesh->add_property( areaStar, "areaStar" );

      mesh->update_normals();

      for ( int i = 0 ; i < _iterations ; ++i )
      {

          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              mesh->property(smoothVector,v_it).vectorize(0.0f);
              mesh->property(areaStar,v_it) = 0;
          }

          for (TriMesh::EdgeIter e_it=mesh->edges_begin(); e_it!=mesh->edges_end(); ++e_it)
          // do something with *e_it, e_it->, or e_it.handle()
          {

              TriMesh::HalfedgeHandle  hh;
              TriMesh::VertexHandle    v0, v1;

              hh = mesh->halfedge_handle(e_it.handle(), 0);
              v0 = mesh->to_vertex_handle(hh);
              hh = mesh->halfedge_handle(e_it.handle(), 1);
              v1 = mesh->to_vertex_handle(hh);


              TriMesh::Normal edgeNormal;
              TriMesh::Normal smoothDirection;

              edgeNormal.vectorize(0.0);

              //printf("normal before %f %f %f \n", edgeNormal[0], edgeNormal[1], edgeNormal[2]);
              double meanCurvature = edgeMeanCurvature(mesh, e_it.handle(), edgeNormal);
             //printf("normal after %f %f %f \n", edgeNormal[0], edgeNormal[1], edgeNormal[2]);

              smoothDirection = (meanCurvature*anisotropicWeight(meanCurvature, lambda, r))*edgeNormal;
              mesh->property(smoothVector, v0) -= smoothDirection;
              mesh->property(smoothVector, v1) -= smoothDirection;

          }

          for (TriMesh::FaceIter f_it=mesh->faces_begin(); f_it!=mesh->faces_end(); ++f_it)
          {
              TriMesh::Scalar area = faceArea(mesh, f_it.handle(), areaStar);




          }

          //last step update all vertices, so the geometry has not changed in the previous calculation
          for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
          {
              TriMesh::Scalar coefficient = 3*step/(2*mesh->property(areaStar, v_it.handle()));
              TriMesh::Normal updateVector = coefficient*mesh->property(smoothVector, v_it.handle());
              mesh->set_point(v_it, mesh->point(v_it) + updateVector);

          }

          mesh->update_normals();


      }// Iterations end


      mesh->update_normals();

      emit updatedObject( o_it->id(), UPDATE_GEOMETRY );
      
      // Create backup
      emit createBackup(o_it->id(), "ExplicitAnisotropicMeanCurvature Smoothing", UPDATE_GEOMETRY );

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

      emit updatedObject( o_it->id() , UPDATE_GEOMETRY);
      
      // Create backup
      emit createBackup(o_it->id(), "ExplicitAnisotropicMeanCurvature Smoothing", UPDATE_GEOMETRY);

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

    double dihedral = _mesh->calc_dihedral_angle_fast(_eh);
    //double dihedral;

    TriMesh::HalfedgeHandle  hh1, hh2;
    TriMesh::FaceHandle    fh1, fh2;

    hh1 = _mesh->halfedge_handle(_eh, 0);
    hh2 = _mesh->halfedge_handle(_eh, 1);

    fh1 = _mesh->face_handle(hh1);
    fh2 = _mesh->face_handle(hh2);

    TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
    TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);

    normal = n1+n2;
    normal /= normal.norm();//or normal.normalize();

    double edgeLength = _mesh->calc_edge_length(_eh);

    #define PI 3.14159265

    //dihedral = PI - acos(n1|n2);

    //printf("dihedral %f \n", dihedral*180/PI);

    return 2*edgeLength*cos(dihedral/2.0);

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

TriMesh::Scalar ExplicitAnisotropicMeanCurvature::faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh, const OpenMesh::VPropHandleT< double > & areaStar)
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

    _mesh->property(areaStar, v1) += area;
    _mesh->property(areaStar, v2) += area;
    _mesh->property(areaStar, v3) += area;

    return area;
}


Q_EXPORT_PLUGIN2( explicitAnisotropicMeanCurvature , ExplicitAnisotropicMeanCurvature );

