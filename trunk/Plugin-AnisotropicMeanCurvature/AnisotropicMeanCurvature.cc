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



#include "AnisotropicMeanCurvature.hh"
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"
#include "ACG/Utils/ColorCoder.hh"


AnisotropicMeanCurvature::AnisotropicMeanCurvature()
{
    smooth_type = PrescribedMeanCurvature::ANISO_MEAN_CURVATURE;
    scheme = PrescribedMeanCurvature::EXPLICIT;
    visualize = PrescribedMeanCurvature::UPDATE_VECTOR;
    pmc.set_lambda(0.8);
}

AnisotropicMeanCurvature::~AnisotropicMeanCurvature()
{

}

void AnisotropicMeanCurvature::initializePlugin()
{

    gui_ = new SmootherToolboxWidget();
    QSize size(100, 100);
    gui_->resize(size);

    //gui_->lineEdit_lambda.

    connect(gui_->pushButton_smooth, SIGNAL(clicked()), this, SLOT(smooth()));
    connect(gui_->comboBox_smooth_type, SIGNAL(currentIndexChanged(int)), this, SLOT(slotModeChanged(int)));
    connect(gui_->comboBox_integration_scheme, SIGNAL(currentIndexChanged(int)), this, SLOT(slotSchemeChanged(int)));
    connect(gui_->comboBox_visualize, SIGNAL(currentIndexChanged(int)), this, SLOT(slotVisualizeChanged(int)));

    QIcon* toolIcon = new QIcon(OpenFlipper::Options::iconDirStr()+OpenFlipper::Options::dirSeparator()+"smoother1.png");
    emit addToolbox( tr("AnisotropicMeanCurvature Smoother") , gui_, toolIcon );
}

void AnisotropicMeanCurvature::pluginsInitialized() {
    
    // Emit slot description
    emit setSlotDescription(tr("smooth(int)"),   tr("Smooth mesh using AnisotropicMeanCurvature."),
                            QStringList(tr("iterations")), QStringList(tr("Number of iterations")));
}


void AnisotropicMeanCurvature::
slotModeChanged(int _idx)
{
    int type = gui_->comboBox_smooth_type->currentIndex();

    if (type == 0)
        smooth_type = PrescribedMeanCurvature::ANISO_MEAN_CURVATURE;
    else if (type == 1)
        smooth_type = PrescribedMeanCurvature::PRESCRIBED_MEAN_CURVATURE;

}


void AnisotropicMeanCurvature::
slotSchemeChanged(int _idx)
{
    int type = gui_->comboBox_integration_scheme->currentIndex();

    if (type == 0)
        scheme = PrescribedMeanCurvature::EXPLICIT;
    else if (type == 1)
        scheme = PrescribedMeanCurvature::IMPLICIT;

}


void AnisotropicMeanCurvature::
slotVisualizeChanged(int _idx)
{
    int type = gui_->comboBox_visualize->currentIndex();

    if (type == 0)
        visualize = PrescribedMeanCurvature::NONE;
    else if (type == 1)
        visualize = PrescribedMeanCurvature::UPDATE_VECTOR;
    else if (type == 2)
        visualize = PrescribedMeanCurvature::COLOR_CODING;


    for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS) ;
          o_it != PluginFunctions::objectsEnd(); ++o_it)
    {

        if ( o_it->dataType( DATA_TRIANGLE_MESH ) )
        {
            TriMeshObject * meshObject = PluginFunctions::triMeshObject(o_it);
            if (visualize == PrescribedMeanCurvature::UPDATE_VECTOR)
            {
                pmc.showLineNode(meshObject);
            }else
            {
                pmc.clearLineNode(meshObject);
                if (visualize == PrescribedMeanCurvature::COLOR_CODING)
                {
                    TriMesh * mesh = meshObject->mesh();
                    ColorCoder coder(0, 10, false);
                    mesh->request_vertex_colors();

                    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
                    {
                        TriMesh::Point current = mesh->point(v_it);
                        TriMesh::Point source = mesh->property(source_points, v_it);
                        ACG::Vec3f color = coder.color_float(source|current)*255;
                        mesh->set_color(v_it, TriMesh::Color(color[0], color[1], color[2], 1));
                        printf("color %f %f %f \n", color[0], color[1], color[2]);
                    }
                    emit updatedObject( o_it->id(), UPDATE_COLOR );
                }

            }
        }
    }


}


/** \brief
 *
 *
 *
 */
void AnisotropicMeanCurvature::smooth() {

    int iterations = 1;
    
    if(!OpenFlipper::Options::nogui()) {
        iterations = gui_->spinBox_iterations->value();
    }

    pmc.set_lambda(gui_->doubleSpinBox->value());

    prescribedMeanCurvature(iterations);
}

/** \brief explicit anisotropic mean curvature
 *
 *
 *
 *
 * @param _iterations Number of iterations
 */
void AnisotropicMeanCurvature::smooth(int _iterations) {
    
    for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS) ;
          o_it != PluginFunctions::objectsEnd(); ++o_it)
    {

        bool selectionExists = false;

        if ( o_it->dataType( DATA_TRIANGLE_MESH ) )
        {

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

            double threshold = pmc.get_feature_threshold(mesh);

            for ( int i = 0 ; i < _iterations ; ++i )
            {

                for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
                {
                    mesh->property(smoothVector,v_it).vectorize(0.0f);
                    mesh->property(areaStar,v_it) = 0;
                    selectionExists |= mesh->status(v_it).selected();
                }


                for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
                {

                    for (TriMesh::VertexEdgeIter ve_it=mesh->ve_iter(v_it.handle()); ve_it; ++ve_it)
                    {

                        TriMesh::Normal edgeNormal;
                        TriMesh::Scalar meanCurvature = pmc.edge_mean_curvature_He_Ne(mesh, ve_it.handle()
                                                                                      , edgeNormal);
                        double weight = pmc.anisotropic_weight(meanCurvature, threshold
                                                               , PrescribedMeanCurvature::R);

                        mesh->property(smoothVector, v_it) += 0.5*meanCurvature*weight*edgeNormal;

                    }

                    for (TriMesh::VertexFaceIter vf_it=mesh->vf_iter(v_it.handle()); vf_it; ++vf_it)
                    {

                        mesh->property(areaStar, v_it) += pmc.face_area(mesh, vf_it.handle());
                    }

                }

                for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
                {
                    if(selectionExists && mesh->status(v_it).selected() == false) {
                        continue;
                    }

                    TriMesh::Scalar area = mesh->property(areaStar, v_it);
                    TriMesh::Normal updateVector = mesh->property(smoothVector, v_it);

                    mesh->set_point(v_it, mesh->point(v_it) - (3*PrescribedMeanCurvature::TIME_STEP/area)*updateVector);
                }

                mesh->update_normals();


            }// Iterations end


            mesh->update_normals();



            pmc.updateLineNode(meshObject, smoothVector, areaStar);

            // Remove the property
            mesh->remove_property( smoothVector );
            mesh->remove_property( areaStar );


            emit updatedObject( o_it->id(), UPDATE_ALL );

            // Create backup
            emit createBackup(o_it->id(), "ExplicitAnisotropicMeanCurvature Smoothing", UPDATE_ALL );


        } else {

            emit log(LOGERR, "DataType not supported.");
        }
    }

    emit updateView();
}



void AnisotropicMeanCurvature::prescribedMeanCurvature(int _iterations)
{

    for ( PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS) ;
          o_it != PluginFunctions::objectsEnd(); ++o_it)
    {

        if ( o_it->dataType( DATA_TRIANGLE_MESH ) )
        {


            TriMeshObject * meshObject = PluginFunctions::triMeshObject(o_it);
            TriMesh* mesh = meshObject->mesh();

            for ( PluginFunctions::ObjectIterator osrc_it(PluginFunctions::SOURCE_OBJECTS) ;
                  osrc_it != PluginFunctions::objectsEnd(); ++osrc_it)
            {
                TriMeshObject * sourceObj = PluginFunctions::triMeshObject(osrc_it);
                TriMesh * sourceMesh = sourceObj->mesh();
                mesh->add_property( source_points, "source_points" );
                for (TriMesh::VertexIter v_it=mesh->vertices_begin(), srcv_it=sourceMesh->vertices_begin();
                     v_it!=mesh->vertices_end(); ++v_it, ++srcv_it)
                {
                    mesh->property(source_points,v_it) = sourceMesh->point(srcv_it);
                }
            }

            if (scheme == PrescribedMeanCurvature::EXPLICIT
                    && smooth_type == PrescribedMeanCurvature::PRESCRIBED_MEAN_CURVATURE)
            {
                pmc.smooth_explicit_pmc(_iterations, meshObject, visualize);
            }else if (scheme == PrescribedMeanCurvature::IMPLICIT
                      && smooth_type == PrescribedMeanCurvature::PRESCRIBED_MEAN_CURVATURE)
            {
                return;
            }else
            {
                pmc.smooth_aniso(_iterations, meshObject, smooth_type, scheme, visualize);
            }

            emit updatedObject( o_it->id(), UPDATE_ALL );

            // Create backup
            emit createBackup(o_it->id(), "ExplicitAnisotropicMeanCurvature Smoothing", UPDATE_ALL );

        }
        else
        {

            emit log(LOGERR, "DataType not supported.");
        }

    }

    emit updateView();
}






Q_EXPORT_PLUGIN2( explicitAnisotropicMeanCurvature , AnisotropicMeanCurvature );

