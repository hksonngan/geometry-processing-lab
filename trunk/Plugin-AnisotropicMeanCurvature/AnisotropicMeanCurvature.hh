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
*   $Revision: 13361 $                                                       *
*   $LastChangedBy: phan-anh $                                                *
*   $Date: 2012-01-12 16:33:16 +0100 (Thu, 12 Jan 2012) $                     *
*                                                                            *
\*===========================================================================*/




#ifndef ANISOTROPICMEANCURVATURE_HH
#define ANISOTROPICMEANCURVATURE_HH

#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/ScriptInterface.hh>
#include <OpenFlipper/BasePlugin/BackupInterface.hh>
#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "prescribedmeancurvature.h"
#include "SmootherToolboxWidget.hh"
#include "implicit_integration.h"
#include <limits>

using namespace std;

class AnisotropicMeanCurvature : public QObject, BaseInterface, ToolboxInterface
        , LoggingInterface, ScriptInterface, BackupInterface
{
    Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(ScriptInterface)
    Q_INTERFACES(BackupInterface)

signals:
    //BaseInterface
    void updateView();
    void updatedObject(int _id, const UpdateType& _type);
    void setSlotDescription(QString     _slotName,   QString     _slotDescription,
                            QStringList _parameters, QStringList _descriptions);

    //LoggingInterface
    void log(Logtype _type, QString _message);
    void log(QString _message);
    
    // ToolboxInterface
    void addToolbox( QString _name  , QWidget* _widget, QIcon* _icon);
    
    // ScriptInterface
    void scriptInfo(QString _functionName);
    
    // BackupInterface
    void createBackup( int _id , QString _name, UpdateType _type = UPDATE_ALL );

public:

    AnisotropicMeanCurvature();
    ~AnisotropicMeanCurvature();

    // BaseInterface
    QString name() { return (QString("Explicit Anisotropic Mean Curvature Smoother")); };
    QString description( ) { return (QString("Smooths the active Mesh")); };






private:

    OpenMesh::VPropHandleT< TriMesh::Point > source_points;

    /// Widget for Toolbox
    SmootherToolboxWidget * gui_;

    double color_range;


    PrescribedMeanCurvature::SmoothingMode smooth_type;
    PrescribedMeanCurvature::IntegrationScheme scheme;
    PrescribedMeanCurvature::VisualizeMode visualize;

    PrescribedMeanCurvature pmc;

    bool attach_source(TriMeshObject * meshObject);
    void recompute_color(TriMeshObject * meshObject, int object_id);
    void add_noise(TriMesh::Point & point, const TriMesh::Normal & normal, double range);
    bool add_noise(TriMesh* mesh);

private slots:
    void smooth();

    void slotModeChanged(int);

    void slotSchemeChanged(int _idx);

    void slot_add_noise();

    void slotVisualizeChanged(int);
    
    void initializePlugin(); // BaseInterface
    
    void pluginsInitialized(); // BaseInterface
    
    // Scriptable functions
public slots:

    void smooth(int _iterations, double time_step);

    void prescribedMeanCurvature(int _iterations, double time_step);

    QString version() { return QString("1.0"); };
};

#endif //ANISOTROPICMEANCURVATURE_HH
