#ifndef PRESCRIBEDMEANCURVATURE_H
#define PRESCRIBEDMEANCURVATURE_H

#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ACG/Scenegraph/LineNode.hh>
#include <ACG/Math/VectorT.hh>

class PrescribedMeanCurvature
{
public:
    PrescribedMeanCurvature();


    static const double R = 10;
    //stable time step for explicit methods
    //for implicit method multiply by a factor of 10 or 100
    //for fine scale noise multiply by a factor of 0.1
    static const double TIME_STEP = 0.00001;

    TriMesh::Scalar faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh);

    void volumeGradient(TriMesh *_mesh, TriMesh::FaceHandle fh, TriMesh::VertexHandle vh, TriMesh::Normal & gradient);

    double edgeMeanCurvature(TriMesh * _mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal);
    double edgeMeanCurvature(PolyMesh * _mesh, PolyMesh::EdgeHandle _eh, TriMesh::Normal & normal);

    double anisotropicWeight(double curvature, double lambda, double r);

    void smooth(int _iterations, TriMeshObject * meshObject);

    void updateLineNode(TriMeshObject * _meshObject, OpenMesh::VPropHandleT< TriMesh::Normal > & anisoMeanCurvature, OpenMesh::VPropHandleT< double >& areaStar);

private:

    typedef ACG::Vec3uc Color;
    typedef ACG::Vec3d  Vec3d;

    double calAngle(TriMesh::Point p, TriMesh::Point q, TriMesh::Point r);
    void smoothAnisotropicMeanCurvature(TriMesh *_mesh, OpenMesh::VPropHandleT< TriMesh::Normal > & anisoMeanCurvature, OpenMesh::VPropHandleT< TriMesh::Normal > & smoothedAMC);

    ACG::SceneGraph::LineNode * getLineNode(TriMeshObject * _meshObject);
    void addLine( ACG::SceneGraph::LineNode * _line_node, Vec3d _p0, Vec3d _p1, Color _color );

};

#endif // PRESCRIBEDMEANCURVATURE_H
