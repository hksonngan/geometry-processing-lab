#ifndef PRESCRIBEDMEANCURVATURE_H
#define PRESCRIBEDMEANCURVATURE_H

#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
//#include "ExplicitAnisotropicMeanCurvature.hh"

class PrescribedMeanCurvature
{
public:
    PrescribedMeanCurvature();


    static const double R = 10;

    TriMesh::Scalar faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh);

    double edgeMeanCurvature(TriMesh * _mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal);
    double edgeMeanCurvature(PolyMesh * _mesh, PolyMesh::EdgeHandle _eh, TriMesh::Normal & normal);

    double anisotropicWeight(double curvature, double lambda, double r);

    void smooth(int _iterations);


};

#endif // PRESCRIBEDMEANCURVATURE_H
