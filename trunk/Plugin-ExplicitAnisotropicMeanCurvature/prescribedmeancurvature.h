#ifndef PRESCRIBEDMEANCURVATURE_H
#define PRESCRIBEDMEANCURVATURE_H

#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

class PrescribedMeanCurvature
{
public:
    PrescribedMeanCurvature();



    double edgeMeanCurvature(TriMesh * _mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal);


    double anisotropicWeight(double curvature, double lambda, double r);





    TriMesh::Scalar faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh);

};

#endif // PRESCRIBEDMEANCURVATURE_H
