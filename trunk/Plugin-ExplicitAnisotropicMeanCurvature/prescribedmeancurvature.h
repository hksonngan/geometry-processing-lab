#ifndef PRESCRIBEDMEANCURVATURE_H
#define PRESCRIBEDMEANCURVATURE_H

#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "ExplicitAnisotropicMeanCurvature.hh"

class PrescribedMeanCurvature: public ExplicitAnisotropicMeanCurvature
{
public:
    PrescribedMeanCurvature();

protected:


    void smooth(int _iterations);


};

#endif // PRESCRIBEDMEANCURVATURE_H
