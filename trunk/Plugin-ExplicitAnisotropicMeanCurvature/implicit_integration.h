#ifndef IMPLICIT_INTEGRATION_H
#define IMPLICIT_INTEGRATION_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
//install cholmod first
#include <unsupported/Eigen/CholmodSupport>
#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "prescribedmeancurvature.h"

/**
 * @brief perform semi-implicit integration scheme for prescribed mean curvature mesh smoothing.
 *
 * vertices of the mesh are aggregated into a global vector as follows:
 * (p11, p12, p13, ... pn1, pn2, pn3) where pi1, pi2, pi3 are the 3 coordinates of the vertex pi.
 * Consequently, the mass matrix will have the size of 3n*3n
 * where each entry is a 3*3 diagonal matrix storing the same value Mij.
 */
class Implicit_Integration
{
public:
    Implicit_Integration();



private:

    void init_vertex_vector(TriMesh * mesh, Eigen::VectorXd & vertices);

    void init_mass_matrix(TriMesh * mesh, PrescribedMeanCurvature * pmc
                          , Eigen::SparseMatrix<double> & mass);

};

#endif // IMPLICIT_INTEGRATION_H
