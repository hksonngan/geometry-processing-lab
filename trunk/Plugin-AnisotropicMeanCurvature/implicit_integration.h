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
#include "Mat3x3.hh"

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


//    static const double time_step = 0.0001;
//    static const double time_step = 0.001;


    /**
     * @brief establish the mass matrix and the matrix for calculating aniso. mean curvature Ka.
     *  solve the system of linear equations.
     *
     * @param mesh          the mesh
     * @param mesh_size     the number of vertices
     * @param area_star     vertex property containing area star
     * @param vertex_id     vertex property containing vertex id
     */
    void compute_semi_implicit_integration(TriMesh * mesh
                                           , unsigned int mesh_size
                                           , const OpenMesh::VPropHandleT< double > & area_star
                                           , const OpenMesh::VPropHandleT< int > & vertex_id
                                           , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex
                                           , double time_step
                                           , bool is_lumped_mass = false
                                           );







    void compute_explicit_integration_with_mass(TriMesh *mesh
                                                , unsigned int mesh_size
                                                , const OpenMesh::VPropHandleT< double > & area_star
                                                , const OpenMesh::VPropHandleT< int > & vertex_id
                                                , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex
                                                , double time_step
                                                , bool is_lumped_mass = false
                                                );


    void compute_taylor_semi_implicit(TriMesh *mesh
                                      , unsigned int mesh_size
                                      , const OpenMesh::VPropHandleT< double > & area_star
                                      , const OpenMesh::VPropHandleT< int > & vertex_id
                                      , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex
                                      , double time_step
                                      , bool is_lumped_mass = false
                                      );



private:

    /**
     * @brief           populate a dense vector of vertices.
     *
     * @param mesh
     * @param vertices  this has the length of 3*n where n is the number of vertices
     */
    void init_vertex_vector(TriMesh * mesh, Eigen::VectorXd & vertices
                            , const OpenMesh::VPropHandleT< int > & vertex_id);

    void init_mass_matrix(TriMesh * mesh, PrescribedMeanCurvature * pmc
                          , const OpenMesh::VPropHandleT< double > & area_star
                          , const OpenMesh::VPropHandleT< int > & vertex_id
                          , Eigen::SparseMatrix<double> & mass);

    void init_amc_matrix(TriMesh * mesh, PrescribedMeanCurvature * pmc
                         , const OpenMesh::VPropHandleT< int > & vertex_id
                         , Eigen::SparseMatrix<double> & amc_matrix);


    void init_lumped_mass_matrix(TriMesh *mesh , PrescribedMeanCurvature * pmc
                     , const OpenMesh::VPropHandleT< double > & area_star
                     , const OpenMesh::VPropHandleT< int > & vertex_id
                     , Eigen::SparseMatrix<double> & mass);


    void init_Jacobian(TriMesh *mesh, PrescribedMeanCurvature *pmc
                    , const OpenMesh::VPropHandleT< int > &vertex_id
                    , Eigen::SparseMatrix<double> &jacobian);


    void init_lumped_mass_matrix_inverted(TriMesh *mesh , PrescribedMeanCurvature * pmc
                     , const OpenMesh::VPropHandleT< double > & area_star
                     , const OpenMesh::VPropHandleT< int > & vertex_id
                     , Eigen::SparseMatrix<double> & mass);


};

#endif // IMPLICIT_INTEGRATION_H
