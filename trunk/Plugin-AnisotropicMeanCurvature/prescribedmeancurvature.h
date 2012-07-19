#ifndef PRESCRIBEDMEANCURVATURE_H
#define PRESCRIBEDMEANCURVATURE_H

#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ACG/Scenegraph/LineNode.hh>
#include <ACG/Math/VectorT.hh>
#include "Mat3x3.hh"



/**
 * @brief   This class implements the smoothing method from the paper:
 *   "Anisotropic Filtering of Non-Linear Surface Features".
 *
 *  for semi implicit method we need to express
 * Ha = Ka^j.P^j+1, where j is time step
 * Ka can be computed using the information from step j
 * but it is applied into the vertices from step j+1 and never into step j
 * Consequently, Ha is semi-implicit.
 * We have to determine which part of Ha is from j (Ka component) and which is from j+1 (P component).
 * and they must not be inter-dependent on each other to make the method linear.
 * Otherwise, we have to use iterative solver.
 * Obviously, we want as many components being j+1 as possible.

 * To build the matrix Ka, we need neighborhood map from a point to its neighbors.
 * Ka has the size of 3n*3n where 3 consecutive rows correspond to 1 vertex.
 * In these rows, we use neighborhood map to refer to its neighbors.
 * Thus column indices refer to its neighbors.
 *
 */
class PrescribedMeanCurvature
{
public:

    PrescribedMeanCurvature();

    //! radius of the smoothness of the weight
    static const double R = 10;

    enum SmoothingMode { ANISO_MEAN_CURVATURE, PRESCRIBED_MEAN_CURVATURE };
    enum IntegrationScheme { INDIVIDUAL_UPDATE, BATCH_UPDATE };
    enum VisualizeMode {NONE, UPDATE_VECTOR, COLOR_CODING};

    /**
     * @brief               put less smoothing weight on feature edge based on edge mean curvature.
     *
     * @param curvature     edge mean curvature He
     * @param lambda        threshold to classify feature edges
     * @param r             radius of the smoothness of the weight
     * @return double       anisotropic weight of the edge mean curvature He
     */
    double anisotropic_weight(double curvature, double lambda, double r);

    void smooth_explicit_pmc(int _iterations, TriMeshObject * meshObject
                             , VisualizeMode visualize, double time_step);

    double get_lambda()
    {
        return m_lambda;
    }

    void set_lambda(double lambda)
    {
        m_lambda = lambda;
    }

    double get_feature_threshold(TriMesh * mesh);

    void updateLineNode(TriMeshObject * _meshObject
                        , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex);

    void clearLineNode(TriMeshObject * _meshObject);

    void showLineNode(TriMeshObject * _meshObject);

    /**
     * @brief                       calculate face area.
     *
     * @param _mesh                 the mesh
     * @param fh                    face handle
     * @return TriMesh::Scalar      face area
     */
    TriMesh::Scalar face_area(TriMesh *_mesh, TriMesh::FaceHandle fh);

    /**
     * @brief           compute edge mean curvature.
     *
     * @param _mesh     the mesh
     * @param _eh       edge handle to extract edge information from the mesh
     * @param normal    edge normal
     * @return double   edge mean curvature value He
     */
    double edge_mean_curvature_He_Ne(TriMesh * _mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal);

protected:

    typedef ACG::Vec3uc Color;
    typedef ACG::Vec3d  Vec3d;

    double m_lambda;



    /**
     * @brief   calculate vertex angle at p.
     *
     * @param p         center vertex
     * @param q         neighbor to p
     * @param r         neighbor to p
     * @return double   vertex angle at p.
     */
    double cal_angle(const TriMesh::Point & p, const TriMesh::Point & q, const TriMesh::Point & r);

    /**
     * @brief   replaces volume gradient of a vertex by its anisotropic one if the vertex is a feature vertex
     *
     * @param _mesh             The mesh.
     * @param amc_Ha            aniso. mean curvature Ha.
     * @param smoothed_amc_Ha   smoothed amc_Ha for temporary use
     * @param volume_gradient_Va    returned volume gradient
     * @param is_feature    define if Ha is different from the mean curvature H computed using cotangent formula
     */
    void smooth_amc(TriMesh *_mesh
                    , const OpenMesh::VPropHandleT< TriMesh::Normal > & amc_Ha
                    , const OpenMesh::VPropHandleT< TriMesh::Normal > & smoothed_amc_Ha
                    , const OpenMesh::VPropHandleT< TriMesh::Normal > & volume_gradient_Va
                    , const OpenMesh::VPropHandleT< bool > & is_feature);

    /**
     * @brief   compute prescribed mean curvature function value f at a vertex p.
     *          f = M^-1*Ha but is approximated by Ha/|Va|
     *
     * @param _mesh     the mesh
     * @param amc_Ha    anisotropic mean curvature Ha
     * @param volume_gradient_Va   volume gradient Va
     * @param is_feature    is a feature vertex
     * @param apmc_function_f   aniso. prescribed mean curvature function f
     * @param smoothed_apmc_function_f  the smoothed version of f
     */
    void compute_apmc_function_f(TriMesh *_mesh
                            , const OpenMesh::VPropHandleT< TriMesh::Normal > & amc_Ha
                            , const OpenMesh::VPropHandleT< TriMesh::Normal > & volume_gradient_Va
                            , const OpenMesh::VPropHandleT< bool > & is_feature
                            , const OpenMesh::VPropHandleT< double > & apmc_function_f
                            , const OpenMesh::VPropHandleT< double > & smoothed_apmc_function_f);

    double area_star_edge(TriMesh *_mesh, TriMesh::EdgeHandle _eh);

    double calculate_cross_matrix_Ax_qjpi(TriMesh *_mesh, const TriMesh::EdgeHandle & _eh
                                          , TriMesh::Scalar & normal_length
                                          , const TriMesh::VertexHandle & _vh, Mat3x3 & cross);

    double cal_edge_norm_deriv_qjpi(TriMesh *_mesh, const TriMesh::HalfedgeHandle & hh1
                                    , TriMesh::Scalar & normal_length
                                    , std::vector<Mat3x3> & derivatives);

    double get_edge_length(TriMesh * mesh);

    /**
     * @brief   circulate around the oriented face using FaceHalfEdgeIter to extract p, q, r in CCW, take (qxr)/6.
     *
     * @param _mesh     the mesh.
     * @param fh        face handle.
     * @param vh        vertex handle for the point p.
     * @param gradient  returned volume gradient.
     */
    void volume_gradient(TriMesh *_mesh, TriMesh::FaceHandle fh
                         , TriMesh::VertexHandle vh, TriMesh::Normal & gradient);

    ACG::SceneGraph::LineNode * getLineNode(TriMeshObject * _meshObject);
    void addLine( ACG::SceneGraph::LineNode * _line_node, Vec3d _p0, Vec3d _p1, Color _color );

};

#endif // PRESCRIBEDMEANCURVATURE_H
