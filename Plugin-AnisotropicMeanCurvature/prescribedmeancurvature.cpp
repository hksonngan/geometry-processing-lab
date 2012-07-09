#include "prescribedmeancurvature.h"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
//install cholmod first
#include <unsupported/Eigen/CholmodSupport>
#include "implicit_integration.h"

PrescribedMeanCurvature::PrescribedMeanCurvature()
{
    m_lambda = 0.8;
}



void PrescribedMeanCurvature::
smooth_explicit_pmc(int _iterations, TriMeshObject * meshObject
                    , VisualizeMode visualize)
{


    bool selectionExists = false;


    TriMesh* mesh = meshObject->mesh();

    // Property for the active mesh to store original point positions
    OpenMesh::VPropHandleT< TriMesh::Normal > amc_Ha;
    //the smoothed_amc_Ha is for temporary use
    OpenMesh::VPropHandleT< TriMesh::Normal > smoothed_amc_Ha;
    OpenMesh::VPropHandleT< double > area_star;
    OpenMesh::VPropHandleT< double > apmc_function_f;
    OpenMesh::VPropHandleT< double > smoothed_apmc_function_f;
    OpenMesh::VPropHandleT< bool > is_feature;
    OpenMesh::VPropHandleT< TriMesh::Normal > volume_gradient_Va;

    // Add a property to the mesh to store mean curvature and area
    mesh->add_property( amc_Ha, "anisotropic_mean_curvature_Ha" );
    mesh->add_property( smoothed_amc_Ha, "smoothed_mnisotropic_mean_curvature_Ha" );
    mesh->add_property( area_star, "area_star" );
    mesh->add_property( apmc_function_f, "aniso_prescribed_mean_curvature_function_f" );
    mesh->add_property( smoothed_apmc_function_f, "smoothed_apmc_function_f" );
    mesh->add_property( is_feature, "is_feature" );
    mesh->add_property( volume_gradient_Va, "volume_gradient_Va" );

    mesh->request_vertex_normals();
    mesh->request_vertex_colors();
    mesh->request_face_normals();

    unsigned int count(meshObject->mesh()->n_vertices());

    mesh->update_normals();

    double threshold = get_feature_threshold(mesh);


    for ( int i = 0 ; i < _iterations ; ++i )
    {

        unsigned int noFeatureVertices = 0;

        for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
        {
            mesh->property(amc_Ha,v_it).vectorize(0.0f);
            mesh->property(smoothed_amc_Ha,v_it).vectorize(0.0f);
            mesh->property(volume_gradient_Va,v_it).vectorize(0.0f);
            mesh->property(area_star,v_it) = 0;
            mesh->property(apmc_function_f,v_it) = 0;
            mesh->property(smoothed_apmc_function_f,v_it) = 0;
            mesh->property(is_feature,v_it) = false;
            selectionExists |= mesh->status(v_it).selected();
        }


        for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
        {
            TriMesh::Normal isotropic;
            isotropic.vectorize(0);

            for (TriMesh::VertexEdgeIter ve_it=mesh->ve_iter(v_it.handle()); ve_it; ++ve_it)
            {

                TriMesh::Normal edgeNormal;
                TriMesh::Scalar meanCurvature = edge_mean_curvature_He_Ne(mesh, ve_it.handle(), edgeNormal);
                double weight = anisotropic_weight(meanCurvature, threshold, R);

                mesh->property(amc_Ha, v_it) += 0.5*meanCurvature*weight*edgeNormal;

                isotropic += 0.5*meanCurvature*edgeNormal;

            }

            if (mesh->property(amc_Ha, v_it) != isotropic)
            {
                mesh->property(is_feature,v_it) = true;
                noFeatureVertices++;
            }

            for (TriMesh::VertexFaceIter vf_it=mesh->vf_iter(v_it.handle()); vf_it; ++vf_it)
            {

                mesh->property(area_star, v_it) += face_area(mesh, vf_it.handle());
                //volume gradient
                TriMesh::Normal volGrad;
                volume_gradient(mesh, vf_it.handle(), v_it, volGrad);
                mesh->property(volume_gradient_Va, v_it) += volGrad;
            }

        }

        smooth_amc(mesh, amc_Ha, smoothed_amc_Ha, volume_gradient_Va, is_feature);

        compute_apmc_function_f(mesh, amc_Ha, volume_gradient_Va
                                , is_feature, apmc_function_f, smoothed_apmc_function_f);

        printf("number of feature vertices: %d in total %d \n", noFeatureVertices, count);

        for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
        {
            if(mesh->is_boundary(v_it.handle())) {
                continue;
            }

            TriMesh::Scalar area = mesh->property(area_star, v_it);
            TriMesh::Normal updateVector = mesh->property(amc_Ha, v_it);

            TriMesh::Normal result = updateVector -
                    mesh->property(smoothed_apmc_function_f, v_it)*mesh->property(volume_gradient_Va, v_it);

            mesh->set_point(v_it, mesh->point(v_it) - (3*TIME_STEP/area)*result);
        }

        mesh->update_normals();


    }// Iterations end


    mesh->update_normals();

    //if (visualize == UPDATE_VECTOR) updateLineNode(meshObject, amc_Ha, area_star);
    updateLineNode(meshObject, amc_Ha, area_star);
    if (visualize != UPDATE_VECTOR) this->clearLineNode(meshObject);


    // Remove the property
    mesh->remove_property( amc_Ha );
    mesh->remove_property( area_star );
    mesh->remove_property( is_feature );
    mesh->remove_property( volume_gradient_Va );
    mesh->remove_property( smoothed_amc_Ha );
    mesh->remove_property( smoothed_apmc_function_f );
    mesh->remove_property( apmc_function_f );


}




void PrescribedMeanCurvature::
smooth_aniso(int _iterations, TriMeshObject * meshObject
            , SmoothingMode smooth_type
            , IntegrationScheme scheme
            , VisualizeMode visualize)
{

    TriMesh* mesh = meshObject->mesh();

    // Property for the active mesh to store original point positions
    OpenMesh::VPropHandleT< TriMesh::Point > old_vertex;
    OpenMesh::VPropHandleT< double > area_star;
    OpenMesh::VPropHandleT< int > vertex_id;

    // Add a property to the mesh to store mean curvature and area
    mesh->add_property( old_vertex, "new_amc_vertex" );
    mesh->add_property( area_star, "area_star" );
    mesh->add_property( vertex_id, "vertex_id" );


    mesh->request_vertex_normals();
    mesh->request_vertex_colors();
    mesh->request_face_normals();

    unsigned int count(meshObject->mesh()->n_vertices());


    mesh->update_normals();




    //initialize vertex id and its neighbors
    int idx = 0;
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        mesh->property(vertex_id, v_it) = idx;
        idx++;
    }


    for ( int i = 0 ; i < _iterations ; ++i )
    {

        for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
        {
            mesh->property(old_vertex,v_it) = TriMesh::Point(0, 0, 0);

            mesh->property(area_star,v_it) = 0;

        }

        for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
        {


            for (TriMesh::VertexFaceIter vf_it=mesh->vf_iter(v_it.handle()); vf_it; ++vf_it)
            {

                mesh->property(area_star, v_it) += face_area(mesh, vf_it.handle());

            }

        }

        Implicit_Integration implicit_solver;
        //implicit_solver.compute_semi_implicit_integration(mesh, count, area_star, vertex_id, old_vertex);

        if (scheme == EXPLICIT
                && smooth_type == ANISO_MEAN_CURVATURE)
        {
            implicit_solver.compute_explicit_integration_with_mass(mesh, count, area_star, vertex_id, old_vertex);
        }

        if (scheme == IMPLICIT
                && smooth_type == ANISO_MEAN_CURVATURE)
        {
            implicit_solver.compute_taylor_semi_implicit(mesh, count, area_star, vertex_id, old_vertex);
        }

        mesh->update_normals();


    }// Iterations end


    mesh->update_normals();
    //if (visualize == UPDATE_VECTOR) updateLineNode(meshObject, old_vertex);
    updateLineNode(meshObject, old_vertex);
    if (visualize != UPDATE_VECTOR) this->clearLineNode(meshObject);


    // Remove the property
    mesh->remove_property( area_star );
    mesh->remove_property( vertex_id );
    mesh->remove_property( old_vertex );

}




double PrescribedMeanCurvature::get_feature_threshold(TriMesh *mesh)
{
    double count = 0;
    double threshold = 0;
    for (TriMesh::EdgeIter e_it=mesh->edges_begin(); e_it!=mesh->edges_end(); ++e_it)
    {
        threshold += mesh->calc_edge_length(e_it);
        count++;
    }
    threshold = threshold*m_lambda*2/count;
    return threshold;
}



double PrescribedMeanCurvature::
cal_edge_norm_deriv_qjpi(TriMesh *_mesh, const TriMesh::HalfedgeHandle & hh1
                         , TriMesh::Scalar & normal_length
                         , std::vector<Mat3x3> & derivatives)
{

    double dihedral = _mesh->calc_dihedral_angle(hh1);

    TriMesh::HalfedgeHandle  hh2;
    TriMesh::FaceHandle    fh1, fh2;
    TriMesh::Normal normal;
    TriMesh::Point qj, pi, q_plus, q_minus;

    hh2 = _mesh->opposite_halfedge_handle(hh1);
    fh1 = _mesh->face_handle(hh1);
    fh2 = _mesh->face_handle(hh2);

    TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
    TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);
    normal = (n1+n2);
    normal_length = normal.norm();

    qj = _mesh->point(_mesh->to_vertex_handle(hh1));
    pi = _mesh->point(_mesh->from_vertex_handle(hh1));
    q_plus = _mesh->point(_mesh->to_vertex_handle(_mesh->next_halfedge_handle(hh1)));
    q_minus = _mesh->point(_mesh->from_vertex_handle(_mesh->prev_halfedge_handle( hh2 ) ) );

    TriMesh::Point cross_A = q_plus - q_minus;
    TriMesh::Point cross_B = qj - pi;

    if (dihedral < 0)
    {
        cross_A *= -1;
        cross_B *= -1;
    }

    TriMesh::Point edge_dev_pi = cross_A;
    TriMesh::Point edge_dev_qj = -cross_A;
    TriMesh::Point edge_dev_q_plus = cross_B;
    TriMesh::Point edge_dev_q_minus = -cross_B;


    derivatives.push_back(Mat3x3(0, -edge_dev_pi[2], edge_dev_pi[1],
                                 edge_dev_pi[2], 0, -edge_dev_pi[0],
                                 -edge_dev_pi[1], edge_dev_pi[0], 0
                                 ));

    derivatives.push_back(Mat3x3(0, -edge_dev_q_minus[2], edge_dev_q_minus[1],
                                 edge_dev_q_minus[2], 0, -edge_dev_q_minus[0],
                                 -edge_dev_q_minus[1], edge_dev_q_minus[0], 0
                                 ));

    derivatives.push_back(Mat3x3(0, -edge_dev_qj[2], edge_dev_qj[1],
                                 edge_dev_qj[2], 0, -edge_dev_qj[0],
                                 -edge_dev_qj[1], edge_dev_qj[0], 0
                                 ));

    derivatives.push_back(Mat3x3(0, -edge_dev_q_plus[2], edge_dev_q_plus[1],
                                 edge_dev_q_plus[2], 0, -edge_dev_q_plus[0],
                                 -edge_dev_q_plus[1], edge_dev_q_plus[0], 0
                                 ));


    double edgeLength = _mesh->calc_edge_length(hh1);
    return 2.0*edgeLength*cos(dihedral/2.0);

}










double PrescribedMeanCurvature::
calculate_cross_matrix_Ax_qjpi(TriMesh *_mesh, const TriMesh::EdgeHandle & _eh, TriMesh::Scalar & normal_length
                               , const TriMesh::VertexHandle & _vh, Mat3x3 & cross)
{


    //dihedral = 0 on the boundary
    double dihedral = 0;

    TriMesh::HalfedgeHandle  hh1, hh2;
    TriMesh::FaceHandle    fh1, fh2;
    TriMesh::Normal normal, the_edge;
    TriMesh::Point qj, pi;

    hh1 = _mesh->halfedge_handle(_eh, 0);
    hh2 = _mesh->halfedge_handle(_eh, 1);

    //normal is the edge vector rotated by 90 degree
    if ( _mesh->is_boundary( hh1 ) )
    {

        fh2 = _mesh->face_handle(hh2);
        TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);
        TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh2)) - _mesh->point(_mesh->from_vertex_handle(hh2));
        //rotation by -90, normal of 2d triangle edge pointing outward
        normal = (edgeVector%n2);
        normal_length = normal.norm();

    }else if ( _mesh->is_boundary( hh2 ) )
    {

        fh1 = _mesh->face_handle(hh1);
        TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
        TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh1)) - _mesh->point(_mesh->from_vertex_handle(hh1));
        normal = (edgeVector%n1);
        normal_length = normal.norm();

    }else
    {

        fh1 = _mesh->face_handle(hh1);
        fh2 = _mesh->face_handle(hh2);

        TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
        TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);

        normal = (n1+n2);
        normal_length = normal.norm();

        dihedral = _mesh->calc_dihedral_angle(_eh);

        if (_vh == _mesh->to_vertex_handle(hh1))
        {
            pi = _mesh->point(_mesh->to_vertex_handle(hh1));
            qj = _mesh->point(_mesh->from_vertex_handle(hh1));
        }else
        {
            pi = _mesh->point(_mesh->from_vertex_handle(hh1));
            qj = _mesh->point(_mesh->to_vertex_handle(hh1));
        }

        the_edge = qj - pi;

        //in both cases, the vertex is pushed in the tangental direction
        //without this flipping of normal, it will get worse
        if (dihedral < 0) the_edge *= -1;

        cross = Mat3x3(0, -the_edge[2], the_edge[1],
                       the_edge[2], 0, -the_edge[0],
                       -the_edge[1], the_edge[0], 0);

    }

    double edgeLength = _mesh->calc_edge_length(_eh);
    return 2.0*edgeLength*cos(dihedral/2.0);

}



double PrescribedMeanCurvature::area_star_edge(TriMesh *_mesh, TriMesh::EdgeHandle _eh)
{

    TriMesh::HalfedgeHandle  hh1, hh2;
    TriMesh::FaceHandle    fh1, fh2;
    double area_star = 0;

    hh1 = _mesh->halfedge_handle(_eh, 0);
    hh2 = _mesh->halfedge_handle(_eh, 1);

    if ( _mesh->is_boundary( hh1 ) )
    {
        fh2 = _mesh->face_handle(hh2);
        area_star = face_area(_mesh, fh2);

    }else if ( _mesh->is_boundary( hh2 ) )
    {
        fh1 = _mesh->face_handle(hh1);
        area_star = face_area(_mesh, fh1);
    }else
    {
        fh1 = _mesh->face_handle(hh1);
        fh2 = _mesh->face_handle(hh2);

        area_star += face_area(_mesh, fh2);
        area_star += face_area(_mesh, fh1);
    }

    return area_star;
}






















void PrescribedMeanCurvature::
compute_apmc_function_f(TriMesh *_mesh
                        , const OpenMesh::VPropHandleT< TriMesh::Normal > & amc_Ha
                        , const OpenMesh::VPropHandleT< TriMesh::Normal > & volume_gradient_Va
                        , const OpenMesh::VPropHandleT< bool > & is_feature
                        , const OpenMesh::VPropHandleT< double > & apmc_function_f
                        , const OpenMesh::VPropHandleT< double > & smoothed_apmc_function_f)
{
    for (TriMesh::VertexIter v_it=_mesh->vertices_begin(); v_it!=_mesh->vertices_end(); ++v_it)
    {
        TriMesh::Normal volGrad = _mesh->property(volume_gradient_Va, v_it);
        TriMesh::Normal anisoMCVec = _mesh->property(amc_Ha, v_it);
        TriMesh::Scalar pmc = anisoMCVec.norm()/volGrad.norm();

        _mesh->property(apmc_function_f, v_it) = pmc;

    }

    //smooth by averaging over the neighborhood
    for (TriMesh::VertexIter v_it=_mesh->vertices_begin(); v_it!=_mesh->vertices_end(); ++v_it)
    {

        bool feature = _mesh->property(is_feature, v_it);
        TriMesh::Scalar avg = 0;
        int i = 0;

        TriMesh::VertexVertexIter vv_iter(*_mesh,v_it);
        for ( ; vv_iter; ++vv_iter )
        {
            //if the center vertex is feature vertex we only average over feature neighbors
            if (feature && !_mesh->property(is_feature, vv_iter)) continue;

            avg += _mesh->property(apmc_function_f, vv_iter);
            i++;

        }

        avg /= i;
        _mesh->property(smoothed_apmc_function_f, v_it) = avg;

    }

}


void PrescribedMeanCurvature::
smooth_amc(TriMesh *_mesh
           , const OpenMesh::VPropHandleT< TriMesh::Normal > & amc_Ha
           , const OpenMesh::VPropHandleT< TriMesh::Normal > & smoothed_amc_Ha
           , const OpenMesh::VPropHandleT< TriMesh::Normal > & volume_grad
           , const OpenMesh::VPropHandleT< bool > & is_feature)
{
    //only replace volume gradient by its aniso counterpart at feature vertices
    for (TriMesh::VertexIter v_it=_mesh->vertices_begin(); v_it!=_mesh->vertices_end(); ++v_it)
    {

        if (!_mesh->property(is_feature, v_it)) continue;

        TriMesh::Point p, q, r1, r2;
        p = _mesh->point(v_it);
        double totalAngle = 0;
        TriMesh::Normal e_Ha;
        TriMesh::Normal volGradient;


        // ( .. for each Incoming halfedge .. )
        TriMesh::VertexIHalfedgeIter vih_it(*_mesh,v_it);
        for ( ; vih_it; ++vih_it )
        {

            double sumAngle = 0;
            TriMesh::VertexHandle qvh = _mesh->from_vertex_handle(vih_it);
            q = _mesh->point(qvh);

            if (!_mesh->is_boundary(vih_it))
            {

                r1 = _mesh->point(_mesh->to_vertex_handle(_mesh->next_halfedge_handle(vih_it)));
                sumAngle += cal_angle(p, q, r1);

            }

            TriMesh::HalfedgeHandle opp = _mesh->opposite_halfedge_handle(vih_it);

            if (!_mesh->is_boundary(opp))
            {
                r2 = _mesh->point(_mesh->from_vertex_handle(_mesh->prev_halfedge_handle(opp)));
                sumAngle += cal_angle(p, q, r2);

            }

            totalAngle += sumAngle;

            _mesh->property(smoothed_amc_Ha, v_it) += sumAngle*_mesh->property(amc_Ha, qvh);

        }

        e_Ha = 0.5*(_mesh->property(amc_Ha, v_it) + _mesh->property(smoothed_amc_Ha, v_it)/totalAngle);
        e_Ha.normalize();

        volGradient = _mesh->property(volume_grad, v_it);
        double sign = volGradient|e_Ha;
        if (sign < 0) e_Ha *= -1;

        _mesh->property(volume_grad, v_it) = e_Ha;

    }
}





double PrescribedMeanCurvature::
cal_angle(const TriMesh::Point & p, const TriMesh::Point & q, const TriMesh::Point & r)
{
    TriMesh::Normal pq = q - p;
    TriMesh::Normal pr = r - p;
    pq.normalize();
    pr.normalize();

    return acos((pq|pr));
}


double PrescribedMeanCurvature::
edge_mean_curvature_He_Ne(TriMesh *_mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal)
{

    //dihedral = 0 on the boundary
    double dihedral = 0;

    TriMesh::HalfedgeHandle  hh1, hh2;
    TriMesh::FaceHandle    fh1, fh2;

    hh1 = _mesh->halfedge_handle(_eh, 0);
    hh2 = _mesh->halfedge_handle(_eh, 1);

    //normal is the edge vector rotated by 90 degree
    if ( _mesh->is_boundary( hh1 ) )
    {

        fh2 = _mesh->face_handle(hh2);
        TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);
        TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh2)) - _mesh->point(_mesh->from_vertex_handle(hh2));
        normal = (edgeVector%n2);
        normal.normalize();

    }else if ( _mesh->is_boundary( hh2 ) )
    {

        fh1 = _mesh->face_handle(hh1);
        TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
        TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh1)) - _mesh->point(_mesh->from_vertex_handle(hh1));
        normal = (edgeVector%n1);
        normal.normalize();

    }else
    {

        fh1 = _mesh->face_handle(hh1);
        fh2 = _mesh->face_handle(hh2);

        TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
        TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);

        normal = (n1+n2);
        normal.normalize();

        dihedral = _mesh->calc_dihedral_angle(_eh);

        //std::cout << "dihedral: " << dihedral << std::endl;

        if (dihedral < 0) normal *= -1;

    }

    double edgeLength = _mesh->calc_edge_length(_eh);
    return 2*edgeLength*cos(dihedral/2.0);

}

double PrescribedMeanCurvature::anisotropic_weight(double curvature, double lambda, double r)
{

    double weight;

    //printf("lambda %f curvature %f \n", lambda, curvature);

    if (fabs(curvature) <= lambda)
    {
        weight = 1;
    }
    else
    {
        weight = (lambda*lambda)/(r*pow(lambda-fabs(curvature), 2) + lambda*lambda);
    }

    return weight;
}

TriMesh::Scalar PrescribedMeanCurvature::face_area(TriMesh *_mesh, TriMesh::FaceHandle fh)
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

    return area;
}


/*
 * circulate around the oriented face using FaceHalfEdgeIter to extract p, q, r in CCW
 *
 * take (qxr)/6
 */
void
PrescribedMeanCurvature::
volume_gradient(TriMesh *_mesh, TriMesh::FaceHandle fh, TriMesh::VertexHandle vh, TriMesh::Normal & gradient)
{

    TriMesh::Point  q, r;
    TriMesh::VertexHandle vq, vr;

    //from p to q to r
    for (TriMesh::FaceHalfedgeIter fhh_it=_mesh->fh_iter(fh); fhh_it; ++fhh_it)
    {
        if (_mesh->from_vertex_handle(fhh_it) == vh)
        {
            vq = _mesh->to_vertex_handle(fhh_it);
            vr = _mesh->to_vertex_handle(_mesh->next_halfedge_handle(fhh_it));
        }
    }

    q = _mesh->point(vq);
    r = _mesh->point(vr);

    gradient = q%r;
    gradient /= 6;

}


void
PrescribedMeanCurvature::
updateLineNode(TriMeshObject * _meshObject
               , const OpenMesh::VPropHandleT< TriMesh::Normal > & anisoMeanCurvature
               , const OpenMesh::VPropHandleT< double >& areaStar)
{
    ACG::SceneGraph::LineNode * node = getLineNode(_meshObject);
    //OpenMesh::VPropHandleT< TriMesh::Normal > smoothVector;

    node->clear();

    for (TriMesh::VertexIter vit = _meshObject->mesh()->vertices_begin();
         vit != _meshObject->mesh()->vertices_end(); ++vit)
    {
        TriMesh::Point  p = _meshObject->mesh()->point(vit);
        TriMesh::Normal n = _meshObject->mesh()->property(anisoMeanCurvature, vit);
        TriMesh::Scalar coefficient = _meshObject->mesh()->property(areaStar, vit);
        coefficient = 3*TIME_STEP/coefficient;
        addLine(node, p, p+coefficient*n*1000, Color(255,0,0) );
    }
}


void
PrescribedMeanCurvature::
updateLineNode(TriMeshObject * _meshObject
               , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex)
{
    ACG::SceneGraph::LineNode * node = getLineNode(_meshObject);

    node->clear();

    double energy_error = 0;

    for (TriMesh::VertexIter vit = _meshObject->mesh()->vertices_begin();
         vit != _meshObject->mesh()->vertices_end(); ++vit)
    {
        TriMesh::Point  p = _meshObject->mesh()->point(vit);
        TriMesh::Point p_old = _meshObject->mesh()->property(old_vertex, vit);

        TriMesh::Normal n = p_old - p;

        energy_error += n.norm();

        addLine(node, p, p_old + n*50, Color(255,0,0) );
    }

    printf("energy error: %f \n", energy_error);

}


ACG::SceneGraph::LineNode *
PrescribedMeanCurvature::
getLineNode(TriMeshObject * _meshObject)
{
    ACG::SceneGraph::LineNode * line_node = 0;

    // get or add line node
    if( !_meshObject->hasAdditionalNode( "NormalEstimationPlugin", "LineNode" ) )
    {
        line_node = new ACG::SceneGraph::LineNode( ACG::SceneGraph::LineNode::LineSegmentsMode, _meshObject->manipulatorNode() );

        if( !_meshObject->addAdditionalNode(line_node, QString("NormalEstimationPlugin"), QString("LineNode") ) )
        {
            std::cerr << "NormalEstimationPlugin::getLineNode(): could not add line node.\n";
            return 0;
        }
        line_node->clear();
        line_node->set_line_width(1.0);
        line_node->set_base_color( ACG::Vec4f(255,0,255,255) );
        line_node->show();
    }
    else
        _meshObject->getAdditionalNode ( line_node, "NormalEstimationPlugin", "LineNode" );

    return line_node;
}

void
PrescribedMeanCurvature::
addLine( ACG::SceneGraph::LineNode * _line_node, Vec3d _p0, Vec3d _p1, Color _color )
{
    _line_node->add_line( _p0, _p1 );
    _line_node->add_color(_color);
}





void
PrescribedMeanCurvature::
clearLineNode(TriMeshObject * _meshObject)
{
    ACG::SceneGraph::LineNode * line_node = 0;

    // get line node
    if( _meshObject->hasAdditionalNode( "NormalEstimationPlugin", "LineNode" ) )
    {
        _meshObject->getAdditionalNode ( line_node, "NormalEstimationPlugin", "LineNode" );

        //line_node->clear();
//        line_node->show();

        line_node->hide();

    }

}




void
PrescribedMeanCurvature::
showLineNode(TriMeshObject * _meshObject)
{
    ACG::SceneGraph::LineNode * line_node = 0;

    // get line node
    if( _meshObject->hasAdditionalNode( "NormalEstimationPlugin", "LineNode" ) )
    {
        _meshObject->getAdditionalNode ( line_node, "NormalEstimationPlugin", "LineNode" );


        line_node->show();



    }

}



