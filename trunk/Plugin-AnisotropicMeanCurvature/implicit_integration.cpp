#include "implicit_integration.h"



Implicit_Integration::Implicit_Integration()
{
}



void Implicit_Integration::
init_vertex_vector(TriMesh *mesh, Eigen::VectorXd &vertices, const OpenMesh::VPropHandleT< int > & vertex_id)
{
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        unsigned int id = mesh->property(vertex_id, v_it);
        TriMesh::Point point = mesh->point(v_it);
        vertices[id*3] = point[0];
        vertices[id*3+1] = point[1];
        vertices[id*3+2] = point[2];
    }
}



void Implicit_Integration::
init_mass_matrix(TriMesh *mesh
                 , const OpenMesh::VPropHandleT< double > & area_star
                 , const OpenMesh::VPropHandleT< int > & vertex_id
                 , Eigen::SparseMatrix<double> & mass)
{
    Eigen::DynamicSparseMatrix<double> mass_dyn(mass.rows(), mass.cols());

    for (TriMesh::EdgeIter e_it=mesh->edges_begin(); e_it!=mesh->edges_end(); ++e_it)
    {
        double edge_area_star = (1.0/12.0)*area_star_edge(mesh, e_it.handle());
        TriMesh::HalfedgeHandle hh = mesh->halfedge_handle(e_it.handle(), 0);

        TriMesh::VertexHandle p, q;

        p = mesh->from_vertex_handle(hh);
        q = mesh->to_vertex_handle(hh);

        int pId = mesh->property(vertex_id, p);
        int qId = mesh->property(vertex_id, q);

        //it's safe to have assignment operator here
        //each row corresponds to one vertex
        //there are at most n-1 edges emanating from one vertex
        mass_dyn.coeffRef(pId*3, qId*3) = edge_area_star;
        mass_dyn.coeffRef(qId*3, pId*3) = edge_area_star;

        mass_dyn.coeffRef(pId*3+1, qId*3+1) = edge_area_star;
        mass_dyn.coeffRef(qId*3+1, pId*3+1) = edge_area_star;

        mass_dyn.coeffRef(pId*3+2, qId*3+2) = edge_area_star;
        mass_dyn.coeffRef(qId*3+2, pId*3+2) = edge_area_star;
    }

    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        double vertex_area_star = (1.0/6.0)*mesh->property(area_star, v_it);

        int pId = mesh->property(vertex_id, v_it.handle());

        mass_dyn.coeffRef(pId*3, pId*3) = vertex_area_star;
        mass_dyn.coeffRef(pId*3+1, pId*3+1) = vertex_area_star;
        mass_dyn.coeffRef(pId*3+2, pId*3+2) = vertex_area_star;
    }

    mass = Eigen::SparseMatrix<double>(mass_dyn);
}



void Implicit_Integration::
init_lumped_mass_matrix(TriMesh *mesh
                 , const OpenMesh::VPropHandleT< double > & area_star
                 , const OpenMesh::VPropHandleT< int > & vertex_id
                 , Eigen::SparseMatrix<double> & mass)
{
    Eigen::DynamicSparseMatrix<double> mass_dyn(mass.rows(), mass.cols());

    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {

        double vertex_area_star = (1.0/3.0)*mesh->property(area_star, v_it);

        //printf("area star %f \n", vertex_area_star);

        int pId = mesh->property(vertex_id, v_it.handle());

        mass_dyn.coeffRef(pId*3, pId*3) = vertex_area_star;
        mass_dyn.coeffRef(pId*3+1, pId*3+1) = vertex_area_star;
        mass_dyn.coeffRef(pId*3+2, pId*3+2) = vertex_area_star;
    }

    mass = Eigen::SparseMatrix<double>(mass_dyn);
}



void Implicit_Integration::
init_lumped_mass_matrix_inverted(TriMesh *mesh
                 , const OpenMesh::VPropHandleT< double > & area_star
                 , const OpenMesh::VPropHandleT< int > & vertex_id
                 , Eigen::SparseMatrix<double> & mass)
{
    Eigen::DynamicSparseMatrix<double> mass_dyn(mass.rows(), mass.cols());

    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        double vertex_area_star = (3.0)/mesh->property(area_star, v_it);

        //printf("area star %f \n", vertex_area_star);

        int pId = mesh->property(vertex_id, v_it.handle());

        mass_dyn.coeffRef(pId*3, pId*3) = vertex_area_star;
        mass_dyn.coeffRef(pId*3+1, pId*3+1) = vertex_area_star;
        mass_dyn.coeffRef(pId*3+2, pId*3+2) = vertex_area_star;
    }

    mass = Eigen::SparseMatrix<double>(mass_dyn);
}



void Implicit_Integration::
init_amc_matrix(TriMesh *mesh
                , const OpenMesh::VPropHandleT< int > &vertex_id
                , Eigen::SparseMatrix<double> &amc_matrix)
{
    double threshold = get_feature_threshold(mesh);

    Eigen::DynamicSparseMatrix<double> amc_dyn(amc_matrix.rows(), amc_matrix.cols());

    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        for (TriMesh::VertexOHalfedgeIter voh_it=mesh->voh_iter(v_it.handle()); voh_it; ++voh_it)
        {
            TriMesh::EdgeHandle eh = mesh->edge_handle(voh_it.handle());
            Mat3x3 the_cross;
            TriMesh::Scalar normalLength;
            TriMesh::Scalar meanCurvature = calculate_cross_matrix_Ax_qjpi(mesh, eh, normalLength
                                                                           , v_it.handle(), the_cross);

            double weight = anisotropic_weight(meanCurvature, threshold
                                                    , PrescribedMeanCurvature::R);

            the_cross *= ((0.5*meanCurvature*weight)/normalLength);

            //identify the vertex qj+1 and qj-1
            //be careful about the haldfedge directions or the normal will point inward
            TriMesh::VertexHandle q_plus, q_minus;

            q_plus = mesh->to_vertex_handle(mesh->next_halfedge_handle(voh_it));
            q_minus = mesh->from_vertex_handle(
                        mesh->prev_halfedge_handle( mesh->opposite_halfedge_handle(voh_it) ) );

            int q_plus_id = mesh->property(vertex_id, q_plus);
            int q_minus_id = mesh->property(vertex_id, q_minus);
            int p_id = mesh->property(vertex_id, v_it);

            //put the_cross in the q_plus_id and -the_cross in the q_minus_id positions
            //remember to sum up the value at each entry

            amc_dyn.coeffRef(p_id*3, q_plus_id*3) += the_cross(0, 0);
            amc_dyn.coeffRef(p_id*3, q_plus_id*3+1) += the_cross(0, 1);
            amc_dyn.coeffRef(p_id*3, q_plus_id*3+2) += the_cross(0, 2);
            amc_dyn.coeffRef(p_id*3+1, q_plus_id*3) += the_cross(1, 0);
            amc_dyn.coeffRef(p_id*3+1, q_plus_id*3+1) += the_cross(1, 1);
            amc_dyn.coeffRef(p_id*3+1, q_plus_id*3+2) += the_cross(1, 2);
            amc_dyn.coeffRef(p_id*3+2, q_plus_id*3) += the_cross(2, 0);
            amc_dyn.coeffRef(p_id*3+2, q_plus_id*3+1) += the_cross(2, 1);
            amc_dyn.coeffRef(p_id*3+2, q_plus_id*3+2) += the_cross(2, 2);

            amc_dyn.coeffRef(p_id*3, q_minus_id*3) -= the_cross(0, 0);
            amc_dyn.coeffRef(p_id*3, q_minus_id*3+1) -= the_cross(0, 1);
            amc_dyn.coeffRef(p_id*3, q_minus_id*3+2) -= the_cross(0, 2);
            amc_dyn.coeffRef(p_id*3+1, q_minus_id*3) -= the_cross(1, 0);
            amc_dyn.coeffRef(p_id*3+1, q_minus_id*3+1) -= the_cross(1, 1);
            amc_dyn.coeffRef(p_id*3+1, q_minus_id*3+2) -= the_cross(1, 2);
            amc_dyn.coeffRef(p_id*3+2, q_minus_id*3) -= the_cross(2, 0);
            amc_dyn.coeffRef(p_id*3+2, q_minus_id*3+1) -= the_cross(2, 1);
            amc_dyn.coeffRef(p_id*3+2, q_minus_id*3+2) -= the_cross(2, 2);
        }
    }

    amc_matrix = Eigen::SparseMatrix<double>(amc_dyn);
}



void Implicit_Integration::
init_Jacobian(TriMesh *mesh
                , const OpenMesh::VPropHandleT< int > &vertex_id
                , Eigen::SparseMatrix<double> &jacobian)
{
    double threshold = get_feature_threshold(mesh);

    Eigen::DynamicSparseMatrix<double> jacobian_dyn(jacobian.rows(), jacobian.cols());

    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        for (TriMesh::VertexOHalfedgeIter voh_it=mesh->voh_iter(v_it.handle()); voh_it; ++voh_it)
        {
            TriMesh::Scalar normalLength;
            std::vector<Mat3x3> derivatives;
            TriMesh::Scalar meanCurvature = cal_edge_norm_deriv_qjpi(mesh, voh_it.handle()
                                                                     , normalLength, derivatives);

            double weight = anisotropic_weight(meanCurvature, threshold
                                                    , PrescribedMeanCurvature::R);

            //the_cross *= ((0.5*meanCurvature*weight)/normalLength);
            weight *= 0.5*meanCurvature/normalLength;

            //identify the vertex qj+1 and qj-1
            //be careful about the haldfedge directions or the normal will point inward
            TriMesh::VertexHandle q_plus, q_minus, qj;

            qj = mesh->to_vertex_handle(voh_it);
            q_plus = mesh->to_vertex_handle(mesh->next_halfedge_handle(voh_it));
            q_minus = mesh->from_vertex_handle(
                        mesh->prev_halfedge_handle( mesh->opposite_halfedge_handle(voh_it) ) );

            int q_plus_id = mesh->property(vertex_id, q_plus);
            int q_minus_id = mesh->property(vertex_id, q_minus);
            int p_id = mesh->property(vertex_id, v_it);
            int qj_id = mesh->property(vertex_id, qj);

            //put derivatives[0] at pi, derivatives[1] at q_minus, derivatives[2] at qj, derivatives[3] at q_plus
            //remember to multiply them by the common weight and to accumulate the value at each entry

            Mat3x3 the_cross;

            the_cross = derivatives[0]*weight;
            jacobian_dyn.coeffRef(p_id*3, p_id*3) += the_cross(0, 0);
            jacobian_dyn.coeffRef(p_id*3, p_id*3+1) += the_cross(0, 1);
            jacobian_dyn.coeffRef(p_id*3, p_id*3+2) += the_cross(0, 2);
            jacobian_dyn.coeffRef(p_id*3+1, p_id*3) += the_cross(1, 0);
            jacobian_dyn.coeffRef(p_id*3+1, p_id*3+1) += the_cross(1, 1);
            jacobian_dyn.coeffRef(p_id*3+1, p_id*3+2) += the_cross(1, 2);
            jacobian_dyn.coeffRef(p_id*3+2, p_id*3) += the_cross(2, 0);
            jacobian_dyn.coeffRef(p_id*3+2, p_id*3+1) += the_cross(2, 1);
            jacobian_dyn.coeffRef(p_id*3+2, p_id*3+2) += the_cross(2, 2);

            the_cross = derivatives[1]*weight;
            jacobian_dyn.coeffRef(p_id*3, q_minus_id*3) += the_cross(0, 0);
            jacobian_dyn.coeffRef(p_id*3, q_minus_id*3+1) += the_cross(0, 1);
            jacobian_dyn.coeffRef(p_id*3, q_minus_id*3+2) += the_cross(0, 2);
            jacobian_dyn.coeffRef(p_id*3+1, q_minus_id*3) += the_cross(1, 0);
            jacobian_dyn.coeffRef(p_id*3+1, q_minus_id*3+1) += the_cross(1, 1);
            jacobian_dyn.coeffRef(p_id*3+1, q_minus_id*3+2) += the_cross(1, 2);
            jacobian_dyn.coeffRef(p_id*3+2, q_minus_id*3) += the_cross(2, 0);
            jacobian_dyn.coeffRef(p_id*3+2, q_minus_id*3+1) += the_cross(2, 1);
            jacobian_dyn.coeffRef(p_id*3+2, q_minus_id*3+2) += the_cross(2, 2);

            the_cross = derivatives[2]*weight;
            jacobian_dyn.coeffRef(p_id*3, qj_id*3) += the_cross(0, 0);
            jacobian_dyn.coeffRef(p_id*3, qj_id*3+1) += the_cross(0, 1);
            jacobian_dyn.coeffRef(p_id*3, qj_id*3+2) += the_cross(0, 2);
            jacobian_dyn.coeffRef(p_id*3+1, qj_id*3) += the_cross(1, 0);
            jacobian_dyn.coeffRef(p_id*3+1, qj_id*3+1) += the_cross(1, 1);
            jacobian_dyn.coeffRef(p_id*3+1, qj_id*3+2) += the_cross(1, 2);
            jacobian_dyn.coeffRef(p_id*3+2, qj_id*3) += the_cross(2, 0);
            jacobian_dyn.coeffRef(p_id*3+2, qj_id*3+1) += the_cross(2, 1);
            jacobian_dyn.coeffRef(p_id*3+2, qj_id*3+2) += the_cross(2, 2);

            the_cross = derivatives[3]*weight;
            jacobian_dyn.coeffRef(p_id*3, q_plus_id*3) += the_cross(0, 0);
            jacobian_dyn.coeffRef(p_id*3, q_plus_id*3+1) += the_cross(0, 1);
            jacobian_dyn.coeffRef(p_id*3, q_plus_id*3+2) += the_cross(0, 2);
            jacobian_dyn.coeffRef(p_id*3+1, q_plus_id*3) += the_cross(1, 0);
            jacobian_dyn.coeffRef(p_id*3+1, q_plus_id*3+1) += the_cross(1, 1);
            jacobian_dyn.coeffRef(p_id*3+1, q_plus_id*3+2) += the_cross(1, 2);
            jacobian_dyn.coeffRef(p_id*3+2, q_plus_id*3) += the_cross(2, 0);
            jacobian_dyn.coeffRef(p_id*3+2, q_plus_id*3+1) += the_cross(2, 1);
            jacobian_dyn.coeffRef(p_id*3+2, q_plus_id*3+2) += the_cross(2, 2);
        }
    }

    jacobian = Eigen::SparseMatrix<double>(jacobian_dyn);
}



void Implicit_Integration::
compute_taylor_semi_implicit(TriMesh *mesh
                             , unsigned int mesh_size
                             , const OpenMesh::VPropHandleT< double > & area_star
                             , const OpenMesh::VPropHandleT< int > & vertex_id
                             , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex
                             , double time_step
                             , bool is_lumped_mass
                             )
{
    mesh_size *= 3;

    Eigen::VectorXd input_vertices(mesh_size);
    this->init_vertex_vector(mesh, input_vertices, vertex_id);

    Eigen::SparseMatrix<double> mass_matrix(mesh_size, mesh_size);

    if (!is_lumped_mass)
    {
        this->init_mass_matrix(mesh, area_star, vertex_id, mass_matrix);
    }else
    {
        this->init_lumped_mass_matrix(mesh, area_star, vertex_id, mass_matrix);
    }

    Eigen::SparseMatrix<double> amc_matrix(mesh_size, mesh_size);
    this->init_amc_matrix(mesh, vertex_id, amc_matrix);

    Eigen::SparseMatrix<double> jacobian(mesh_size, mesh_size);
    this->init_Jacobian(mesh, vertex_id, jacobian);

    Eigen::SparseMatrix<double> A(mesh_size, mesh_size);
    A = mass_matrix + time_step*jacobian;

    Eigen::VectorXd b(mesh_size);
    b = mass_matrix*input_vertices - time_step*amc_matrix*input_vertices
            + time_step*jacobian*input_vertices;

    Eigen::VectorXd x(mesh_size);
    Eigen::SparseLLT<Eigen::SparseMatrix<double>, Eigen::Cholmod> cholmoDec;
    cholmoDec.compute(A);
    x = cholmoDec.solve(b);

    //std::cout << "residual: " << (A * x - b).norm() << std::endl;

    //convert the result back to vertex and update the mesh
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {

        TriMesh::Point old_point = mesh->point(v_it);
        mesh->property(old_vertex, v_it) = old_point;

        TriMesh::Point point = TriMesh::Point(0,0,0);
        int idx = mesh->property(vertex_id, v_it);

        point[0] = x[idx*3];
        point[1] = x[idx*3+1];
        point[2] = x[idx*3+2];

        mesh->set_point(v_it, point);
    }
}



void Implicit_Integration::
smooth_aniso(int _iterations, TriMeshObject * meshObject
            , SmoothingMode smooth_type
            , IntegrationScheme scheme
            , VisualizeMode visualize
            , double time_step)
{
    TriMesh* mesh = meshObject->mesh();

    // Property for the active mesh to store original point positions
    OpenMesh::VPropHandleT< TriMesh::Point > old_vertex;
    OpenMesh::VPropHandleT< double > area_star;
    OpenMesh::VPropHandleT< int > vertex_id;

    // Add a property to the mesh to store mean curvature and area
    mesh->add_property( old_vertex, "old_amc_vertex" );
    mesh->add_property( area_star, "area_star" );
    mesh->add_property( vertex_id, "vertex_id" );

    // Property for the active mesh to store original point positions
    OpenMesh::VPropHandleT< TriMesh::Normal > amc_Ha;
    //the smoothed_amc_Ha is for temporary use
    OpenMesh::VPropHandleT< TriMesh::Normal > smoothed_amc_Ha;
    OpenMesh::VPropHandleT< double > apmc_function_f;
    OpenMesh::VPropHandleT< double > smoothed_apmc_function_f;
    OpenMesh::VPropHandleT< bool > is_feature;
    OpenMesh::VPropHandleT< TriMesh::Normal > volume_gradient_Va;

    // Add a property to the mesh to store mean curvature and area
    mesh->add_property( amc_Ha, "anisotropic_mean_curvature_Ha" );
    mesh->add_property( smoothed_amc_Ha, "smoothed_anisotropic_mean_curvature_Ha" );
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

            mesh->property(amc_Ha,v_it).vectorize(0.0f);
            mesh->property(smoothed_amc_Ha,v_it).vectorize(0.0f);
            mesh->property(volume_gradient_Va,v_it).vectorize(0.0f);
            mesh->property(apmc_function_f,v_it) = 0;
            mesh->property(smoothed_apmc_function_f,v_it) = 0;
            mesh->property(is_feature,v_it) = false;
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
                //noFeatureVertices++;
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

        if (scheme == BATCH_UPDATE
                && smooth_type == ANISO_MEAN_CURVATURE)
        {
            compute_explicit_integration_with_mass(mesh, count, area_star, vertex_id, old_vertex, time_step);
        }

        if(scheme == BATCH_UPDATE && smooth_type == PRESCRIBED_MEAN_CURVATURE)
        {
            compute_pmc(mesh, time_step);
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

    mesh->remove_property( amc_Ha );
    mesh->remove_property( is_feature );
    mesh->remove_property( volume_gradient_Va );
    mesh->remove_property( smoothed_amc_Ha );
    mesh->remove_property( smoothed_apmc_function_f );
    mesh->remove_property( apmc_function_f );
}



void Implicit_Integration::
compute_semi_implicit_integration(TriMesh *mesh
                                  , unsigned int mesh_size
                                  , const OpenMesh::VPropHandleT< double > & area_star
                                  , const OpenMesh::VPropHandleT< int > & vertex_id
                                  , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex
                                  , double time_step
                                  , bool is_lumped_mass
                                  )
{
    mesh_size *= 3;

    Eigen::VectorXd input_vertices(mesh_size);
    this->init_vertex_vector(mesh, input_vertices, vertex_id);

    Eigen::SparseMatrix<double> mass_matrix(mesh_size, mesh_size);

    if (!is_lumped_mass)
    {
        this->init_mass_matrix(mesh, area_star, vertex_id, mass_matrix);
    }else
    {
        this->init_lumped_mass_matrix(mesh, area_star, vertex_id, mass_matrix);
    }

    Eigen::SparseMatrix<double> amc_matrix(mesh_size, mesh_size);
    this->init_amc_matrix(mesh, vertex_id, amc_matrix);

    Eigen::SparseMatrix<double> A(mesh_size, mesh_size);
    A = mass_matrix + time_step*amc_matrix;

    Eigen::VectorXd b(mesh_size);
    b = mass_matrix*input_vertices;

    Eigen::VectorXd x(mesh_size);
    Eigen::SparseLLT<Eigen::SparseMatrix<double>, Eigen::Cholmod> cholmoDec;
    cholmoDec.compute(A);
    x = cholmoDec.solve(b);

    //std::cout << "residual: " << (A * x - b).norm() << std::endl;

    //convert the result back to vertex and update the mesh
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {

        TriMesh::Point old_point = mesh->point(v_it);
        mesh->property(old_vertex, v_it) = old_point;

        if(mesh->is_boundary(v_it.handle())) {
            continue;
        }

        TriMesh::Point point = TriMesh::Point(0,0,0);
        int idx = mesh->property(vertex_id, v_it);

        point[0] = x[idx*3];
        point[1] = x[idx*3+1];
        point[2] = x[idx*3+2];

        mesh->set_point(v_it, point);
    }
}



void Implicit_Integration::
compute_explicit_integration_with_mass(TriMesh *mesh
                                       , unsigned int mesh_size
                                       , const OpenMesh::VPropHandleT< double > & area_star
                                       , const OpenMesh::VPropHandleT< int > & vertex_id
                                       , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex
                                       , double time_step
                                       , bool is_lumped_mask
                                       )
{
    mesh_size *= 3;

    Eigen::VectorXd x(mesh_size);
    Eigen::VectorXd b(mesh_size);

    Eigen::VectorXd input_vertices(mesh_size);
    this->init_vertex_vector(mesh, input_vertices, vertex_id);

    Eigen::SparseMatrix<double> amc_matrix(mesh_size, mesh_size);
    this->init_amc_matrix(mesh, vertex_id, amc_matrix);

    Eigen::VectorXd Ha(mesh_size);
    Ha = amc_matrix*input_vertices;

    //Eigen::SparseMatrix<double> mass_matrix_inverted(mesh_size, mesh_size);
    Eigen::SparseMatrix<double> mass_matrix(mesh_size, mesh_size);

    if (!is_lumped_mask)
    {
        this->init_mass_matrix(mesh, area_star, vertex_id, mass_matrix);
        Eigen::SparseLLT<Eigen::SparseMatrix<double>, Eigen::Cholmod> cholmoDec;
        cholmoDec.compute(mass_matrix);
        x = cholmoDec.solve(Ha);
        b = time_step*x;
    }else
    {
        this->init_lumped_mass_matrix_inverted(mesh, area_star, vertex_id, mass_matrix);
        b = time_step*mass_matrix*Ha;
    }

    //convert the result back to vertex and update the mesh
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        TriMesh::Point old_point = mesh->point(v_it);
        mesh->property(old_vertex, v_it) = old_point;

        if(mesh->is_boundary(v_it.handle())) {
            continue;
        }

        TriMesh::Point point = TriMesh::Point(0,0,0);
        int idx = mesh->property(vertex_id, v_it);

        point[0] = old_point[0] - b[idx*3];
        point[1] = old_point[1] - b[idx*3+1];
        point[2] = old_point[2] - b[idx*3+2];

        mesh->set_point(v_it, point);
    }
}



void Implicit_Integration::init_pmc(TriMesh * mesh
              , const OpenMesh::VPropHandleT< int > & vertex_id
              , const OpenMesh::VPropHandleT< double > & smoothed_apmc_function_f
              , const OpenMesh::VPropHandleT< TriMesh::Normal > & volume_gradient_Va
              , Eigen::VectorXd & pmc_vector)
{
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        unsigned int id = mesh->property(vertex_id, v_it);

        TriMesh::Normal volumeGrad = mesh->property(smoothed_apmc_function_f, v_it)
                *mesh->property(volume_gradient_Va, v_it);

        pmc_vector[id*3] = volumeGrad[0];
        pmc_vector[id*3+1] = volumeGrad[1];
        pmc_vector[id*3+2] = volumeGrad[2];
    }
}



void Implicit_Integration::
compute_pmc(TriMesh *mesh, double time_step)
{

    OpenMesh::VPropHandleT< double > area_star;
    OpenMesh::VPropHandleT< double > smoothed_apmc_function_f;
    OpenMesh::VPropHandleT< TriMesh::Normal > volume_gradient_Va;
    OpenMesh::VPropHandleT< int > vertex_id;
    OpenMesh::VPropHandleT< TriMesh::Point > old_vertex;

    unsigned int mesh_size = mesh->n_vertices();
    mesh->get_property_handle(area_star, "area_star");
    mesh->get_property_handle(smoothed_apmc_function_f, "smoothed_apmc_function_f");
    mesh->get_property_handle(volume_gradient_Va, "volume_gradient_Va");
    mesh->get_property_handle( vertex_id, "vertex_id" );
    mesh->get_property_handle( old_vertex, "old_amc_vertex" );

    mesh_size *= 3;

    Eigen::VectorXd x(mesh_size);
    Eigen::VectorXd b(mesh_size);

    Eigen::VectorXd input_vertices(mesh_size);
    this->init_vertex_vector(mesh, input_vertices, vertex_id);

    Eigen::SparseMatrix<double> amc_matrix(mesh_size, mesh_size);
    this->init_amc_matrix(mesh, vertex_id, amc_matrix);

    Eigen::VectorXd fVa(mesh_size);
    this->init_pmc(mesh, vertex_id, smoothed_apmc_function_f, volume_gradient_Va, fVa);

    Eigen::VectorXd Ha(mesh_size);
    Ha = (amc_matrix*input_vertices) - fVa*MODULATOR;

    Eigen::SparseMatrix<double> mass_matrix(mesh_size, mesh_size);

    this->init_mass_matrix(mesh, area_star, vertex_id, mass_matrix);
    Eigen::SparseLLT<Eigen::SparseMatrix<double>, Eigen::Cholmod> cholmoDec;
    cholmoDec.compute(mass_matrix);
    x = cholmoDec.solve(Ha);
    b = time_step*x;

    //convert the result back to vertex and update the mesh
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {

        TriMesh::Point old_point = mesh->point(v_it);
        mesh->property(old_vertex, v_it) = old_point;

        if(mesh->is_boundary(v_it.handle())) {
            continue;
        }

        TriMesh::Point point = TriMesh::Point(0,0,0);
        int idx = mesh->property(vertex_id, v_it);

        //TriMesh::Scalar area = mesh->property(area_star, v_it);

        point[0] = old_point[0] - b[idx*3];
        point[1] = old_point[1] - b[idx*3+1];
        point[2] = old_point[2] - b[idx*3+2];

        mesh->set_point(v_it, point);
    }
}


















