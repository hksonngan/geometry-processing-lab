#include "implicit_integration.h"

Implicit_Integration::Implicit_Integration()
{
}


void Implicit_Integration::
init_vertex_vector(TriMesh *mesh, Eigen::VectorXd &vertices)
{
    unsigned int id = 0;
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        TriMesh::Point point = mesh->point(v_it);
        vertices[id*3] = point[0];
        vertices[id*3+1] = point[1];
        vertices[id*3+2] = point[2];
        id++;
    }
}


void Implicit_Integration::
init_mass_matrix(TriMesh *mesh , PrescribedMeanCurvature * pmc
                 , const OpenMesh::VPropHandleT< double > & area_star
                 , const OpenMesh::VPropHandleT< int > & vertex_id
                 , Eigen::SparseMatrix<double> & mass)
{

    Eigen::DynamicSparseMatrix<double> mass_dyn(mass.rows(), mass.cols());

    for (TriMesh::EdgeIter e_it=mesh->edges_begin(); e_it!=mesh->edges_end(); ++e_it)
    {
        double edge_area_star = (1.0/12.0)*pmc->area_star_edge(mesh, e_it.handle());
        TriMesh::HalfedgeHandle hh = mesh->halfedge_handle(e_it.handle(), 0);

        TriMesh::VertexHandle p, q;

        p = mesh->from_vertex_handle(hh);
        q = mesh->to_vertex_handle(hh);

        int pId = mesh->property(vertex_id, p);
        int qId = mesh->property(vertex_id, q);

        //it's safe to have assignment operator here
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
init_amc_matrix(TriMesh *mesh, PrescribedMeanCurvature *pmc
                , const OpenMesh::VPropHandleT< int > &vertex_id
                , Eigen::SparseMatrix<double> &amc_matrix)
{

    Eigen::DynamicSparseMatrix<double> amc_dyn(amc_matrix.rows(), amc_matrix.cols());

    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {

        for (TriMesh::VertexOHalfedgeIter voh_it=mesh->voh_iter(v_it.handle()); voh_it; ++voh_it)
        {

            TriMesh::EdgeHandle eh = mesh->edge_handle(voh_it.handle());
            Mat3x3 the_cross;
            TriMesh::Scalar normalLength;
            TriMesh::Scalar meanCurvature = pmc->calculate_cross_matrix_Ax_qjpi(mesh, eh
                                                                           , normalLength
                                                                           , v_it.handle(), the_cross);


            double weight = pmc->anisotropic_weight(meanCurvature, pmc->get_lambda()
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
compute_semi_implicit_integration(TriMesh *mesh
                                  , unsigned int mesh_size
                                  , const OpenMesh::VPropHandleT< double > & area_star
                                  , const OpenMesh::VPropHandleT< int > & vertex_id
                                  , const OpenMesh::VPropHandleT< TriMesh::Point > & old_vertex)
{

    printf("entering implicit step \n");

    PrescribedMeanCurvature pmc;
    mesh_size *= 3;

    Eigen::VectorXd input_vertices(mesh_size);
    this->init_vertex_vector(mesh, input_vertices);
    printf("done init vertex vector \n");

    Eigen::SparseMatrix<double> mass_matrix(mesh_size, mesh_size);
    this->init_mass_matrix(mesh, &pmc, area_star, vertex_id, mass_matrix);
    printf("done init mass matrix \n");

    Eigen::SparseMatrix<double> amc_matrix(mesh_size, mesh_size);
    this->init_amc_matrix(mesh, &pmc, vertex_id, amc_matrix);
    printf("done init amc matrix \n");

    Eigen::SparseMatrix<double> A(mesh_size, mesh_size);
    A = mass_matrix + IMPLICIT_TIME_FACTOR*PrescribedMeanCurvature::TIME_STEP*amc_matrix;

    Eigen::VectorXd b(mesh_size);
    b = mass_matrix*input_vertices;

    Eigen::VectorXd x(mesh_size);
    Eigen::SparseLLT<Eigen::SparseMatrix<double>, Eigen::Cholmod> cholmoDec;
    printf("start computing the matrix \n");
    cholmoDec.compute(A);
    printf("start solving the equations \n");
    x = cholmoDec.solve(b);
    printf("start updating the mesh \n");

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














