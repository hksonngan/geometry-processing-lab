#include "implicit_integration.h"

Implicit_Integration::Implicit_Integration()
{
}


void Implicit_Integration::
init_vertex_vector(TriMesh *mesh, Eigen::VectorXd &vertices)
{
    unsigned int count = 0;
    for (TriMesh::VertexIter v_it=mesh->vertices_begin(); v_it!=mesh->vertices_end(); ++v_it)
    {
        TriMesh::Point point = mesh->point(v_it);
        vertices[count] = point[0];
        vertices[count+1] = point[1];
        vertices[count+2] = point[2];
        count += 3;
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

}








