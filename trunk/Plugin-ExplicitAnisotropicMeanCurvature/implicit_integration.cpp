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
                 , Eigen::SparseMatrix<double> & mass)
{

    Eigen::DynamicSparseMatrix<double> mass_dyn(mass.rows(), mass.cols());


}
