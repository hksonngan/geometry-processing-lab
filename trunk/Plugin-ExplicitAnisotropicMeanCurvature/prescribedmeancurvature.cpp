#include "prescribedmeancurvature.h"

PrescribedMeanCurvature::PrescribedMeanCurvature()
{


    double PrescribedMeanCurvature::edgeMeanCurvature(TriMesh *_mesh, TriMesh::EdgeHandle _eh, TriMesh::Normal & normal)
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

        }else if ( _mesh->is_boundary( hh2 ) )
        {

            fh1 = _mesh->face_handle(hh1);
            TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
            TriMesh::Normal edgeVector = _mesh->point(_mesh->to_vertex_handle(hh1)) - _mesh->point(_mesh->from_vertex_handle(hh1));
            normal = (edgeVector%n1);

        }else
        {

            fh1 = _mesh->face_handle(hh1);
            fh2 = _mesh->face_handle(hh2);

            TriMesh::Normal n1 = _mesh->calc_face_normal(fh1);
            TriMesh::Normal n2 = _mesh->calc_face_normal(fh2);

            normal = (n1+n2);
            normal.normalize();

            dihedral = _mesh->calc_dihedral_angle(_eh);

            if (dihedral < 0) normal *= -1;

        }

        double edgeLength = _mesh->calc_edge_length(_eh);
        return 2*edgeLength*cos(dihedral/2.0);

    }

    double PrescribedMeanCurvature::anisotropicWeight(double curvature, double lambda, double r)
    {

        double weight;

        //printf("lambda curvature %f %f \n", lambda, curvature);

        if (fabs(curvature) <= lambda)
        {
            weight = 1;
            //printf("weight %f \n", weight);
        }
        else
        {
            weight = (lambda*lambda)/(r*pow(lambda-fabs(curvature), 2) + lambda*lambda);
        }

        //printf("weight %f \n", weight);

        return weight;
    }

    TriMesh::Scalar PrescribedMeanCurvature::faceArea(TriMesh *_mesh, TriMesh::FaceHandle fh)
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


}
