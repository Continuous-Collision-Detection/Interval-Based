#include <iostream>


#include <interval_ccd/interval_ccd.hpp>
#include <interval_ccd/interval_root_finder.hpp>

#include <vector>

namespace intervalccd
{

double compute_face_vertex_tolerance_1d(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    const double input_tol)
{
    double dl = std::max(
        std::max(
            (vertex_end - vertex_start).norm(),
            (face_vertex0_end - face_vertex0_start).norm()),
        std::max(
            (face_vertex1_end - face_vertex1_start).norm(),
            (face_vertex2_end - face_vertex2_start).norm()));

    return input_tol / dl;
}

Eigen::Vector3I cross(const Eigen::Vector3I& a, const Eigen::Vector3I& b)
{
    Eigen::Vector3I c;
    c(0) = a(1) * b(2) - a(2) * b(1);
    c(1) = a(2) * b(0) - a(0) * b(2);
    c(2) = a(0) * b(1) - a(1) * b(0);
    return c;
}

bool are_points_on_same_side_of_edge(
    const Eigen::Vector3I& p1,
    const Eigen::Vector3I& p2,
    const Eigen::Vector3I& a,
    const Eigen::Vector3I& b)
{
    Eigen::Vector3I cp1 = cross(b - a, p1 - a); //(b - a).cross(p1 - a);
    Eigen::Vector3I cp2 = cross(b - a, p2 - a); //(b - a).cross(p2 - a);
    return cp1.dot(cp2).upper() >= 0;
}

bool is_point_inside_triangle_(
    const Eigen::Vector3I& point,
    const Eigen::Vector3I& triangle_vertex0,
    const Eigen::Vector3I& triangle_vertex1,
    const Eigen::Vector3I& triangle_vertex2)
{
    return are_points_on_same_side_of_edge(
               point, triangle_vertex0, triangle_vertex1, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex1, triangle_vertex0, triangle_vertex2)
        && are_points_on_same_side_of_edge(
               point, triangle_vertex2, triangle_vertex0, triangle_vertex1);
}

inline Eigen::Vector3I triangle_normal(
    const Eigen::Vector3I& face_vertex0,
    const Eigen::Vector3I& face_vertex1,
    const Eigen::Vector3I& face_vertex2)
{
    return cross(
        face_vertex1 - face_vertex0,
        face_vertex2
            - face_vertex0); //(face_vertex1 - face_vertex0).cross(face_vertex2
                             //- face_vertex0);
}

bool vertexFaceCCD_Redon(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    double& toi,
    const double &input_tol)
{
    const auto distance = [&](const Interval& t) {
        // Get the vertex at time t
        Eigen::Vector3I vertex
            = (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t
            + vertex_start.cast<Interval>();

        // Get the vertex of the face at time t
        Eigen::Vector3I face_vertex0 = (face_vertex0_end.cast<Interval>()
                                        - face_vertex0_start.cast<Interval>())
                * t
            + face_vertex0_start.cast<Interval>();
        Eigen::Vector3I face_vertex1 = (face_vertex1_end.cast<Interval>()
                                        - face_vertex1_start.cast<Interval>())
                * t
            + face_vertex1_start.cast<Interval>();
        Eigen::Vector3I face_vertex2 = (face_vertex2_end.cast<Interval>()
                                        - face_vertex2_start.cast<Interval>())
                * t
            + face_vertex2_start.cast<Interval>();

        return (vertex - face_vertex0)
            .dot(triangle_normal(face_vertex0, face_vertex1, face_vertex2));
    };

    const auto is_point_inside_triangle = [&](const Interval& t) {
        // Get the vertex at time t
        Eigen::Vector3I vertex
            = (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t
            + vertex_start.cast<Interval>();

        // Get the vertex of the face at time t
        Eigen::Vector3I face_vertex0 = (face_vertex0_end.cast<Interval>()
                                        - face_vertex0_start.cast<Interval>())
                * t
            + face_vertex0_start.cast<Interval>();
        Eigen::Vector3I face_vertex1 = (face_vertex1_end.cast<Interval>()
                                        - face_vertex1_start.cast<Interval>())
                * t
            + face_vertex1_start.cast<Interval>();
        Eigen::Vector3I face_vertex2 = (face_vertex2_end.cast<Interval>()
                                        - face_vertex2_start.cast<Interval>())
                * t
            + face_vertex2_start.cast<Interval>();

        return is_point_inside_triangle_(
            vertex, face_vertex0, face_vertex1, face_vertex2);
    };

    double tol = compute_face_vertex_tolerance_1d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end,input_tol);

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_point_inside_triangle, Interval(0, 1), tol, toi_interval);

    // Return a conservative time-of-impact
    toi = toi_interval.lower();

    return is_impacting;
}

////////////////////////////////////////////////////////////////////////////////
// Edge-Edge

double compute_edge_edge_tolerance_1d(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    const double input_tol)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        std::max(
            (edge0_vertex0_end - edge0_vertex0_start).norm(),
            (edge0_vertex1_end - edge0_vertex1_start).norm()),
        std::max(
            (edge1_vertex0_end - edge1_vertex0_start).norm(),
            (edge1_vertex1_end - edge1_vertex1_start).norm()));

    return input_tol / dl;
}

bool are_segments_intersecting(
    const Eigen::Vector3I& edge0_vertex0,
    const Eigen::Vector3I& edge0_vertex1,
    const Eigen::Vector3I& edge1_vertex0,
    const Eigen::Vector3I& edge1_vertex1)
{
    const Eigen::Vector3I& a = edge0_vertex0;
    const Eigen::Vector3I& b = edge0_vertex1;
    const Eigen::Vector3I& c = edge1_vertex0;
    const Eigen::Vector3I& d = edge1_vertex1;

    Eigen::Vector3I n = cross(b - a, c - a);

    return cross(b - d, c - d).dot(n).upper() >= 0
        && cross(c - d, a - d).dot(n).upper() >= 0
        && cross(a - d, b - d).dot(n).lower() <= 0;
}

bool edgeEdgeCCD_Redon(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    double& toi,
    const double &input_tol)
{
    const auto distance = [&](const Interval& t) {
        Eigen::Vector3I edge0_vertex0 = (edge0_vertex0_end.cast<Interval>()
                                         - edge0_vertex0_start.cast<Interval>())
                * t
            + edge0_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex1 = (edge0_vertex1_end.cast<Interval>()
                                         - edge0_vertex1_start.cast<Interval>())
                * t
            + edge0_vertex1_start.cast<Interval>();

        Eigen::Vector3I edge1_vertex0 = (edge1_vertex0_end.cast<Interval>()
                                         - edge1_vertex0_start.cast<Interval>())
                * t
            + edge1_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex1 = (edge1_vertex1_end.cast<Interval>()
                                         - edge1_vertex1_start.cast<Interval>())
                * t
            + edge1_vertex1_start.cast<Interval>();

        return cross(
                   edge0_vertex1 - edge0_vertex0, edge1_vertex1 - edge1_vertex0)
            .dot(edge1_vertex0 - edge0_vertex0);
    };

    const auto is_intersection_inside_edges = [&](const Interval& t) -> bool {
        Eigen::Vector3I edge0_vertex0 = (edge0_vertex0_end.cast<Interval>()
                                         - edge0_vertex0_start.cast<Interval>())
                * t
            + edge0_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex1 = (edge0_vertex1_end.cast<Interval>()
                                         - edge0_vertex1_start.cast<Interval>())
                * t
            + edge0_vertex1_start.cast<Interval>();

        Eigen::Vector3I edge1_vertex0 = (edge1_vertex0_end.cast<Interval>()
                                         - edge1_vertex0_start.cast<Interval>())
                * t
            + edge1_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex1 = (edge1_vertex1_end.cast<Interval>()
                                         - edge1_vertex1_start.cast<Interval>())
                * t
            + edge1_vertex1_start.cast<Interval>();

        return are_segments_intersecting(
            edge0_vertex0, edge0_vertex1, edge1_vertex0, edge1_vertex1);
    };

    double tol = compute_edge_edge_tolerance_1d(
        edge0_vertex0_start, edge0_vertex1_start, edge1_vertex0_start,
        edge1_vertex1_start, edge0_vertex0_end, edge0_vertex1_end,
        edge1_vertex0_end, edge1_vertex1_end, input_tol);

    Interval toi_interval;
    bool is_impacting = interval_root_finder(
        distance, is_intersection_inside_edges, Interval(0, 1), tol,
        toi_interval);

    // Return a conservative time-of-impact
    toi = toi_interval.lower();

    return is_impacting;
}

Eigen::Vector3d compute_face_vertex_tolerance_3d(
    const Eigen::Vector3d& vs,
    const Eigen::Vector3d& f0s,
    const Eigen::Vector3d& f1s,
    const Eigen::Vector3d& f2s,
    const Eigen::Vector3d& ve,
    const Eigen::Vector3d& f0e,
    const Eigen::Vector3d& f1e,
    const Eigen::Vector3d& f2e,
    const double input_tol)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        std::max((ve - vs).norm(), (f0e - f0s).norm()),
        std::max((f1e - f1s).norm(), (f2e - f2s).norm()));

    double edge0_length = std::max((f1s - f0s).norm(), (f1e - f0e).norm());
    double edge1_length = std::max((f2s - f0s).norm(), (f2e - f0e).norm());
    // double edge_length = std::max(edge0_length, edge1_length);

    return Eigen::Vector3d(
        input_tol / dl, input_tol / edge0_length,
        input_tol / edge1_length);
}

bool vertexFaceCCD_Interval(
    const Eigen::Vector3d& vertex_start,
    const Eigen::Vector3d& face_vertex0_start,
    const Eigen::Vector3d& face_vertex1_start,
    const Eigen::Vector3d& face_vertex2_start,
    const Eigen::Vector3d& vertex_end,
    const Eigen::Vector3d& face_vertex0_end,
    const Eigen::Vector3d& face_vertex1_end,
    const Eigen::Vector3d& face_vertex2_end,
    double& toi,
    const double &input_tol)
{
    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 3);
        Interval t = params(0);
        Interval u = params(1);
        Interval v = params(2);

        Eigen::Vector3I vertex = vertex_start.cast<Interval>()
            + (vertex_end.cast<Interval>() - vertex_start.cast<Interval>()) * t;

        Eigen::Vector3I face0 = face_vertex0_start.cast<Interval>()
            + (face_vertex0_end.cast<Interval>()
               - face_vertex0_start.cast<Interval>())
                * t;

        Eigen::Vector3I face1 = face_vertex1_start.cast<Interval>()
            + (face_vertex1_end.cast<Interval>()
               - face_vertex1_start.cast<Interval>())
                * t;

        Eigen::Vector3I face2 = face_vertex2_start.cast<Interval>()
            + (face_vertex2_end.cast<Interval>()
               - face_vertex2_start.cast<Interval>())
                * t;

        Eigen::Vector3I face_vertex = face0.cast<Interval>()
            + u * (face1.cast<Interval>() - face0.cast<Interval>())
            + v * (face2.cast<Interval>() - face0.cast<Interval>());

        return (vertex - face_vertex).eval();
    };

    Eigen::Vector3d tol = compute_face_vertex_tolerance_3d(
        vertex_start, face_vertex0_start, face_vertex1_start,
        face_vertex2_start, vertex_end, face_vertex0_end, face_vertex1_end,
        face_vertex2_end,input_tol);
    Eigen::VectorX3I toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,
        true);

    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;

    return 0;
}

    
    Eigen::Vector3d compute_edge_edge_tolerance(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    const double input_tol)
{
    // Compute the maximum trajectory length of all the vertices
    double dl = std::max(
        std::max(
            (edge0_vertex0_end - edge0_vertex0_start).norm(),
            (edge0_vertex1_end - edge0_vertex1_start).norm()),
        std::max(
            (edge1_vertex0_end - edge1_vertex0_start).norm(),
            (edge1_vertex1_end - edge1_vertex1_start).norm()));

    double edge0_length = std::max(
        (edge0_vertex1_start - edge0_vertex0_start).norm(),
        (edge0_vertex1_end - edge0_vertex0_end).norm());
    double edge1_length = std::max(
        (edge1_vertex1_start - edge1_vertex0_start).norm(),
        (edge1_vertex1_end - edge1_vertex0_end).norm());
    // double edge_length = std::max(edge0_length, edge1_length);

    return Eigen::Vector3d(
        input_tol / dl, input_tol / edge0_length,
        input_tol / edge1_length);
}


    bool edgeEdgeCCD_Interval(
    const Eigen::Vector3d& edge0_vertex0_start,
    const Eigen::Vector3d& edge0_vertex1_start,
    const Eigen::Vector3d& edge1_vertex0_start,
    const Eigen::Vector3d& edge1_vertex1_start,
    const Eigen::Vector3d& edge0_vertex0_end,
    const Eigen::Vector3d& edge0_vertex1_end,
    const Eigen::Vector3d& edge1_vertex0_end,
    const Eigen::Vector3d& edge1_vertex1_end,
    double& toi,
    const double &input_tol)
{
    const auto distance = [&](const Eigen::VectorX3I& params) {
        assert(params.size() == 3);
        Interval t = params(0);
        Interval edge0_alpha = params(1);
        Interval edge1_alpha = params(2);

        Eigen::Vector3I edge0_vertex0 = (edge0_vertex0_end.cast<Interval>()
                                         - edge0_vertex0_start.cast<Interval>())
                * t
            + edge0_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex1 = (edge0_vertex1_end.cast<Interval>()
                                         - edge0_vertex1_start.cast<Interval>())
                * t
            + edge0_vertex1_start.cast<Interval>();
        Eigen::Vector3I edge0_vertex
            = (edge0_vertex1 - edge0_vertex0) * edge0_alpha + edge0_vertex0;

        Eigen::Vector3I edge1_vertex0 = (edge1_vertex0_end.cast<Interval>()
                                         - edge1_vertex0_start.cast<Interval>())
                * t
            + edge1_vertex0_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex1 = (edge1_vertex1_end.cast<Interval>()
                                         - edge1_vertex1_start.cast<Interval>())
                * t
            + edge1_vertex1_start.cast<Interval>();
        Eigen::Vector3I edge1_vertex
            = (edge1_vertex1 - edge1_vertex0) * edge1_alpha + edge1_vertex0;

        return (edge1_vertex - edge0_vertex).eval();
    };

    Eigen::Vector3d tol = compute_edge_edge_tolerance(
        edge0_vertex0_start, edge0_vertex1_start, edge1_vertex0_start,
        edge1_vertex1_start, edge0_vertex0_end, edge0_vertex1_end,
        edge1_vertex0_end, edge1_vertex1_end,input_tol);

    Eigen::VectorX3I toi_interval;
    bool is_impacting = interval_root_finder(
        distance, Eigen::Vector3I::Constant(Interval(0, 1)), tol, toi_interval,
        false);

    // Return a conservative time-of-impact
    if (is_impacting) {
        toi = toi_interval(0).lower();
    }
    // This time of impact is very dangerous for convergence
    // assert(!is_impacting || toi > 0);
    return is_impacting;
}
    

    
} // namespace inclusion_ccd
