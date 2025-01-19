//
// Created by Andrew Astakhov on 18.01.25.
//

#ifndef MOLLER_H
#define MOLLER_H
#include <cmath>
#include <src/defines.h>
#include <src/vec.h>

namespace bvh {
    namespace detail {
        template<typename T>
        inline bool is_close(T a, T b, T eps = CMMCORE_EPSILON<T>()) noexcept {
            return std::abs(a - b) < eps;
        };

        template<typename T>
        inline bool point_equals(vec<T, 3> &p, vec<T, 3> &q, T eps = CMMCORE_EPSILON<T>()) noexcept {
            return (is_close(p[0], q[0], eps) and is_close(p[1], q[1], eps) and is_close(p[2], q[2], eps));
        }

        template<typename T>
        inline bool point_on_edge(vec<T, 3> &P, vec<T, 3> &V0, vec<T, 3> &V1,
                                  const T eps = CMMCORE_EPSILON<T>()) noexcept {
            T cross_x, cross_y, cross_z, len_sq, dot_prod;
            T vec_px = P[0] - V0[0];
            T vec_py = P[1] - V0[1];
            T vec_pz = P[2] - V0[2];
            T edge_x = V1[0] - V0[0];
            T edge_y = V1[1] - V0[1];
            T edge_z = V1[2] - V0[2];

            // Compute cross product to check collinearity
            cross_x = vec_py * edge_z - vec_pz * edge_y;
            cross_y = vec_pz * edge_x - vec_px * edge_z;
            cross_z = vec_px * edge_y - vec_py * edge_x;

            // If cross product is not near zero, P is not on the line
            bool cond = is_close(cross_x, 0.0, eps) and is_close(cross_y, 0.0, eps) and is_close(cross_z, 0.0, eps);
            if (!cond) {
                return 0;
            }

            // Compute dot product to check if P is between V0 and V1
            dot_prod = vec_px * edge_x + vec_py * edge_y + vec_pz * edge_z;
            len_sq = edge_x * edge_x + edge_y * edge_y + edge_z * edge_z;

            if (dot_prod < -eps or dot_prod > len_sq + eps) {
                return 0;
            }

            return 1;
        };
    };


    template<typename T>
    inline int classify_intersection(vec<T, 3> &P, vec<T, 3> &V0, vec<T, 3> &V1, vec<T, 3> &V2,
                                     const T eps = CMMCORE_EPSILON<T>()) noexcept {
        /*
        Classify the intersection point P with the triangle defined by V0, V1, V2.
        Returns:
            1: Intersection within the interior
            2: Intersection at a vertex
            3: Intersection along an edge
        */
        if (point_equals(P, V0, eps)) {
            return 2;
        }
        if (point_equals(P, V1, eps)) {
            return 2;
        }
        if (point_equals(P, V2, eps)) {
            return 2;
        }
        if (point_on_edge(P, V0, V1, eps)) {
            return 3;
        }
        if (point_on_edge(P, V1, V2, eps)) {
            return 3;
        }
        if (point_on_edge(P, V2, V0, eps)) {
            return 3;
        }

        return 1;
    }


    template<typename T>
    inline int intersect_triangle_segment(
        const vec<T, 3> &V0, const vec<T, 3> &V1, const vec<T, 3> &V2,
        const vec<T, 3> &S, const vec<T, 3> &E,
        vec<T, 3> &I, const T eps = CMMCORE_EPSILON<T>()) noexcept {
        /*
        Möller–Trumbore intersection algorithm adapted for segments.
        Returns:
            1 if intersection exists and is within [0,1] segment parameters, else 0.
            If intersection exists, fills I with intersection point.
        */
        T dir_x, dir_y, dir_z;
        T edge1_x, edge1_y, edge1_z;
        T edge2_x, edge2_y, edge2_z;
        T h_x, h_y, h_z;
        T a, f, u, v;
        T s_x, s_y, s_z;
        T q_x, q_y, q_z;
        T t;

        // Compute direction vector of segment
        dir_x = E[0] - S[0];
        dir_y = E[1] - S[1];
        dir_z = E[2] - S[2];

        // Find vectors for two edges sharing V0
        edge1_x = V1[0] - V0[0];
        edge1_y = V1[1] - V0[1];
        edge1_z = V1[2] - V0[2];

        edge2_x = V2[0] - V0[0];
        edge2_y = V2[1] - V0[1];
        edge2_z = V2[2] - V0[2];

        // Begin calculating determinant - also used to calculate u parameter
        h_x = dir_y * edge2_z - dir_z * edge2_y;
        h_y = dir_z * edge2_x - dir_x * edge2_z;
        h_z = dir_x * edge2_y - dir_y * edge2_x;

        a = edge1_x * h_x + edge1_y * h_y + edge1_z * h_z;

        if ((-eps < a) and (a < eps)) {
            return 0; // This means parallel
        }

        f = 1.0 / a;
        s_x = S[0] - V0[0];
        s_y = S[1] - V0[1];
        s_z = S[2] - V0[2];

        u = f * (s_x * h_x + s_y * h_y + s_z * h_z);
        if (u < -eps or u > 1.0 + eps) {
            return 0;
        }

        q_x = s_y * edge1_z - s_z * edge1_y;
        q_y = s_z * edge1_x - s_x * edge1_z;
        q_z = s_x * edge1_y - s_y * edge1_x;

        v = f * (dir_x * q_x + dir_y * q_y + dir_z * q_z);
        if (v < -eps or (u + v) > (1.0 + eps)) {
            return 0;
        }

        // At this stage we can compute t to find out where the intersection point is on the line
        t = f * (edge2_x * q_x + edge2_y * q_y + edge2_z * q_z);

        if (t < -eps or t > 1.0 + eps) {
            return 0; // Intersection not within the segment
        };
        // Compute intersection point
        I[0] = S[0] + t * dir_x;
        I[1] = S[1] + t * dir_y;
        I[2] = S[2] + t * dir_z;

        return 1;
    };
}
#endif //MOLLER_H
