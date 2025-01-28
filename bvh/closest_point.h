//
// Created by Andrew Astakhov on 29.01.25.
//

#ifndef CLOSEST_POINT_H
#define CLOSEST_POINT_H
#include <cmath>
#include "vec.h"
#include "distance.h"
namespace bvh {
template<typename Vec, typename T=typename Vec::value_type>
// Function to compute the closest point on a line segment AB to the origin
constexpr inline Vec closestPointOnLine(const Vec& a, const Vec& b) {
    Vec ab = b - a;
    T t = -(a.dot(ab)) / ab.sqLength();
    t = clamp(t, 0.0, 1.0);
    return a + ab * t;
}
    template<typename Vec, typename T=typename Vec::value_type>
// Function to compute the closest point on a triangle ABC to the origin
constexpr inline Vec closestPointOnTriangle(const Vec& a, const Vec& b, const Vec& c) {
    // Barycentric coordinates method
    Vec ab = b - a;
    Vec ac = c - a;
    Vec ao = -a;

    T d1 = ab.dot(ao);
    T d2 = ac.dot(ao);
    if (d1 <= 0 && d2 <= 0) return a; // Vertex region A

    Vec bo = -b;
    T d3 = ab.dot(bo);
    T d4 = ac.dot(bo);
    if (d3 >= 0 && d4 <= d3) return b; // Vertex region B

    T vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        T v = d1 / (d1 - d3);
        return a + ab * v; // Edge region AB
    }

    Vec co = -c;
    T d5 = ab.dot(co);
    T d6 = ac.dot(co);
    if (d6 >= 0 && d5 <= d6) return c; // Vertex region C

    T vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        T w = d2 / (d2 - d6);
        return a + ac * w; // Edge region AC
    }

    T va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w; // Edge region BC
    }

    // Inside face region
    T denom = 1.0 / (va + vb + vc);
    T v = vb * denom;
    T w = vc * denom;
    return a + ab * v + ac * w;
}
    template<typename Vec, typename T=typename Vec::value_type>
constexpr inline Vec closestPointOnTriangle(const Vec p, const Vec a, const Vec b, const Vec c)
{
    // Check if P in vertex region outside A
    Vec ab = b - a;
    Vec ac = c - a;
    Vec ap = p - a;
    T d1 = ab.dot( ap);
    T d2 = ac.dot( ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a; // barycentric coordinates (1,0,0)
    // Check if P in vertex region outside B
    Vec bp = p - b;
    T d3 = ab.dot( bp);
    T d4 = ac.dot( bp);
    if (d3 >= 0.0 && d4 <= d3) return b; // barycentric coordinates (0,1,0)
    // Check if P in edge region of AB, if so return projection of P onto AB
    T vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        T v = d1 / (d1 - d3);
        return a + v * ab; // barycentric coordinates (1-v,v,0)
    }
    // Check if P in vertex region outside C
    Vec cp = p - c;
    T d5 = ab.dot( cp);
    T d6 = ac.dot( cp);
    if (d6 >= 0.0 && d5 <= d6) return c; // barycentric coordinates (0,0,1)
    // Check if P in edge region of AC, if so return projection of P onto AC
    T vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        T w = d2 / (d2 - d6);
        return a + w * ac; // barycentric coordinates (1-w,0,w)
    }
    // Check if P in edge region of BC, if so return projection of P onto BC
    T va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b); // barycentric coordinates (0,1-w,w)
    }
    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    T denom = 1.0 / (va + vb + vc);
    T v = vb * denom;
    T w = vc * denom;
    return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0 - v - w
}
}
#endif //CLOSEST_POINT_H
