//
// Created by Andrew Astakhov on 28.01.25.
//

#ifndef DISTANCE_H
#define DISTANCE_H
#include <cmath>
namespace bvh {
    namespace detail {
        #define sqr(x) ((x)*(x))
        template<typename Vec, typename T=typename Vec::value_type, size_t N=Vec::dim>
        constexpr inline T mag2(const Vec &a)
        {
            T l=std::sqrt(a[0]);
            for(unsigned int i=1; i<N; ++i) {l+=sqr(a[i]);};
            return l;
        }

    }
    template<typename T>
    constexpr inline T sign(const T val) {
        return (val > 0) ? 1 : ((val < 0) ? -1 : 0);
    }


    #define clamp(val, minVal, maxVal) std::max((minVal), std::min((val), (maxVal)))
    #define min3(a, b, c) std::min((a), std::min((b), (c)))
    #define DOT(a, b) ((a).dot((b)))
    /************************************
     *  Signed Distance (from original code)
     ************************************/


    template<typename Vec, typename T=typename Vec::value_type>
    constexpr inline T point_segment_distance(const Vec &x0, const Vec &x1, const Vec &x2)
    {
        Vec dx(x2-x1);
        double m2=detail::mag2(dx);
        // find parameter value of closest point on segment
        T s12=DOT(x2-x0 , dx)/m2;
        if(s12<0){
            s12=0;
        }else if(s12>1){
            s12=1;
        }
        // and find the distance
        return (x0-(x1*s12+x2*(1-s12))).length();
    }
    template<typename Vec, typename T=typename Vec::value_type>
    constexpr inline T point_triangle_distance(const Vec& p, const Vec& a, const Vec& b, const Vec& c) {
        Vec ba = b - a;
        Vec pa = p - a;
        Vec cb = c - b;
        Vec pb = p - b;
        Vec ac = a - c;
        Vec pc = p - c;
        Vec nor = ba.cross(ac);

        T signSum = sign(ba.cross(nor).dot( pa)) +
                        sign(cb.cross(nor).dot(  pb) )+
                        sign(ac.cross(nor).dot( pc));

        if (signSum < 2) {
            // Compute distances to the edges
            // Clamp the projection of p onto each edge to the segment [0,1]
            T t1 = clamp(DOT(ba, pa) /ba.sqLength(), 0, 1);
            Vec edge1 = ba * t1 - pa;
            T dist1 =edge1.dot(edge1);

            T t2 = clamp(DOT(cb, pb) /cb.sqLength(), 0, 1);
            Vec edge2 = cb * t2 - pb;
            T dist2 = edge2.dot( edge2);

            T t3 = clamp(DOT(ac, pc) / ac.sqLength(), 0, 1);
            Vec edge3 = ac * t3 - pc;
            T dist3 = edge3.dot( edge3);

            // Return the square root of the minimum distance squared
            T minDistSq = min3(dist1, dist2, dist3);
            return std::sqrt(minDistSq);
        } else {
            // Compute the distance to the plane of the triangle
            T dist = DOT(nor, pa);
            return std::sqrt((dist * dist) / nor.sqLength());
        }
    }
}
#endif //DISTANCE_H
