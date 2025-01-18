
#ifndef AABB_HEADER_H
#define AABB_HEADER_H

#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <ostream>
#include <type_traits>
#include <algorithm> // for std::min/max


namespace bvh {
namespace detail {

template<typename Vec, typename T=typename Vec::value_type>
T sdBBox(const Vec &p,const Vec &b) {
        Vec d = p.abs() - b;
        return std::min(d.max_val(),(T)0)+d.max((T)0).length();
    }

}
/*****************************************************************************************
 *  Primary Template AABB< vec<T,N> >
 *****************************************************************************************/
template <typename T,size_t N>
class AABB
{
public:
    using vec_type    = vec<T,N>;                     // e.g. vec<T,N>
    using scalar_type = T;

    vec_type min;
    vec_type max;

    vec_type centroid;

public:
    /************************************
     * Constructors
     ************************************/
     AABB()
    {
        // Set min to +max, max to +lowest => "reversed" bounding box
        // so any real point expansions will fix them properly.
        // We'll do dimension-based initialization:
        //for(std::size_t i=0; i<sizeof...(scalar_type); /* not correct, we must do N from Vec. */) {}

       
        for(std::size_t i=0; i<N; ++i) {
            min[i] = std::numeric_limits<scalar_type>::max();
            max[i] = std::numeric_limits<scalar_type>::lowest();
        }
        centroid.set(0);
    }

    // Construct from explicit min, max
     AABB(const vec<scalar_type,N>& vmin, const vec<scalar_type,N>& vmax)
        : min(vmin), max(vmax)
    {
        updateCentroid();
    }

    // Construct from list of points
     AABB(const std::vector<vec<scalar_type,N>>& pts)
    {
        if(pts.empty()) {
            // define an "empty" box
            for(std::size_t i=0; i<N; ++i) {
                min[i] = std::numeric_limits<scalar_type>::max();
                max[i] = std::numeric_limits<scalar_type>::lowest();
            }
            updateCentroid();
            return;
        }
        min = pts[0];
        max = pts[0];
        for(std::size_t i=1; i<pts.size(); ++i) {
            expand(pts[i]);
        }
    }

    /************************************
     *  Accessors
     ************************************/




    /************************************
     *  Intersection Check
     ************************************/
    // Overlaps if each dimension's intervals overlap
     bool intersects(const AABB& other) const
    {
        for(std::size_t i=0; i<N; ++i) {
            if(min[i] > other.max[i] || max[i] < other.min[i]) {
                return false;
            }
        }
        return true;
    }

    // Intersection region, if any
     bool intersection(const AABB& other, AABB& result) const
    {
        if(!intersects(other)) {
            return false;
        }
        vec_type newMin, newMax;
        for(std::size_t i=0; i<N; ++i) {
            newMin[i] = (min[i] > other.min[i]) ? min[i] : other.min[i];
            newMax[i] = (max[i] < other.max[i]) ? max[i] : other.max[i];
        }
        result = AABB(newMin, newMax);
        return true;
    }

    /************************************
     *  Volume
     ************************************/
    // volume = product of (max[i] - min[i]) for i in [0..N)
    // If any dimension is negative (min>max), the volume is 0.
     scalar_type volume() const
    {
        scalar_type vol = scalar_type(1);
        for(std::size_t i=0; i<N; ++i) {
            scalar_type d = max[i] - min[i];
            if(d < scalar_type(0)) {
                return scalar_type(0); // invalid or collapsed box
            }
            vol *= d;
        }
        return vol;
    }

    /************************************
     *  Merge
     ************************************/
    // Merged bounding box of this and 'other'
     AABB merge(const AABB& other) const
    {

        vec_type mmin = min;
        vec_type mmax = max;
        for(std::size_t i=0; i<N; ++i) {
            mmin[i] = (mmin[i] < other.min[i]) ? mmin[i] : other.min[i];
            mmax[i] = (mmax[i] > other.max[i]) ? mmax[i] : other.max[i];
        }
        return AABB(mmin, mmax);
    }

    /************************************
     *  Expand
     ************************************/
    // Expand by a point
     void expand(const vec_type& pt)
    {
        for(std::size_t i=0; i<N; ++i) {
            if(pt[i] < min[i]) min[i] = pt[i];
            if(pt[i] > max[i]) max[i] = pt[i];
        }
        updateCentroid();
    }

    // Expand by another AABB
     void expand(const AABB& bb)
    {
        for(std::size_t i=0; i<N; ++i) {
            if(bb.min[i] < min[i]) min[i] = bb.min[i];
            if(bb.max[i] > max[i]) max[i] = bb.max[i];
        }
        updateCentroid();
    }

    T sd(const vec_type& pt) const
     {
         vec_type cnt=(min+max)/
         2;return detail::sdBBox<vec_type>(pt-cnt,max-cnt);
     }

    // Recompute centroid = (min+max)/2
     void updateCentroid()
    {
        centroid = (min + max) * scalar_type(0.5);
    }
};

/*****************************************************************************************
 *  Partial Specialization for AABB< vec<T,2> > (2D)
 *  (If we want special extra performance or special checks)
 *****************************************************************************************/
template <typename T>
class AABB<T,2>
{
public:
    using vec_type    = vec<T,2>;
    using scalar_type = T;


    vec_type min, max, centroid;

public:
     AABB()
    {
        min.x = std::numeric_limits<T>::max();
        min.y = std::numeric_limits<T>::max();
        max.x = std::numeric_limits<T>::lowest();
        max.y = std::numeric_limits<T>::lowest();
        updateCentroid();
    }
     AABB(const vec_type& vmin, const vec_type& vmax)
        : min(vmin), max(vmax)
    {
        updateCentroid();
    }
     AABB(const std::vector<vec_type>& pts)
    {
        if(pts.empty()) {
            min.x = std::numeric_limits<T>::max();
            min.y = std::numeric_limits<T>::max();
            max.x = std::numeric_limits<T>::lowest();
            max.y = std::numeric_limits<T>::lowest();
            updateCentroid();
            return;
        }
        min = pts[0];
        max = pts[0];
        for(std::size_t i=1; i<pts.size(); ++i) {
            expand(pts[i]);
        }
    }




     bool intersects(const AABB& other) const
    {
        // overlap if [minX <= other.maxX, maxX >= other.minX] etc
        if(min.x > other.max.x) return false;
        if(max.x < other.min.x) return false;
        if(min.y > other.max.y) return false;
        if(max.y < other.min.y) return false;
        return true;
    }
     bool intersection(const AABB& other, AABB& result) const
    {
        if(!intersects(other)) {
            return false;
        }
        result.min.x = (min.x > other.min.x)? min.x: other.min.x;
        result.min.y = (min.y > other.min.y)? min.y: other.min.y;
        result.max.x = (max.x < other.max.x)? max.x: other.max.x;
        result.max.y = (max.y < other.max.y)? max.y: other.max.y;
        result.updateCentroid();
        return true;
    }

     scalar_type volume() const
    {
        // 2D "volume" => area
        T dx = max.x - min.x;
        T dy = max.y - min.y;
        if(dx < T(0) || dy < T(0)) {
            return T(0);
        }
        return dx * dy;
    }
     AABB merge(const AABB& other) const
    {
        AABB ret(*this);
        ret.expand(other);
        return ret;
    }
     void expand(const vec_type& pt)
    {
        if(pt.x < min.x) min.x = pt.x;
        if(pt.y < min.y) min.y = pt.y;
        if(pt.x > max.x) max.x = pt.x;
        if(pt.y > max.y) max.y = pt.y;
        updateCentroid();
    }
     void expand(const AABB& bb)
    {
        if(bb.min.x < min.x) min.x = bb.min.x;
        if(bb.min.y < min.y) min.y = bb.min.y;
        if(bb.max.x > max.x) max.x = bb.max.x;
        if(bb.max.y > max.y) max.y = bb.max.y;
        updateCentroid();
    }



     void updateCentroid()
    {
        centroid.x = (min.x + max.x)*T(0.5);
        centroid.y = (min.y + max.y)*T(0.5);
    }

};

/*****************************************************************************************
 *  Partial Specialization for AABB< vec<T,3> > (3D)
 *****************************************************************************************/
template <typename T>
class AABB< T,3 >
{
public:
    using vec_type    = vec<T,3>;
    using scalar_type = T;



    vec_type min, max, centroid;

public:
     AABB()
    {
        min.x = std::numeric_limits<T>::max();
        min.y = std::numeric_limits<T>::max();
        min.z = std::numeric_limits<T>::max();
        max.x = std::numeric_limits<T>::lowest();
        max.y = std::numeric_limits<T>::lowest();
        max.z = std::numeric_limits<T>::lowest();
        updateCentroid();
    }
     AABB(const vec_type& vmin, const vec_type& vmax)
        : min(vmin.x,vmin.y,vmin.z), max(vmax.x,vmax.y,vmax.z)
    {
        updateCentroid();
    }
     AABB(const std::vector<vec_type>& pts)
    {
        if(pts.empty()) {
            min.x = std::numeric_limits<T>::max();
            min.y = std::numeric_limits<T>::max();
            min.z = std::numeric_limits<T>::max();
            max.x = std::numeric_limits<T>::lowest();
            max.y = std::numeric_limits<T>::lowest();
            max.z = std::numeric_limits<T>::lowest();
            updateCentroid();
            return;
        }
        min = pts[0];
        max = pts[0];
        for(std::size_t i=1; i<pts.size(); ++i) {
            expand(pts[i]);
        }
    }

    AABB(const AABB< T,3> &bb) {
         min=bb.min;
         max=bb.max;
         updateCentroid();

     }



    // Overlap test
     bool intersects(const AABB& other) const
    {
        if(min.x > other.max.x) return false;
        if(max.x < other.min.x) return false;
        if(min.y > other.max.y) return false;
        if(max.y < other.min.y) return false;
        if(min.z > other.max.z) return false;
        if(max.z < other.min.z) return false;
        return true;
    }
    // Intersection AABB
     bool intersection(const AABB& other, AABB& result) const
    {
        if(!intersects(other)) {
            return false;
        }
        vec_type rmin, rmax;
        rmin.x = (min.x > other.min.x)? min.x: other.min.x;
        rmin.y = (min.y > other.min.y)? min.y: other.min.y;
        rmin.z = (min.z > other.min.z)? min.z: other.min.z;

        rmax.x = (max.x < other.max.x)? max.x: other.max.x;
        rmax.y = (max.y < other.max.y)? max.y: other.max.y;
        rmax.z = (max.z < other.max.z)? max.z: other.max.z;

        result = AABB(rmin, rmax);
        return true;
    }
    // Volume = (dx * dy * dz)
     scalar_type volume() const
    {
        T dx = max.x - min.x;
        T dy = max.y - min.y;
        T dz = max.z - min.z;
        if(dx < T(0) || dy < T(0) || dz < T(0)) {
            return T(0);
        }
        return dx*dy*dz;
    }

    // Merge
     AABB merge(const AABB<T,3>& other) const
    {
        auto ret=AABB(min, max);
        ret.expand(other);
        return ret;
    }

    // Expand by a point
     void expand(const vec_type& pt)
    {
        if(pt.x < min.x) min.x = pt.x;
        if(pt.y < min.y) min.y = pt.y;
        if(pt.z < min.z) min.z = pt.z;
        if(pt.x > max.x) max.x = pt.x;
        if(pt.y > max.y) max.y = pt.y;
        if(pt.z > max.z) max.z = pt.z;
        updateCentroid();
    }
    // Expand by another AABB
     void expand(const AABB<T,3>& bb)
    {
        if(bb.min.x <= min.x) min.x = bb.min.x;
        if(bb.min.y <= min.y) min.y = bb.min.y;
        if(bb.min.z <= min.z) min.z = bb.min.z;
        if(bb.max.x > max.x) max.x = bb.max.x;
        if(bb.max.y > max.y) max.y = bb.max.y;
        if(bb.max.z > max.z) max.z = bb.max.z;
        updateCentroid();
    }




     void updateCentroid()
    {
        centroid.x = (min.x + max.x)*T(0.5);
        centroid.y = (min.y + max.y)*T(0.5);
        centroid.z = (min.z + max.z)*T(0.5);
    }
};

/*****************************************************************************************
 *  Partial Specialization for AABB< vec<T,4> > (4D)
 *****************************************************************************************/
template <typename T>
class AABB< T,4 >
{
public:
    using vec_type    = vec<T,4>;
    using scalar_type = T;


    vec_type min, max, centroid;

public:
     AABB()
    {
        // We'll treat w as just another dimension here,
        // but typically 4D might be used for something else (homogeneous coords).
        // We'll just do the standard approach:
        min.x = std::numeric_limits<T>::max();
        min.y = std::numeric_limits<T>::max();
        min.z = std::numeric_limits<T>::max();
        min.w = std::numeric_limits<T>::max();
        max.x = std::numeric_limits<T>::lowest();
        max.y = std::numeric_limits<T>::lowest();
        max.z = std::numeric_limits<T>::lowest();
        max.w = std::numeric_limits<T>::lowest();
        updateCentroid();
    }
     AABB(const vec_type& vmin, const vec_type& vmax)
        : min(vmin), max(vmax)
    {
        updateCentroid();
    }
     AABB(const std::vector<vec_type>& pts)
    {
        if(pts.empty()) {
            min.x = std::numeric_limits<T>::max();
            min.y = std::numeric_limits<T>::max();
            min.z = std::numeric_limits<T>::max();
            min.w = std::numeric_limits<T>::max();
            max.x = std::numeric_limits<T>::lowest();
            max.y = std::numeric_limits<T>::lowest();
            max.z = std::numeric_limits<T>::lowest();
            max.w = std::numeric_limits<T>::lowest();
            updateCentroid();
            return;
        }
        min = pts[0];
        max = pts[0];
        for(std::size_t i=1; i<pts.size(); ++i) {
            expand(pts[i]);
        }
    }




    // Overlap test in 4D
     bool intersects(const AABB& other) const
    {
        if(min.x > other.max.x) return false;
        if(max.x < other.min.x) return false;
        if(min.y > other.max.y) return false;
        if(max.y < other.min.y) return false;
        if(min.z > other.max.z) return false;
        if(max.z < other.min.z) return false;
        if(min.w > other.max.w) return false;
        if(max.w < other.min.w) return false;
        return true;
    }
     bool intersection(const AABB& other, AABB& result) const
    {
        if(!intersects(other)) {
            return false;
        }
        vec_type rmin, rmax;
        rmin.x = (min.x > other.min.x)? min.x: other.min.x;
        rmin.y = (min.y > other.min.y)? min.y: other.min.y;
        rmin.z = (min.z > other.min.z)? min.z: other.min.z;
        rmin.w = (min.w > other.min.w)? min.w: other.min.w;

        rmax.x = (max.x < other.max.x)? max.x: other.max.x;
        rmax.y = (max.y < other.max.y)? max.y: other.max.y;
        rmax.z = (max.z < other.max.z)? max.z: other.max.z;
        rmax.w = (max.w < other.max.w)? max.w: other.max.w;

        result = AABB(rmin, rmax);
        return true;
    }

    // "Volume" in 4D => hyper-volume
     scalar_type volume() const
    {
        T dx = max.x - min.x;
        T dy = max.y - min.y;
        T dz = max.z - min.z;
        T dw = max.w - min.w;
        if(dx < T(0) || dy < T(0) || dz < T(0) || dw < T(0)) {
            return T(0);
        }
        return dx*dy*dz*dw;
    }

     AABB merge(const AABB& other) const
    {
        AABB ret(*this);
        ret.expand(other);
        return ret;
    }
     void expand(const vec_type& pt)
    {
        if(pt.x < min.x) min.x = pt.x;
        if(pt.y < min.y) min.y = pt.y;
        if(pt.z < min.z) min.z = pt.z;
        if(pt.w < min.w) min.w = pt.w;
        if(pt.x > max.x) max.x = pt.x;
        if(pt.y > max.y) max.y = pt.y;
        if(pt.z > max.z) max.z = pt.z;
        if(pt.w > max.w) max.w = pt.w;
        updateCentroid();
    }
     void expand(const AABB& bb)
    {
        if(bb.min.x < min.x) min.x = bb.min.x;
        if(bb.min.y < min.y) min.y = bb.min.y;
        if(bb.min.z < min.z) min.z = bb.min.z;
        if(bb.min.w < min.w) min.w = bb.min.w;
        if(bb.max.x > max.x) max.x = bb.max.x;
        if(bb.max.y > max.y) max.y = bb.max.y;
        if(bb.max.z > max.z) max.z = bb.max.z;
        if(bb.max.w > max.w) max.w = bb.max.w;
        updateCentroid();
    }




     void updateCentroid()
    {
        centroid.x = (min.x + max.x)*T(0.5);
        centroid.y = (min.y + max.y)*T(0.5);
        centroid.z = (min.z + max.z)*T(0.5);
        centroid.w = (min.w + max.w)*T(0.5);
    }
};

/*****************************************************************************************
 *  Operator<< for AABB
 *****************************************************************************************/
template <typename T,size_t N>
inline std::ostream& operator<<(std::ostream& os, const AABB<T,N>& obj)
{
    os << "[" << obj.min << ", " << obj.max << "]";
    return os;
}

using AABB3d = AABB<double,3>;

} // end namespace bvh

#endif // AABB_HEADER_H
