#ifndef AABB_HEADER_H
#define AABB_HEADER_H

#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <ostream>
#include <type_traits>
#include <algorithm>
#include "prims.h"
namespace bvh {
namespace detail {

/**
 * Helper for "signed distance" to box, used by sd(...) method:
 * This is not necessary for slab intersection,
 * but was in your original code, so we keep it.
 */
template<typename Vec, typename T=typename Vec::value_type>
T sdBBox(const Vec &p,const Vec &b) {
        Vec d = p.abs() - b;
    return std::min(d.max_val(), T(0)) + d.max(T(0)).length();
    }

} // end namespace detail

/*****************************************************************************************
 *  Primary Template AABB< vec<T,N> >
 *****************************************************************************************/
template <typename T,size_t N>
/**
 * Axis-Aligned Bounding Box (AABB) class.
 *
 * Represents a bounding box defined by its minimum and maximum corners in N-dimensional space.
 * Used for various geometric computations, intersection tests, and bounding operations.
 *
 * Template parameters:
 * - T: Scalar type used for coordinates (e.g., float, double).
 * - N: Number of dimensions of the bounding box.
 */
class AABB {
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
        // Initialize "empty" box so any real expansions fix them
        for(std::size_t i=0; i<N; ++i) {
            min[i] = std::numeric_limits<scalar_type>::max();
            max[i] = std::numeric_limits<scalar_type>::lowest();
        }
        centroid.set(0);
    }

    // Construct from explicit min, max
    AABB(const vec_type& vmin, const vec_type& vmax)
        : min(vmin), max(vmax)
    {
        updateCentroid();
    }

    // Construct from list of points
    AABB(const std::vector<vec_type>& pts)
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
     *  Intersection Check (AABB vs AABB)
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
            if(other.min[i] < mmin[i]) mmin[i] = other.min[i];
            if(other.max[i] > mmax[i]) mmax[i] = other.max[i];
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

    /************************************
     *  Signed Distance (not strictly needed for intersection, but from original code)
     ************************************/
    T sd(const vec_type& pt) const
    {
        vec_type cnt = (min + max) / 2;
        return detail::sdBBox<vec_type>(pt - cnt, max - cnt);
    }

    /************************************
     *  Update centroid
     ************************************/
    void updateCentroid()
    {
        centroid = (min + max) * scalar_type(0.5);
    }

    /************************************
     *  Ray vs. AABB Intersection (Generic N-D Slab Test)
     *
     *  Returns true if intersection occurs; tEnter/tExit define param range.
     ************************************/
    bool intersectRay(const Ray<vec_type>& ray, T& tEnter, T& tExit) const
    {
        tEnter = -std::numeric_limits<T>::infinity();
        tExit  =  std::numeric_limits<T>::infinity();

        for (size_t i = 0; i < N; i++)
        {
            T dir = ray.direction[i];
            T start = ray.start[i];

            // If direction is almost zero, check if we lie within the slab
            if (std::fabs(dir) < std::numeric_limits<T>::epsilon())
            {
                if (start < min[i] || start > max[i]) {
                    return false;
                }
            }
            else
            {
                T invDir = T(1) / dir;
                T t0 = (min[i] - start) * invDir;
                T t1 = (max[i] - start) * invDir;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }
        return true;
    }
    /**
  * Clips the forward ray (t>=0) to the portion inside the AABB.
  * Returns true if the clipped segment is nonempty, false otherwise.
  *
  * On success, segm will contain the portion of the ray inside the box.
  */
    bool clip(const Ray<vec_type> &ray, Segm<vec_type> &segm) const
    {
        T tEnter = -std::numeric_limits<T>::infinity();
        T tExit  =  std::numeric_limits<T>::infinity();

        // For each dimension, clip the ray against [min[i], max[i]]
        for (size_t i = 0; i < N; ++i)
        {
            T startCoord = ray.start[i];
            T dir        = ray.direction[i];
            T minCoord   = min[i];
            T maxCoord   = max[i];

            // If direction is ~0, we must be within the slab for that axis
            if (std::fabs(dir) < std::numeric_limits<T>::epsilon())
            {
                // Check if outside slab
                if (startCoord < minCoord || startCoord > maxCoord)
                    return false; // No intersection
            }
            else
            {
                // Compute intersection parameters with the bounding planes
                T invD = T(1) / dir;
                T t0 = (minCoord - startCoord) * invD;
                T t1 = (maxCoord - startCoord) * invD;
                if (t0 > t1) std::swap(t0, t1);

                // Shrink our [tEnter, tExit] interval
                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;

                // If the interval is empty, no intersection
                if (tEnter > tExit)
                    return false;
            }
        }

        // If the entire intersection is behind the start of the ray, no intersection
        if (tExit < 0)
            return false;

        // If the portion starts behind, clamp it to 0 so we only keep t >= 0
        if (tEnter < 0)
            tEnter = 0;

        // Now build the segment
        segm.start = ray.start + ray.direction * tEnter;
        segm.end   = ray.start + ray.direction * tExit;
        return true;
    }


};



/*****************************************************************************************
 *  Partial Specialization for AABB< vec<T,2> >
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

    // Slab-based Ray vs. AABB intersection (2D)
    bool intersectRay(const Ray<vec_type>& ray, T& tEnter, T& tExit) const
    {
        tEnter = -std::numeric_limits<T>::infinity();
        tExit  =  std::numeric_limits<T>::infinity();

        // X dimension
        if (std::fabs(ray.direction.x) < std::numeric_limits<T>::epsilon())
        {
            if (ray.start.x < min.x || ray.start.x > max.x)
                return false;
        }
        else
        {
            T invDx = T(1) / ray.direction.x;
            T t0 = (min.x - ray.start.x) * invDx;
            T t1 = (max.x - ray.start.x) * invDx;
            if (t0 > t1) std::swap(t0, t1);
            if (t0 > tEnter) tEnter = t0;
            if (t1 < tExit)  tExit  = t1;
            if (tEnter > tExit) return false;
        }

        // Y dimension
        if (std::fabs(ray.direction.y) < std::numeric_limits<T>::epsilon())
        {
            if (ray.start.y < min.y || ray.start.y > max.y)
                return false;
        }
        else
        {
            T invDy = T(1) / ray.direction.y;
            T t0 = (min.y - ray.start.y) * invDy;
            T t1 = (max.y - ray.start.y) * invDy;
            if (t0 > t1) std::swap(t0, t1);
            if (t0 > tEnter) tEnter = t0;
            if (t1 < tExit)  tExit  = t1;
            if (tEnter > tExit) return false;
        }

        return true;
    }

};

/*****************************************************************************************
 *  Partial Specialization for AABB< vec<T,3> >
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

    AABB(const AABB& bb)
    {
         min=bb.min;
         max=bb.max;
         updateCentroid();

     }

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

    AABB merge(const AABB& other) const
    {
        auto ret=AABB(min, max);
        ret.expand(other);
        return ret;
    }

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

    void expand(const AABB& bb)
    {
        if(bb.min.x < min.x) min.x = bb.min.x;
        if(bb.min.y < min.y) min.y = bb.min.y;
        if(bb.min.z < min.z) min.z = bb.min.z;
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

    // Slab-based Ray vs. AABB intersection (3D)
    bool intersectRay(const Ray<vec_type>& ray, T& tEnter, T& tExit) const
    {
        tEnter = -std::numeric_limits<T>::infinity();
        tExit  =  std::numeric_limits<T>::infinity();

        // X
        if (std::fabs(ray.direction.x) < std::numeric_limits<T>::epsilon())
        {
            if (ray.start.x < min.x || ray.start.x > max.x)
                return false;
        }
        else
        {
            T invDx = T(1) / ray.direction.x;
            T t0 = (min.x - ray.start.x) * invDx;
            T t1 = (max.x - ray.start.x) * invDx;
            if (t0 > t1) std::swap(t0, t1);
            if (t0 > tEnter) tEnter = t0;
            if (t1 < tExit)  tExit  = t1;
            if (tEnter > tExit) return false;
        }

        // Y
        if (std::fabs(ray.direction.y) < std::numeric_limits<T>::epsilon())
        {
            if (ray.start.y < min.y || ray.start.y > max.y)
                return false;
        }
        else
        {
            T invDy = T(1) / ray.direction.y;
            T t0 = (min.y - ray.start.y) * invDy;
            T t1 = (max.y - ray.start.y) * invDy;
            if (t0 > t1) std::swap(t0, t1);
            if (t0 > tEnter) tEnter = t0;
            if (t1 < tExit)  tExit  = t1;
            if (tEnter > tExit) return false;
        }

        // Z
        if (std::fabs(ray.direction.z) < std::numeric_limits<T>::epsilon())
        {
            if (ray.start.z < min.z || ray.start.z > max.z)
                return false;
        }
        else
        {
            T invDz = T(1) / ray.direction.z;
            T t0 = (min.z - ray.start.z) * invDz;
            T t1 = (max.z - ray.start.z) * invDz;
            if (t0 > t1) std::swap(t0, t1);
            if (t0 > tEnter) tEnter = t0;
            if (t1 < tExit)  tExit  = t1;
            if (tEnter > tExit) return false;
        }

        return true;
    }

    /**
     * 3D-specific clip method that clamps a forward ray (t >= 0)
     * to the portion inside this AABB. Returns true if the
     * clipped segment is non-empty. segm will hold the resulting segment.
     */
    bool clip(const Ray<vec_type> &ray, Segm<vec_type> &segm) const
    {
        T tEnter = -std::numeric_limits<T>::infinity();
        T tExit  =  std::numeric_limits<T>::infinity();

        //
        // --- X dimension ---
        //
        {
            const T startX = ray.start.x;
            const T dirX   = ray.direction.x;

            if (std::fabs(dirX) < std::numeric_limits<T>::epsilon())
            {
                // Ray is parallel in X; must be within [min.x, max.x]
                if (startX < min.x || startX > max.x)
                    return false;
            }
            else
            {
                const T invDx = T(1) / dirX;
                T t0 = (min.x - startX) * invDx;
                T t1 = (max.x - startX) * invDx;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }

        //
        // --- Y dimension ---
        //
        {
            const T startY = ray.start.y;
            const T dirY   = ray.direction.y;

            if (std::fabs(dirY) < std::numeric_limits<T>::epsilon())
            {
                // Ray is parallel in Y; must be within [min.y, max.y]
                if (startY < min.y || startY > max.y)
                    return false;
            }
            else
            {
                const T invDy = T(1) / dirY;
                T t0 = (min.y - startY) * invDy;
                T t1 = (max.y - startY) * invDy;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }

        //
        // --- Z dimension ---
        //
        {
            const T startZ = ray.start.z;
            const T dirZ   = ray.direction.z;

            if (std::fabs(dirZ) < std::numeric_limits<T>::epsilon())
            {
                // Ray is parallel in Z; must be within [min.z, max.z]
                if (startZ < min.z || startZ > max.z)
                    return false;
            }
            else
            {
                const T invDz = T(1) / dirZ;
                T t0 = (min.z - startZ) * invDz;
                T t1 = (max.z - startZ) * invDz;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }

        // If the entire intersection is behind the origin, no intersection
        if (tExit < T(0))
            return false;

        // If intersection starts behind origin, clamp it to 0 for forward ray
        if (tEnter < T(0))
            tEnter = T(0);

        // Build the resulting segment
        segm.start = ray.start + ray.direction * tEnter;
        segm.end   = ray.start + ray.direction * tExit;

        return true;
    }
/**
     * Clips the input segment 'inp' (s -> e) to the portion inside this AABB (3D).
     * Returns true if the clipped segment is non-empty.
     * 'segm' will hold the result on success.
     */
    bool clip(const Segm<vec_type>& inp, Segm<vec_type>& segm) const
    {
        // Param range [tEnter, tExit] within [0,1]
        T tEnter = T(0);
        T tExit  = T(1);

        // s + dir * t
        vec_type s   = inp.start;
        vec_type dir = inp.end - inp.start;

        // ---- X slab ----
        {
            const T startX = s.x;
            const T dirX   = dir.x;
            const T minX   = min.x;
            const T maxX   = max.x;

            if (std::fabs(dirX) < std::numeric_limits<T>::epsilon())
            {
                if (startX < minX || startX > maxX)
                    return false;
            }
            else
            {
                T invDx = T(1) / dirX;
                T t0 = (minX - startX) * invDx;
                T t1 = (maxX - startX) * invDx;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }

        // ---- Y slab ----
        {
            const T startY = s.y;
            const T dirY   = dir.y;
            const T minY   = min.y;
            const T maxY   = max.y;

            if (std::fabs(dirY) < std::numeric_limits<T>::epsilon())
            {
                if (startY < minY || startY > maxY)
                    return false;
            }
            else
            {
                T invDy = T(1) / dirY;
                T t0 = (minY - startY) * invDy;
                T t1 = (maxY - startY) * invDy;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }

        // ---- Z slab ----
        {
            const T startZ = s.z;
            const T dirZ   = dir.z;
            const T minZ   = min.z;
            const T maxZ   = max.z;

            if (std::fabs(dirZ) < std::numeric_limits<T>::epsilon())
            {
                if (startZ < minZ || startZ > maxZ)
                    return false;
            }
            else
            {
                T invDz = T(1) / dirZ;
                T t0 = (minZ - startZ) * invDz;
                T t1 = (maxZ - startZ) * invDz;
                if (t0 > t1) std::swap(t0, t1);

                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
        }

        // Must be overlapping [0,1]
        if (tExit < T(0) || tEnter > T(1))
            return false;

        // Clamp
        if (tEnter < T(0)) tEnter = T(0);
        if (tExit  > T(1)) tExit  = T(1);

        if (tEnter > tExit)
            return false;

        // Build clipped segment
        segm.start = s + dir * tEnter;
        segm.end   = s + dir * tExit;
        return true;
    }

};

/*****************************************************************************************
 *  Partial Specialization for AABB< vec<T,4> >
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

    // Slab-based Ray vs. AABB intersection (4D)
    bool intersectRay(const Ray<vec_type>& ray, T& tEnter, T& tExit) const
    {
        tEnter = -std::numeric_limits<T>::infinity();
        tExit  =  std::numeric_limits<T>::infinity();

        // helper to do each dimension
        auto checkAxis = [&](T rayStart, T rayDir, T minVal, T maxVal) -> bool
        {
            if (std::fabs(rayDir) < std::numeric_limits<T>::epsilon())
            {
                // Ray parallel in this dimension: must be within slab
                if (rayStart < minVal || rayStart > maxVal)
                    return false;
            }
            else
            {
                T invD = T(1) / rayDir;
                T t0 = (minVal - rayStart) * invD;
                T t1 = (maxVal - rayStart) * invD;
                if (t0 > t1) std::swap(t0, t1);
                if (t0 > tEnter) tEnter = t0;
                if (t1 < tExit)  tExit  = t1;
                if (tEnter > tExit) return false;
            }
            return true;
        };

        if(!checkAxis(ray.start.x, ray.direction.x, min.x, max.x)) return false;
        if(!checkAxis(ray.start.y, ray.direction.y, min.y, max.y)) return false;
        if(!checkAxis(ray.start.z, ray.direction.z, min.z, max.z)) return false;
        if(!checkAxis(ray.start.w, ray.direction.w, min.w, max.w)) return false;

        return true;
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

// Example typedef
using AABB3d = AABB<double,3>;

} // end namespace bvh

#endif // AABB_HEADER_H
