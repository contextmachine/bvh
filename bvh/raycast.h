//
// Created by Andrew Astakhov on 15.01.25.
//

#ifndef RAYCAST_H
#define RAYCAST_H

#include "moller.h"
#include "aabb.h"
#include "bvh.h"
#include "prims.h"

namespace bvh {
    namespace detail {
        inline void raycast_internal(const Segm<vec<double, 3> > &inp, const BVH &bvh, size_t node_idx,
                                     std::vector<vec<double, 3> > &hits, const std::vector<Tri<vec3d> > &primitives) {
            auto &node = bvh.nodes[node_idx];
            Segm<vec<double, 3> > out;
            bool result = node.bbox.clip(inp, out);

            if (result) {
                if (node.isLeaf()) {
                    const Tri<vec<double, 3> > &primitive = primitives[node.object];
                    vec<double, 3> hit;
                    int res = intersect_triangle_segment<double>(primitive.a, primitive.b, primitive.c, inp.start,
                                                                 out.end, hit,
                                                                 1e-8);
                    if (res > 0) { hits.push_back(hit); }
                } else {

                    raycast_internal(inp, bvh, node.left, hits, primitives);
                    raycast_internal(inp, bvh, node.right, hits, primitives);
                }
            }
        }
        inline void raycast_internal(const Segm<vec<double, 3> > &inp, const BVH &bvh, size_t node_idx,
                                      const std::vector<Tri<vec3d> > &primitives, size_t &counts) {
            auto &node = bvh.nodes[node_idx];
            Segm<vec<double, 3> > out;
            bool result = node.bbox.clip(inp, out);

            if (result) {
                if (node.isLeaf()) {
                    const Tri<vec<double, 3> > &primitive = primitives[node.object];
                    vec<double, 3> hit;
                    int res = intersect_triangle_segment<double>(primitive.a, primitive.b, primitive.c, inp.start,
                                                                 out.end, hit,
                                                                 1e-8);
                    if (res > 0) { counts+=1; }
                } else {
                    raycast_internal(inp, bvh, node.left, primitives,counts);
                    raycast_internal(inp, bvh, node.right, primitives,counts);
                }
            }
        }


    }


    /**
     * Performs raycasting for a given ray against a BVH (Bounding Volume Hierarchy) and a set of triangular primitives.
     * This method computes the intersection points of the ray with the triangular primitives while utilizing the BVH
     * structure for efficient spatial query.
     *
     * @param ray The ray to be tested for intersections, specified by an origin and direction.
     * @param bvh The bounding volume hierarchy used for spatial partitioning and efficient intersection testing.
     * @param primitives A vector of triangular primitives against which the ray will be tested.
     * @param hits A vector to store the resulting intersection points. If the ray intersects any triangles, the points of intersection are appended to this vector.
     */
    inline void raycast(const Ray<vec<double, 3> > &ray,
                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives,
                        std::vector<vec3d> &hits) {
        double tmin = 0;
        double tmax = 0;
        bool res = bvh.nodes[0].bbox.intersectRay(ray, tmin, tmax);

        if (res == 0) { return; };

        Segm<vec<double, 3> > inp = {ray.start, ray.start + ray.direction * tmax};

        detail::raycast_internal(inp, bvh, 0, hits, primitives);

    }


    /**
     * Executes a raycasting operation to determine intersections of a given ray with a BVH (Bounding Volume Hierarchy) and a set of triangular primitives.
     * Uses the BVH to efficiently traverse the spatial structure and calculate potential intersections.
     *
     * @param ray The ray to be tested for intersections, defined by its starting point and direction.
     * @param bvh The bounding volume hierarchy used for organizing and accelerating intersection computations.
     * @param primitives A collection of triangular primitives against which the ray will be tested for intersection.
     * @param counts A reference to a variable tracking the total count of intersections detected during the raycasting process.
     */
    inline void raycast(const Ray<vec<double, 3> > &ray,
                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives,
                        size_t &counts) {
        double tmin = 0;
        double tmax = 0;
        bool res = bvh.nodes[0].bbox.intersectRay(ray, tmin, tmax);

        if (res == 0) { return; };

        Segm<vec<double, 3> > inp = {ray.start, ray.start + ray.direction * tmax};

        detail::raycast_internal(inp, bvh, 0,  primitives,counts);
    }
    /**
     * Performs raycasting for a segment against a BVH (Bounding Volume Hierarchy) and a set of primitives.
     * This method determines the intersection points where the segment intersects the primitives.
     *
     * @param segm The segment to be tested for intersections.
     * @param bvh The bounding volume hierarchy used for efficient intersection testing.
     * @param primitives A vector of triangular primitives against which the segment will be tested.
     * @param hits A vector to store the intersection points. The results of the segment intersections are appended to this vector.
     */
    inline void raycast(const Segm<vec<double, 3> > &segm,
                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives, std::vector<vec3d> &hits) {
        detail::raycast_internal(segm, bvh, 0, hits, primitives);
    }

    /**
     * Performs raycasting for a collection of rays against a BVH (Bounding Volume Hierarchy) and a set of primitives.
     * This method finds intersection points where rays intersect the primitives, calculates the hit points, and counts the number of hits per ray.
     *
     * @param rays A vector of rays to be tested for intersections.
     * @param bvh The bounding volume hierarchy used for efficient intersection testing.
     * @param primitives A vector of triangular primitives against which the rays will be tested.
     * @param hits A vector to store the intersection points. The results of the ray intersections are appended to this vector.
     * @param counts A vector storing the count of intersection hits for each ray. Its size will be equal to the size of the input ray vector.
     */
    inline void raycast(const std::vector<Ray<vec<double, 3> > > &rays,

                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives,
                        std::vector<vec3d> &hits,
                        std::vector<size_t> &counts) {
        hits.reserve(rays.size() * 2);
        counts.resize(rays.size());

        for (size_t i = 0; i < rays.size(); ++i) {
            const size_t size1 = hits.size();
            raycast(rays[i], bvh, primitives, hits);

            counts[i] = hits.size() - size1;
        }
    }


    /**
     * Performs batched raycasting for a collection of rays against a BVH (Bounding Volume Hierarchy) and a set of triangular primitives.
     * For each ray in the input collection, this method computes the number of intersections with the triangular primitives
     * while leveraging the BVH structure for efficient queries.
     *
     * @param rays A vector of rays, each specified by an origin and a direction, to be tested for intersections.
     * @param bvh The bounding volume hierarchy utilized for spatial partitioning and efficient intersection testing.
     * @param primitives A vector of triangular primitives against which each ray in the collection will be tested.
     * @param counts A vector to store the intersection counts for each ray. For each ray in the input, the corresponding entry in this vector
     *        is updated with the number of primitives that the ray intersects.
     */
    inline void raycast(const std::vector<Ray<vec<double, 3> > > &rays,

                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives,
                        std::vector<size_t> &counts) {
        counts.resize(rays.size());

        for (size_t i = 0; i < rays.size(); ++i) {
            raycast(rays[i], bvh, primitives, counts[i]);
        }
    }


}
#endif //RAYCAST_H
