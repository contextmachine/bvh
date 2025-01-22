//
// Created by Andrew Astakhov on 15.01.25.
//

#ifndef RAYCAST_H
#define RAYCAST_H

#include "moller.h"
#include "aabb.h"
#include "bvh.h"
#include "prims.h"
#include <omp.h> // Include this for OpenMP
namespace bvh {
    struct Hit {

        double first;
        vec<double,3> second;
        size_t prim;

    };
    namespace detail {

        inline void raycast_internal(const Segm<vec<double, 3> > &inp, const vec<double,3> &dir, const BVH &bvh, size_t node_idx,
                                     std::vector<Hit> &hits, const std::vector<Tri<vec3d> > &primitives) {
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

                    if (res > 0) {

                        hits.push_back( {dir.dot(hit-inp.start),hit,(size_t)node.object});
                    }
                } else {

                    raycast_internal(inp,dir, bvh, node.left, hits, primitives);
                    raycast_internal(inp, dir,bvh, node.right, hits, primitives);
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
        inline bool raycast_internal(const Segm<vec<double, 3> > &inp, const BVH &bvh, size_t node_idx,
                                         const std::vector<Tri<vec3d> > &primitives) {
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
                    if (res > 0) { return true;}
                } else {
                    bool res=raycast_internal(inp, bvh, node.left, primitives);
                    if (res) {
                        return true;
                    }

                    res=raycast_internal(inp, bvh, node.right, primitives);
                    if (res) {
                        return true;
                    }
                }


            }
            return false;
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
                        std::vector<Hit> &hits) {
        double tmin = 0;
        double tmax = 0;
        bool res = bvh.nodes[0].bbox.intersectRay(ray, tmin, tmax);

        if (res == 0) { return; };

        Segm<vec<double, 3> > inp = {ray.start, ray.start + ray.direction * tmax};

        detail::raycast_internal(inp,ray.direction, bvh, 0, hits, primitives);
        std::sort(hits.begin(), hits.end(),[](auto const& lhs, auto const& rhs) {
            return lhs.first < rhs.first;
        });







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
                        const std::vector<Tri<vec3d> > &primitives, std::vector<Hit> &hits) {
        detail::raycast_internal(segm,segm.end-segm.start, bvh, 0, hits, primitives);
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
                        std::vector<Hit> &hits,
                        std::vector<size_t> &counts) {
        hits.reserve(rays.size() * 2);
        counts.resize(rays.size());

        for (size_t i = 0; i < rays.size(); ++i) {
            std::vector<Hit> tmp;
            double tmin = 0;
            double tmax = 0;
            auto& ray=rays[i];
            bool res = bvh.nodes[0].bbox.intersectRay(ray, tmin, tmax);

            if (res == 0) { continue; };

            Segm<vec<double, 3> > inp = {ray.start, ray.start + ray.direction * tmax};

            detail::raycast_internal(inp,ray.direction, bvh, 0, tmp, primitives);
            counts[i] = tmp.size();

            std::sort(tmp.begin(), tmp.end(),[](auto const& lhs, auto const& rhs) {
                return lhs.first < rhs.first;
            });

            for (int j = 0; j <   counts[i]; ++j) {
                hits.push_back(tmp[j]);
            }

        }
    }
    inline void raycast_first(const std::vector<Ray<vec<double, 3> > > &rays,

                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives,
                        std::vector<Hit> &hits,
                        std::vector<size_t> &counts) {
        hits.resize(rays.size());
        counts.resize(rays.size());

#pragma omp parallel for
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(rays.size()); ++i) {
            std::vector<Hit> tmp;
            double tmin = 0;
            double tmax = 0;
            auto& ray=rays[i];
            bool res = bvh.nodes[0].bbox.intersectRay(ray, tmin, tmax);

            if (res == 0) { continue; };

            Segm<vec<double, 3> > inp = {ray.start, ray.start + ray.direction * tmax};

            detail::raycast_internal(inp,ray.direction, bvh, 0, tmp, primitives);
            counts[i] = tmp.size();
            if (counts[i]==0) {
                continue;
            }
            std::sort(tmp.begin(), tmp.end(),[](auto const& lhs, auto const& rhs) {
                return lhs.first < rhs.first;
            });
            hits[i]=tmp[0];
            


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
                        std::vector<size_t> &counts
                        ) {

        counts.clear();
        counts.resize(rays.size());
        const size_t rootIndex = bvh.root_index;
        // Early exit if BVH has no nodes
        if (bvh.nodes.empty()) {
            return;
        }

        #pragma omp parallel for
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(rays.size()); ++i) {
            double tmin, tmax;
            bool intersectsBBox = bvh.nodes[rootIndex].bbox.intersectRay(rays[i], tmin, tmax);

            if (!intersectsBBox) {
                counts[i]=0;
                continue;
            }

            // Clip the ray to the bounding box for consistent intersection checks
            Segm<vec<double, 3>> segment;
            segment.start = rays[i].start;
            segment.end   = rays[i].start + rays[i].direction * tmax;

            // Collect intersections in allHits[i]
            detail::raycast_internal(segment, bvh, rootIndex, primitives,counts[i]);

        }
    }
#include <vector>
#include <omp.h> // Include this for OpenMP
#include "vec.h"
#include "moller.h"
#include "aabb.h"
#include "bvh.h"
#include "prims.h"


    /**
     * @brief Parallel raycast multiple rays against the BVH and collect intersections.
     *
     * @param rays         A vector of rays to be tested against the BVH.
     * @param bvh          The BVH data structure (already built).
     * @param primitives   The triangles (or other primitives) referenced by the BVH.
     * @param allHits      A vector of hit-lists, one entry per ray in `rays`.
     */
    inline void raycast_single_omp(const std::vector<Ray<vec<double, 3>>> &rays,
                            const BVH &bvh,
                            const std::vector<Tri<vec3d>> &primitives,
                            std::vector<bool> &mask)
    {
        // Resize and clear output
        mask.clear();
        mask.resize(rays.size());

        // Early exit if BVH has no nodes
        if (bvh.nodes.empty()) {
            return;
        }

        const size_t rootIndex = bvh.root_index;

        // Parallel loop over all rays
        #pragma omp parallel for
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(rays.size()); ++i) {
            double tmin, tmax;
            bool intersectsBBox = bvh.nodes[rootIndex].bbox.intersectRay(rays[i], tmin, tmax);

            if (!intersectsBBox) {
                mask[i]=false;
                continue;
            }

            // Clip the ray to the bounding box for consistent intersection checks
            Segm<vec<double, 3>> segment;
            segment.start = rays[i].start;
            segment.end   = rays[i].start + rays[i].direction * tmax;

            // Collect intersections in allHits[i]
            mask[i]=detail::raycast_internal(segment, bvh, rootIndex, primitives);

        }
    }

} // namespace bvh


#endif //RAYCAST_H
