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
    }


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

    inline void raycast(const Segm<vec<double, 3> > &segm,
                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives, std::vector<vec3d> &hits) {
        detail::raycast_internal(segm, bvh, 0, hits, primitives);
    }

    inline void raycast(const std::vector<Ray<vec<double, 3> > > &rays,

                        const BVH &bvh,
                        const std::vector<Tri<vec3d> > &primitives,
                        std::vector<vec3d> &hits,
                        std::vector<size_t> &offsets) {
        hits.reserve(rays.size() * 2);
        offsets.resize(rays.size());

        for (size_t i = 0; i < rays.size(); ++i) {
            const size_t size1 = hits.size();
            raycast(rays[i], bvh, primitives, hits);
            offsets[i] = hits.size() - size1;
        }
    }


}
#endif //RAYCAST_H
