//
// Created by Andrew Astakhov on 27.01.25.
//

#ifndef RAYCAST_EIGEN_H
#define RAYCAST_EIGEN_H
#include "vec.h"
#include "prims.h"
#include "bvh.h"
#include "moller.h"
#include "Eigen/Core"
#include "Eigen/Sparse"
namespace bvh {
    namespace detail {

        // Then use it in your map:

        inline void raycast_bvh_internal(const size_t ray_id,
            const Ray<vec<double, 3> > &inp, const BVH &bvh, const std::vector<Tri<vec<double,3>>> &primitives,const size_t node_idx,
                                    Eigen::SparseMatrix<Eigen::Vector4d> &prims_rays_hits, const double eps=1e-8) {


            //bool result = node.bbox.clip(inp, out);
            double start,end;
            bool result =bvh.nodes[node_idx].bbox.intersectRay(inp,start,end);


            if (result) {
                auto &node = bvh.nodes[node_idx];
                if (node.isLeaf()) {

                    if (start>eps) {
                        auto& prim=primitives[node.object];
                        auto end_pt=inp.start+(inp.direction*end);
                        vec3d resv;
                        int res=intersect_triangle_segment(prim.a,prim.b,prim.c,inp.start,end_pt, resv,eps);
                        if (res>0) {

                            prims_rays_hits.insert(node.object,ray_id)={
                                resv.x,
                                resv.y,
                                resv.z,
                                inp.direction.dot(resv-inp.start)
                            };

                        }
                    }

                } else {

                    raycast_bvh_internal(ray_id,inp, bvh, primitives, bvh.nodes[node_idx].left, prims_rays_hits, eps);
                    raycast_bvh_internal(ray_id,inp,bvh,  primitives,bvh.nodes[node_idx].right, prims_rays_hits, eps);
                }
            }
        }

}




}
#endif //RAYCAST_EIGEN_H
