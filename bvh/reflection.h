//
// Created by Andrew Astakhov on 20.01.25.
//

#ifndef REFLECTION_H
#define REFLECTION_H
#include "vec.h"
#include "aabb.h"
#include "bvh.h"
#include "interval.h"
#include "raycast.h"
#include "prims.h"
#include <vector>
namespace bvh {
    namespace detail {
        inline void tri_normal(const Tri<vec3d> &tri, vec3d& normal
                                     ) {
            auto v1 = tri.b - tri.a;
            auto v2 = tri.c - tri.a;
            normal = v1.cross( v2);

            normal.unitize();


        }
        inline void reflect_one(const vec3d &direction , const vec3d &normal,vec3d &res) {
            res= direction - normal * 2.0 * direction.dot( normal) ;
        }
        // Parallel prefix sum (inclusive scan) of 'in' into 'out'
        inline void parallel_prefix_sum(const std::vector<size_t>& in, std::vector<size_t>& out)
        {
            out.resize(in.size());
            if (in.empty()) return;

            const size_t n = in.size();

            // Choose a block size (tune for your system/cache)
            const size_t blockSize = 1024;

            // Number of blocks
            size_t numBlocks = (n + blockSize - 1) / blockSize;

            // This will hold the sum of each block
            std::vector<size_t> blockSums(numBlocks);

            // 1) Compute partial sums within each block in parallel
#pragma omp parallel for
            for (size_t b = 0; b < numBlocks; ++b) {
                size_t start = b * blockSize;
                size_t end   = std::min(start + blockSize, n);

                size_t sum = 0;
                for (size_t i = start; i < end; ++i) {
                    sum += in[i];
                    out[i] = sum;
                }
                blockSums[b] = sum;
            }

            // 2) Compute the prefix sums of the block sums (sequentially or in parallel with a smaller scan)
            for (size_t b = 1; b < numBlocks; ++b) {
                blockSums[b] += blockSums[b - 1];
            }

            // 3) Add the block offset to each element in block b (for b > 0) in parallel
#pragma omp parallel for
            for (size_t b = 1; b < numBlocks; ++b) {
                size_t offset = blockSums[b - 1];
                size_t start  = b * blockSize;
                size_t end    = std::min(start + blockSize, n);

                for (size_t i = start; i < end; ++i) {
                    out[i] += offset;
                }
            }
        }
    }
inline void reflect(const std::vector<Ray<vec<double, 3> > > &rays,
    const BVH &bvh,
    const std::vector<Tri<vec3d> > &primitives,
    std::vector<Ray<vec3d > > &reflected,
    std::vector<bool> &mask) {

        std::vector<Hit> hits;
        std::vector<size_t> counts;
        auto normals=std::vector<vec3d>(primitives.size(), {0,0,0});
        mask.resize(rays.size(),true);


        reflected.resize(rays.size(), {{0.,0.,0.},{0.,0.,0.}});
        raycast_first(rays,bvh,primitives,hits,counts);
        std::vector<size_t> offsets(counts.size());


        for (size_t i = 0; i < rays.size(); ++i)
            {


                if (counts[i]==0) {
                    mask[i]=false;
                    continue;

                }

                

                const auto& hit=hits[i];
                reflected[i].start=hit.second;
                const auto& prim =primitives[hit.prim];
                auto& normal =normals[hit.prim];
                if (normal.zero()) {
                    detail::tri_normal(prim,normal);

                }
                detail::reflect_one(rays[i].direction,normal,reflected[i].direction);


            }
        }
}
#endif //REFLECTION_H
