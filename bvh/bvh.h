#ifndef BVH_BUILDER_H
#define BVH_BUILDER_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <limits>
#include "vec.h"    // Must include the provided "vec.h"
#include "aabb.h"   // Must include the provided "aabb.h"
#include "prims.h"

namespace bvh {
    /*****************************************************************************************
     *  Using AABB for 3D (common case)
     *****************************************************************************************/
    using AABB3d = AABB<vec3d>;



    /*****************************************************************************************
     *  BVH Node
     *****************************************************************************************/
    struct BVHNode {
        AABB3d bbox;
        int left = -1;
        int right = -1;
        int object = -1; // index of the object (leaf), -1 if it's an internal node

        bool isLeaf() const {
            return object != -1;
        }
    };

    /*****************************************************************************************
     *  Helpers for node indexing in a binary-heap layout
     *****************************************************************************************/
    constexpr inline size_t leftChild(const size_t i) noexcept { return 2 * i + 1; }
    constexpr inline size_t rightChild(const size_t i) noexcept { return 2 * i + 2; }


    template<typename Iter>
    Iter split_objects(Iter begin, Iter end) {
        if (begin == end) {
            throw std::runtime_error("split_objects: empty range");
        }

        // Compute the bounding of centroids
        auto minC = begin->second.centroid;
        auto maxC = begin->second.centroid;
        for (auto it = std::next(begin); it != end; ++it) {
            const auto &c = it->second.centroid;
            if (c.x < minC.x) minC.x = c.x;
            if (c.y < minC.y) minC.y = c.y;
            if (c.z < minC.z) minC.z = c.z;

            if (c.x > maxC.x) maxC.x = c.x;
            if (c.y > maxC.y) maxC.y = c.y;
            if (c.z > maxC.z) maxC.z = c.z;
        }

        bvh::vec3 diff = maxC - minC;
        int axis = max_axis(diff); // 0 => x, 1 => y, 2 => z

        // Sort range [begin, end) by centroid[axis]
        std::sort(begin, end, [axis](auto &a, auto &b) {
            return a.second.centroid[axis] < b.second.centroid[axis];
        });

        // Return the midpoint
        auto total = std::distance(begin, end);
        auto midOffset = total / 2;
        return std::next(begin, midOffset);
    }


    /*
     *  BVH class represents a Bounding Volume Hierarchy (BVH) structure used to efficiently
     *  organize and query spatial data, such as bounding boxes (AABB3d).
    */
    class BVH {
    public:
        std::vector<BVHNode> nodes;
        size_t root_index = 0; // index of root node in 'nodes'

        BVH() = default;

        BVHNode &getRoot() {
            return nodes[root_index];
        }

        AABB3d bbox() const {
            return getRoot().bbox;
        }

        const BVHNode &getRoot() const {
            return nodes[root_index];
        }
        bool empty() const {
            return nodes.empty();
            }



    private:
        // _build_bvh_internal with iterators
        template<typename Iter>
        size_t _build_bvh_internal(Iter begin, Iter end, const size_t currentIndex) {
            auto count = std::distance(begin, end);
            if (count <= 0) {
                throw std::runtime_error("_build_bvh_internal: zero objects in range");
            }

            // Ensure 'nodes' is large enough
            if (currentIndex >= nodes.size()) {
                nodes.resize(currentIndex + 1);
            }
            //BVHNode& node = nodes[currentIndex];

            if (count == 1) {
                auto &node = nodes[currentIndex];
                // Leaf node
                node.bbox = begin->second;
                node.object = begin->first;
                node.left = -1;
                node.right = -1;
                return currentIndex;
            }
            // Otherwise, we split
            auto mid = split_objects(begin, end);

            // Build children
            nodes[currentIndex].left = _build_bvh_internal(begin, mid, leftChild(currentIndex));
            nodes[currentIndex].right = _build_bvh_internal(mid, end, rightChild(currentIndex));
            auto &node = nodes[currentIndex];
            // Merge bounding boxes from children
            node.bbox = nodes[nodes[currentIndex].left].bbox.merge(nodes[nodes[currentIndex].right].bbox);
            node.object = -1; // internal node

            return currentIndex;
        }

        inline std::ostream &print_internal(std::ostream &os, const size_t node_index = 0) const {
            auto &node = nodes[node_index];
            os << node.bbox << ",";
            if (node.isLeaf()) {
                return os;
            }

            print_internal(os, node.left);
            print_internal(os, node.right);
            return os;
        }

    public:
        /*
         * Builds a Bounding Volume Hierarchy (BVH) from a list of axis-aligned bounding boxes (AABB).
         *
         * @param bboxes A vector of 3D axis-aligned bounding boxes (AABB3d) representing geometric objects.
         *               Each AABB corresponds to the bounding box of an individual object.
         * @return A reference to the BVH object after constructing and organizing the tree structure.
         *
         * If the input vector of bounding boxes is empty, the BVH is cleared,
         * and the root index is reset to zero.
         * The method reserves space for all potential BVH nodes and uses a recursive function to construct
         * the hierarchy by partitioning the AABBs into spatial subdivisions.
         */
        BVH &build(const std::vector<AABB3d> &bboxes) {
            if (bboxes.empty()) {
                nodes.clear();
                root_index = 0;
                return *this;
            }
            // Reserve space for up to (2*N - 1) nodes
            std::size_t n = bboxes.size();
            nodes.resize(2 * n - 1);

            // Create temporary array of (objectIndex, AABB3d)
            std::vector<std::pair<int, AABB3d> > objects;
            objects.reserve(n);
            for (std::size_t i = 0; i < n; ++i) {
                objects.push_back({static_cast<int>(i), bboxes[i]});
            }

            // Build recursively
            root_index = _build_bvh_internal(objects.begin(), objects.end(), 0);
            return *this;
        }


        template<typename TP>
        BVH &build(const std::vector<TP > &primitives) {

            std::vector<AABB3d> bboxes(primitives.size());

            for (size_t i = 0; i < primitives.size(); ++i) {
                auto& prim = primitives[i];
                prim.bbox(bboxes[i]);

            }
            return build(bboxes);
        }

        BVH &build(const std::vector<Tri<vec3d> > &primitives) {
            std::vector<AABB3d> bboxes(primitives.size());
            for (size_t i = 0; i < primitives.size(); ++i) {
                bboxes[i].expand(primitives[i].a);
                bboxes[i].expand(primitives[i].b);
                bboxes[i].expand(primitives[i].c);
            }
            return build(bboxes);
        }

        inline friend std::ostream &operator<<(std::ostream &os, const BVH &obj) {
            os << "[";
            obj.print_internal(os, 0);
            os << "]";
            return os;
        }
        
    };
} // end namespace bvh
#endif
