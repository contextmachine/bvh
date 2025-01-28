/*
Copyright (c) 2024 Andrew Astakhov <aa@contextmachine.ru>. All rights reserved.

This software and associated documentation files (the "Software") are
proprietary and confidential. Unauthorized copying, modification, distribution,
or usage is strictly prohibited unless expressly permitted by the author or
organization.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef MESH_H
#define MESH_H



#include <memory>
#include "vec.h"
#include "aabb.h"
#include "closest_point.h"
#include "bvh.h"


#include "raycast.h"
#include <array>
#include <assert.h>
#include <map>
#include <iostream>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <cstdint>
#include <stdexcept>


#include <vector>

#include "defines.h"

namespace bvh {
#define BVH_MESH_COLLISION_EPS 1e-8


  inline void calculateEdgesFromFaces(const std::vector<std::array<int,3>>& faces, std::vector<std::array<int,2>>& edges) {
    edges.clear();
    edges.reserve(faces.size() * 3); // Reserve to avoid multiple allocations.

    // Collect all edges
    for (const auto& face : faces) {
      std::array<int,2> e1 = {std::min(face[0], face[1]), std::max(face[0], face[1])};
      std::array<int,2> e2 = {std::min(face[1], face[2]), std::max(face[1], face[2])};
      std::array<int,2> e3 = {std::min(face[2], face[0]), std::max(face[2], face[0])};

      edges.push_back(e1);
      edges.push_back(e2);
      edges.push_back(e3);
    }

    // Sort edges lexicographically
    std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
        return (a[0] < b[0]) || (a[0] == b[0] && a[1] < b[1]);
    });

    // Remove duplicates
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  }

  using Face = std::array<int, 3>;

  using Vertex = vec;
  template <typename Vec> class Mesh;

  namespace detail {
  }
  template <typename Vec>
/**
 * @class Mesh
 * @brief Represents a 3D mesh with vertices, faces, and edges.
 *
 * The Mesh class is designed to handle 3D geometric data structures. It contains
 * a set of vertices, faces, and edges that constitute a mesh. The class supports
 * operations such as mesh construction, copying, moving, serialization, deserialization,
 * and calculation of geometric features like edges and Euler characteristics.
 */
class Mesh {
  public:
    using vec_type=Vec;
    using value_type = typename Vec::value_type;
    constexpr static size_t dim=Vec::dim;
    std::vector<vec_type> vertices_;
    std::vector<Face> faces_;
  private:
    struct FaceProxy {
      const Mesh* owner ;
      const std::array<int,3>& face;

      FaceProxy(const Mesh* own,const size_t num): owner(own), face(owner->faces_[num]) {

      }
      ~FaceProxy()=default;

      void bbox(AABB<vec_type>& bb) {
        auto& v1=owner->vertices_[face[0]];
        auto& v2=owner->vertices_[face[1]];
        auto& v3=owner->vertices_[face[2]];
        bb.expand(v1);
        bb.expand(v2);
        bb.expand(v3);
      }


    };
    inline void raycast_internal(const Ray<vec_type > &inp,  size_t const node_idx,size_t &counts, const value_type eps=CMMCORE_EPSILON<value_type>())const {
      auto &node = bvh.nodes[node_idx];
      value_type tmin,tmax;
      bool result = node.bbox.intersectRay(inp, tmin,tmax);

      if (result) {
        if (node.isLeaf()) {
          auto& face=faces_[node.object];
          auto& v0=vertices_[face[0]];
          auto& v1=vertices_[face[1]];
          auto& v2=vertices_[face[2]];

          vec_type hit;
          int res = intersect_triangle_segment<double>(v0, v1, v2,
            inp.start,
              inp.start+inp.direction*(tmax+eps),
                                                       hit,
                                                       eps);

          if (res > 0) {
            counts+=1;


          }


        } else {

          raycast_internal(inp,  node.left,counts, eps);

          raycast_internal(inp, node.right,counts, eps);





        }
      }

    }

  public:



    std::vector<std::array<int,2>> edges_{};
    int euler = -1;
    BVH<Vec> bvh;
    int id = 0;

    Mesh() =default;


    Mesh(const Mesh &other) : vertices_(other.vertices_), faces_(other.faces_), edges_(other.edges_), euler(other.euler),
                              id(other.id) {
      if (other.bvh) {
        bvh = std::make_unique<BVHNode>(*other.bvh);
      }
    }

    Mesh(std::vector<Vec> &vxs, std::vector<Face> &fcs): bvh(), vertices_(vxs), faces_(fcs)  {
    }

    // Constructor accepting positions and indices (indexed geometry)
    Mesh(const std::vector<value_type> &positions, const std::vector<int> &indices)
      : bvh() {
      // Build the vertex list from positions
      size_t numVertices = positions.size() / 3;
      vertices_.reserve(numVertices);
      for (size_t i = 0; i < positions.size(); i += 3) {
        vertices_.emplace_back(positions[i], positions[i + 1], positions[i + 2]);
      }

      // Build the face list from indices
      size_t numFaces = indices.size() / 3;
      faces_.reserve(numFaces);
      for (size_t i = 0; i < indices.size(); i += 3) {
        faces_.emplace_back(Face{indices[i], indices[i + 1], indices[i + 2]});
      }
    }





    Mesh &operator=(const Mesh &other) {
      if (this != &other) {
        vertices_ = other.vertices_;
        faces_ = other.faces_;
        edges_ = other.edges_;
        euler = other.euler;
        id = other.id;
        if (other.bvh) {
          bvh = std::make_unique<BVHNode>(*other.bvh);
        } else {
          bvh.reset();
        }
      }
      return *this;
    }


    Mesh(Mesh &&other) noexcept : vertices_(std::move(other.vertices_)),
                                  faces_(std::move(other.faces_)),
                                  edges_(std::move(other.edges_)),
                                  euler(other.euler),
                                  bvh(std::move(other.bvh)),
                                  id(other.id) {
    }

    Mesh &operator=(Mesh &&other) noexcept {
      if (this != &other) {
        vertices_ = std::move(other.vertices_);
        faces_ = std::move(other.faces_);
        edges_ = std::move(other.edges_);
        euler = other.euler;
        id = other.id;
        bvh = std::move(other.bvh);
      }
      return *this;
    }

    bool checkGeometryEquality(const Mesh &other) const {
      if (vertices_.size() != other.vertices_.size()) {
        return false;
      }
      for (size_t i = 0; i < vertices_.size(); ++i) {
        if (!(vertices_[i] == other.vertices_[i])) {
          return false;
        }
      }
      return true;
    }

    void rebuildFeatures() {
      rebuildEdges();
      rebuildEulerCharacteristic();
      rebuildBVH();
    }

    void rebuildEdges() {

      calculateEdgesFromFaces(faces_, edges_);
    }

    void rebuildEulerCharacteristic() {
      euler = vertices_.size() + faces_.size() - edges_.size();
    }

    /**
     * @brief Determines if the mesh is a closed manifold.
     *
     * This method checks if the mesh is a closed manifold by verifying
     * its Euler characteristic and ensuring the edges are correctly constructed.
     * If the edges are insufficient or uninitialized, it invokes a rebuild
     * operation on both the edges and the Euler characteristic.
     *
     * @return true if the mesh is a closed manifold, false otherwise.
     */
    bool isClosedManifold() {

      if (euler == -1) {
        rebuildEulerCharacteristic();
      }
      return (euler == 2);
    }

    // Get a reference to the set of unique vertices (without copying)
    const std::vector<vec_type> &getVertices() const { return vertices_; }

    // Get a reference to the set of faces (without copying)
    const std::vector<Face> &getFaces() const { return faces_; }



    AABB<vec_type>& getAABB()  const {
      if (hasBVH()) {
        rebuildBVH();
      }
      return bvh.bbox();
    }

    size_t serialize(char *buffer) const {
      char *ptr = buffer;


      // Placeholder for size (will be filled later)
      ptr += sizeof(uint64_t);

      // Serialize number of vertices
      uint64_t numVertices = vertices_.size();
      std::memcpy(ptr, &numVertices, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      // Serialize vertex data
      for (const auto &vertex: vertices_) {
        std::memcpy(ptr, &vertex.x, sizeof(value_type) * 3);
        ptr += sizeof(vec_type::value_type) * 3;
      }

      // Serialize number of faces
      uint64_t numFaces = faces_.size();
      std::memcpy(ptr, &numFaces, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      // Serialize face data
      for (const auto &face: faces_) {
        std::memcpy(ptr, &face[0], sizeof(int) * 3);
        ptr += sizeof(int) * 3;
      }

      // Calculate total bytes written (excluding size placeholder)
      size_t totalSize = ptr - buffer - sizeof(uint64_t);

      // Now go back and write the totalSize at the beginning
      std::memcpy(buffer, &totalSize, sizeof(uint64_t));

      // Return total bytes written including the size field
      return ptr - buffer;
    }

    // Deserialization method (updated to read size prefix)
    size_t deserialize(const char *buffer) {
      const char *ptr = buffer;

      // Read the size of the Mesh data
      uint64_t meshSize = 0;
      std::memcpy(&meshSize, ptr, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      const char *meshDataEnd = ptr + meshSize;

      // Deserialize number of vertices
      uint64_t numVertices;
      std::memcpy(&numVertices, ptr, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      // Resize and deserialize vertex data
      this->vertices_.resize(numVertices);
      for (size_t i = 0; i < numVertices; ++i) {
        std::memcpy(&this->vertices_[i].x, ptr, sizeof(value_type) * 3);
        ptr += sizeof(value_type) * 3;
      }

      // Deserialize number of faces
      uint64_t numFaces;
      std::memcpy(&numFaces, ptr, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      // Resize and deserialize face data
      this->faces_.resize(numFaces);
      for (size_t i = 0; i < numFaces; ++i) {
        std::memcpy(&(this->faces_[i])[0], ptr, sizeof(int) * 3);
        ptr += sizeof(int) * 3;
      }

      // Return total bytes read including the size field

      return ptr - buffer;
    }

    // Helper method to calculate serialized size (updated)
    size_t getSerializedSize() const {
      size_t size = 0;
      size += sizeof(uint64_t); // Number of vertices
      size += vertices_.size() * sizeof(value_type) * 3; // vec_type data
      size += sizeof(uint64_t); // Number of faces
      size += faces_.size() * sizeof(int) * 3; // Face data
      return size;
    }


    bool buildBVH() {
      if (bvh.empty()) {
        rebuildBVH();
        return true;
      }

      return false;
    }

    void rebuildBVH() {
      bvh.reset();
      std::vector<FaceProxy> bvh_objects;
      bvh_objects.resize(faces_.size());
      for (size_t i = 0; i <faces_.size(); ++i) {
        auto& face=faces_[i];
        bvh_objects.emplace_back(*this, face);

      }
      bvh.build(bvh_objects);
    }
    bool hasBVH() const {
      return !bvh.empty();
    }

    size_t getBVHBoxes(std::vector<AABB<vec_type>>& output) const  {
      if (hasBVH()) {
        bvh.get_bboxes(
          output);

        return output.size();
      }
      return 0;

    }

    value_type distanceSq(const Vec &point) const {
      // Ensure BVH is built
      if (!bvh) {
        const_cast<Mesh*>(this)->rebuildBVH();
      }

      value_type bestDist = std::numeric_limits<value_type>::max();
      std::stack<const size_t> stack;
      stack.push(bvh.root_index);

      while (!stack.empty()) {

        const size_t node_index = stack.top();
        auto& node=bvh.nodes[node_index];
        stack.pop();

        // Compute the squared distance from the point to this node's bounding box
        value_type boxDist = distanceSqToAABB(point, node.bbox);
        if (boxDist >= bestDist) {
          // If the closest point in this bounding box is already worse than the best found so far,
          // skip this node (and its children).
          continue;
        }

        // If this is a leaf node, compute the actual distance to the face

        if (node.isLeaf()) {
          const size_t faceIndex = node.object;
          const Face &face = faces_[faceIndex];
          const Vec &v0 = vertices_[face[0]];
          const Vec &v1 = vertices_[face[1]];
          const Vec &v2 = vertices_[face[2]];

          Vec cp = closestPointOnTriangle(point, v0, v1, v2);
          value_type dist = (cp - point).sqLength();
          if (dist < bestDist) {
            bestDist = dist;
          }
        } else {
          // Internal node: push children
          if ((node.left!=-1)) stack.push(node.left);
          if ((node.right!=-1)) stack.push(node.right);
        }
      }

      return bestDist;
    }

    value_type distance(const Vec &point) const {
      // First ensure we have a BVH built

      return std::sqrt(distanceSq(point));
    }


    /**
     * @brief Determines if a point is inside a 3D mesh.
     *
     * This method checks whether a given point is inside the 3D mesh by casting rays
     * in different directions and counting intersections with mesh faces. It uses a BVH
     * (Bounding Volume Hierarchy) for acceleration and ensures that the point is sufficiently
     * away from the mesh surface based on a specified epsilon.
     *
     * @param point The 3D point to be tested for containment within the mesh.
     * @param eps The epsilon value used for collision detection tolerance, default is BVH_MESH_COLLISION_EPS.
     * @return Returns true if the point is inside the mesh, false otherwise.
     */
    bool pointInside(const Vec &point, const value_type eps = 1e-8) const {


      // Directions for casting rays
      Vec directions[2] = {{1.f, 0.f, 0.f}, {0.f, 0.f, 1.f}};
      size_t intersections[2] = {0, 0};

      // Compute the (unsigned) distance from the point to the mesh surface
      if (!bvh.bbox().inside(point)) {
        return false;
      }

      // For each direction, cast a ray and count intersections with mesh faces
      for (int ray = 0; ray < 2; ++ray) {
        // We'll use a stack to traverse the BVH
        std::stack<const size_t> stack;
        stack.push(bvh.root_index);

        // Ray origin and direction for this pass
        const Vec &ray_origin = point;
        const Vec &ray_dir = directions[ray];
        const Ray rayy={ray_origin, ray_dir};


        raycast_internal(rayy,bvh.root_index, intersections[ray],eps);

      }


      // After checking two directions, if both have odd intersections, the point is inside
      return ((intersections[0] % 2) == 1 && (intersections[1] % 2) == 1);

    }

    /**
     * @brief Computes the signed distance from a point to the mesh.
     *
     * This method calculates the signed distance between a given point and the mesh.
     * A negative distance indicates the point is inside the mesh, while a positive
     * distance indicates the point is outside. The calculation takes into account
     * a specified epsilon for collision detection.
     *
     * @param point The point for which the signed distance is calculated.
     * @return The signed distance from the point to the mesh. A negative value
     *         indicates the point is inside the mesh, and a positive value indicates
     *         it is outside.
     */
    value_type sd(const Vec &point) const {

      value_type dist=distance(point);
      if (pointInside(point, BVH_MESH_COLLISION_EPS)) {
        return -dist;
      }
      return dist;

    }

  };

  template<typename Vec>
  /**
   * @brief Serializes a vector of Mesh objects into a contiguous buffer.
   *
   * This function allocates a buffer containing the serialized data of all Mesh objects
   * in the provided vector. The serialization format is as follows:
   * - First 8 bytes: uint64_t representing the number of Mesh objects.
   * - For each Mesh:
   *   - 8 bytes: uint64_t representing the size of the serialized Mesh data (including the size prefix).
   *   - Serialized Mesh data.
   *
   * @param meshVector The vector of Mesh objects to serialize.
   * @param buffer A reference to a std::vector<char> that will hold the serialized data.
   *               The function will resize this vector to fit the serialized data.
   * @return size_t The total number of bytes written into the buffer.
   *
   * @throws std::runtime_error If serialization fails due to insufficient buffer space.
   */
  inline size_t serializeMeshVectorToBuffer(const std::vector<Mesh<Vec>> &meshVector, std::vector<char> &buffer) {
    // Step 1: Calculate the total buffer size required
    size_t totalSize = sizeof(uint64_t); // For numMeshes

    for (const auto &mesh: meshVector) {
      size_t meshDataSize = mesh.getSerializedSize();
      size_t meshSizeWithPrefix = sizeof(uint64_t) + meshDataSize; // size prefix + mesh data
      totalSize += meshSizeWithPrefix;
    }

    // Step 2: Allocate the buffer with the calculated size
    buffer.resize(totalSize);

    // Step 3: Write the number of Mesh objects
    char *ptr = buffer.data();
    uint64_t numMeshes = meshVector.size();
    std::memcpy(ptr, &numMeshes, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    // Step 4: Serialize each Mesh object
    for (const auto &mesh: meshVector) {
      // Get the serialized size of the mesh data
      size_t meshDataSize = mesh.getSerializedSize();
      uint64_t meshSizeWithPrefix = sizeof(uint64_t) + meshDataSize;

      // Write the mesh size prefix
      std::memcpy(ptr, &meshSizeWithPrefix, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      // Serialize the mesh data into the buffer
      mesh.serialize(ptr);
      ptr += meshDataSize;
    }

    return totalSize;
  }

  template<typename Vec>
  inline void deserializeMeshVectorFromBuffer(const char *buffer, std::vector<Mesh<Vec>> &meshVector) {
    const char *ptr = buffer;

    // Step 1: Read the number of Mesh objects
    uint64_t numMeshes = 0;
    std::memcpy(&numMeshes, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    // Resize the vector to hold all Mesh objects
    meshVector.resize(numMeshes);

    // Step 2: Deserialize each Mesh object
    for (uint64_t i = 0; i < numMeshes; ++i) {
      // Read the size of the current Mesh data
      uint64_t meshSize = 0;
      std::memcpy(&meshSize, ptr, sizeof(uint64_t));
      ptr += sizeof(uint64_t);

      // Allocate a temporary buffer to hold the Mesh data (including size prefix)
      std::vector<char> meshBuffer(meshSize + sizeof(uint64_t));

      // Copy the size prefix into the temporary buffer
      std::memcpy(meshBuffer.data(), &meshSize, sizeof(uint64_t));

      // Copy the Mesh data into the temporary buffer
      std::memcpy(meshBuffer.data() + sizeof(uint64_t), ptr, meshSize);
      ptr += meshSize;

      // Step 3: Deserialize the Mesh object from the temporary buffer
      meshVector[i].deserialize(meshBuffer.data());
    }
  }

}


#endif //MESH_H


