
/**
 * @file serialization.h
 * @brief Complete implementations of serialize/deserialize functions for
 *        bvh::vec, AABB, Ray, Segm, and Tri types.
 *
 * The encoding format used for all types is:
 *   1) A 32-bit unsigned integer type code (unique per type).
 *   2) A 32-bit unsigned integer for the dim (Vec::dim).
 *   3) A 64-bit unsigned integer for the number of elements in the vector.
 *   4) The raw data for each element.
 *
 * For each element type, the raw data layout is:
 *   - vec<T,N>:        N scalars
 *   - AABB<vec<T,N>>:  2*N scalars (min[], max[])
 *   - Ray<vec<T,N>>:   2*N scalars (start[], direction[])
 *   - Segm<vec<T,N>>:  2*N scalars (start[], end[])
 *   - Tri<vec<T,N>>:   3*N scalars (a[], b[], c[])
 *
 * The deserialize functions read back these prefixes and data. If the type code or
 * dim do not match what is expected for that template instantiation, a
 * std::runtime_error is thrown. The deserialize function resizes the output vector
 * to the decoded length.
 *
 *
 */

#ifndef SERIALIZATION_H
#define SERIALIZATION_H


#include "aabb.h"
#include "prims.h"
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstdlib>
#include "typecodes.h"
#include <string.h>
namespace bvh {



/**********************************************************************
 *  Helper: Write/Read primitives from the buffer
 **********************************************************************/
inline void writePodToBuffer(char*& dest, const void* src, std::size_t size)
{
    std::memcpy(dest, src, size);
    dest += size;
}

inline void readPodFromBuffer(const char*& src, void* dest, std::size_t size)
{
    std::memcpy(dest, src, size);
    src += size;
}

/**********************************************************************
 *  Helper: ReadPod from a double pointer (needed for Tri's signature)
 **********************************************************************/
inline void readPodFromBuffer(const char** src, void* dest, std::size_t size)
{
    std::memcpy(dest, *src, size);
    *src += size;
}


template<typename Vec>
inline char* serialize_vec(const std::vector<Vec> &items, size_t& bytes_written)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t type_id  = TYPE_VEC;
    constexpr std::uint32_t dim      = static_cast<std::uint32_t>(Vec::dim);
    const std::uint64_t t_size       = sizeof(T);
    const std::uint64_t count        = items.size();
    // Each Vec has "dim" components of type T
    const std::uint64_t components_per_item = dim;
    const std::uint64_t total_data_bytes    = count * components_per_item * t_size;

    // The header has: type_id(4 bytes), dim(4 bytes), t_size(8 bytes), count(8 bytes)
    const size_t header_size = sizeof(type_id)
                             + sizeof(dim)
                             + sizeof(t_size)
                             + sizeof(count);
    const size_t total_size  = header_size + total_data_bytes;
    bytes_written = total_size;

    char* buffer = new char[total_size];
    char* p      = buffer;

    // 1) Write the header
    writePodToBuffer(p, &type_id,  sizeof(type_id));
    writePodToBuffer(p, &dim,      sizeof(dim));
    writePodToBuffer(p, &t_size,   sizeof(t_size));
    writePodToBuffer(p, &count,    sizeof(count));

    // 2) Write the data
    for (std::size_t i = 0; i < count; ++i) {
        // We need to copy each component in Vec
        // Typically, we'd do items[i][0], items[i][1], etc.
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i][c];
            writePodToBuffer(p, &val, sizeof(T));
        }
    }
    return buffer;
}

template<typename Vec>
inline void deserialize_vec(const char* buffer, std::vector<Vec> &items)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t expected_type_id = TYPE_VEC;
    constexpr std::uint32_t expected_dim     = Vec::dim;
    const char* p = buffer;

    // 1) Read the header
    std::uint32_t type_id;
    std::uint32_t dim;
    std::uint64_t t_size;
    std::uint64_t count;

    readPodFromBuffer(p, &type_id,  sizeof(type_id));
    readPodFromBuffer(p, &dim,      sizeof(dim));
    readPodFromBuffer(p, &t_size,   sizeof(t_size));
    readPodFromBuffer(p, &count,    sizeof(count));

    // 2) Check
    if (type_id != expected_type_id) {
        throw std::runtime_error("deserialize<Vec>: type_id mismatch");
    }
    if (dim != expected_dim) {
        throw std::runtime_error("deserialize<Vec>: dim mismatch");
    }
    if (t_size != sizeof(T)) {
        throw std::runtime_error("deserialize<Vec>: scalar size mismatch");
    }

    // 3) Read items
    items.resize(count);
    for(std::size_t i = 0; i < count; ++i) {
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i][c] = val;
        }
    }
}




/**********************************************************************
 *  2) Serialize/Deserialize AABB
 *     Each AABB<Vec> stores: min, max, centroid (3 * dim values total)
 **********************************************************************/



template<typename Vec>
inline char* serialize_aabb(const std::vector<AABB<Vec>> &items, size_t& bytes_written)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t type_id  = TYPE_AABB;
    constexpr std::uint32_t dim      = static_cast<std::uint32_t>(Vec::dim);
    const std::uint64_t t_size       = sizeof(T);
    const std::uint64_t count        = items.size();

    // AABB has: min (dim), max (dim), centroid (dim) => 3*dim
    const std::uint64_t components_per_item = 3 * dim;
    const std::uint64_t total_data_bytes    = count * components_per_item * t_size;

    // Header
    const size_t header_size = sizeof(type_id)
                             + sizeof(dim)
                             + sizeof(t_size)
                             + sizeof(count);
    const size_t total_size  = header_size + total_data_bytes;
    bytes_written = total_size;

    char* buffer = new char[total_size];
    char* p      = buffer;

    // Write header
    writePodToBuffer(p, &type_id, sizeof(type_id));
    writePodToBuffer(p, &dim,     sizeof(dim));
    writePodToBuffer(p, &t_size,  sizeof(t_size));
    writePodToBuffer(p, &count,   sizeof(count));

    // Write data
    // each item => min[dim] + max[dim] + centroid[dim]
    for (std::size_t i = 0; i < count; ++i) {
        // min
        for (std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].min[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
        // max
        for (std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].max[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
        // centroid
        for (std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].centroid[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
    }
    return buffer;
}

template<typename Vec>
inline void deserialize_aabb(const char* buffer, std::vector<AABB<Vec>> &items)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t expected_type_id = TYPE_AABB;
    constexpr std::uint32_t expected_dim     = Vec::dim;
    const char* p = buffer;

    // Read header
    std::uint32_t type_id;
    std::uint32_t dim;
    std::uint64_t t_size;
    std::uint64_t count;

    readPodFromBuffer(p, &type_id, sizeof(type_id));
    readPodFromBuffer(p, &dim,     sizeof(dim));
    readPodFromBuffer(p, &t_size,  sizeof(t_size));
    readPodFromBuffer(p, &count,   sizeof(count));

    // Check
    if (type_id != expected_type_id) {
        throw std::runtime_error("deserialize<AABB<Vec>>: type_id mismatch");
    }
    if (dim != expected_dim) {
        throw std::runtime_error("deserialize<AABB<Vec>>: dim mismatch");
    }
    if (t_size != sizeof(T)) {
        throw std::runtime_error("deserialize<AABB<Vec>>: scalar size mismatch");
    }

    // Read items
    items.resize(count);
    for (std::size_t i = 0; i < count; ++i) {
        // min
        for (std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].min[c] = val;
        }
        // max
        for (std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].max[c] = val;
        }
        // centroid
        for (std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].centroid[c] = val;
        }
    }
}

/**********************************************************************
 *  3) Serialize/Deserialize Ray
 *     Each Ray<Vec> = start (dim), direction (dim) => 2*dim
 **********************************************************************/
template<typename Vec>
inline char* serialize_ray(const std::vector<bvh::Ray<Vec>> &items, size_t& bytes_written)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t type_id  = TYPE_RAY;
    constexpr std::uint32_t dim      = static_cast<std::uint32_t>(Vec::dim);
    const std::uint64_t t_size       = sizeof(T);
    const std::uint64_t count        = items.size();

    // Ray has 2*dim components
    const std::uint64_t components_per_item = 2 * dim;
    const std::uint64_t total_data_bytes    = count * components_per_item * t_size;

    // Header
    const size_t header_size = sizeof(type_id)
                             + sizeof(dim)
                             + sizeof(t_size)
                             + sizeof(count);
    const size_t total_size  = header_size + total_data_bytes;
    bytes_written = total_size;

    char* buffer = new char[total_size];
    char* p      = buffer;

    // Write header
    writePodToBuffer(p, &type_id, sizeof(type_id));
    writePodToBuffer(p, &dim,     sizeof(dim));
    writePodToBuffer(p, &t_size,  sizeof(t_size));
    writePodToBuffer(p, &count,   sizeof(count));

    // Write data: start + direction
    for (std::size_t i = 0; i < count; ++i) {
        // start
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].start[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
        // direction
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].direction[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
    }
    return buffer;
}

template<typename Vec>
inline void deserialize_ray(const char* buffer, std::vector<bvh::Ray<Vec>> &items)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t expected_type_id = TYPE_RAY;
    constexpr std::uint32_t expected_dim     = Vec::dim;
    const char* p = buffer;

    // Read header
    std::uint32_t type_id;
    std::uint32_t dim;
    std::uint64_t t_size;
    std::uint64_t count;

    readPodFromBuffer(p, &type_id, sizeof(type_id));
    readPodFromBuffer(p, &dim,     sizeof(dim));
    readPodFromBuffer(p, &t_size,  sizeof(t_size));
    readPodFromBuffer(p, &count,   sizeof(count));

    // Check
    if (type_id != expected_type_id) {
        throw std::runtime_error("deserialize<Ray<Vec>>: type_id mismatch");
    }
    if (dim != expected_dim) {
        throw std::runtime_error("deserialize<Ray<Vec>>: dim mismatch");
    }
    if (t_size != sizeof(T)) {
        throw std::runtime_error("deserialize<Ray<Vec>>: scalar size mismatch");
    }

    // Read items
    items.resize(count);
    for(std::size_t i = 0; i < count; ++i) {
        // start
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].start[c] = val;
        }
        // direction
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].direction[c] = val;
        }
    }
}

/**********************************************************************
 *  4) Serialize/Deserialize Segm
 *     Each Segm<Vec> = start (dim), end (dim) => 2*dim
 **********************************************************************/
template<typename Vec>
inline char* serialize_segm(const std::vector<bvh::Segm<Vec>> &items, size_t& bytes_written)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t type_id  = TYPE_SEGM;
    constexpr std::uint32_t dim      = static_cast<std::uint32_t>(Vec::dim);
    const std::uint64_t t_size       = sizeof(T);
    const std::uint64_t count        = items.size();

    // Segm has 2*dim components
    const std::uint64_t components_per_item = 2 * dim;
    const std::uint64_t total_data_bytes    = count * components_per_item * t_size;

    // Header
    const size_t header_size = sizeof(type_id)
                             + sizeof(dim)
                             + sizeof(t_size)
                             + sizeof(count);
    const size_t total_size  = header_size + total_data_bytes;
    bytes_written = total_size;

    char* buffer = new char[total_size];
    char* p      = buffer;

    // Write header
    writePodToBuffer(p, &type_id, sizeof(type_id));
    writePodToBuffer(p, &dim,     sizeof(dim));
    writePodToBuffer(p, &t_size,  sizeof(t_size));
    writePodToBuffer(p, &count,   sizeof(count));

    // Write data: start + end
    for (std::size_t i = 0; i < count; ++i) {
        // start
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].start[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
        // end
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].end[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
    }
    return buffer;
}

template<typename Vec>
inline void deserialize_segm(const char* buffer, std::vector<bvh::Segm<Vec>> &items)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t expected_type_id = TYPE_SEGM;
    constexpr std::uint32_t expected_dim     = Vec::dim;
    const char* p = buffer;

    // Read header
    std::uint32_t type_id;
    std::uint32_t dim;
    std::uint64_t t_size;
    std::uint64_t count;

    readPodFromBuffer(p, &type_id, sizeof(type_id));
    readPodFromBuffer(p, &dim,     sizeof(dim));
    readPodFromBuffer(p, &t_size,  sizeof(t_size));
    readPodFromBuffer(p, &count,   sizeof(count));

    // Check
    if (type_id != expected_type_id) {
        throw std::runtime_error("deserialize<Segm<Vec>>: type_id mismatch");
    }
    if (dim != expected_dim) {
        throw std::runtime_error("deserialize<Segm<Vec>>: dim mismatch");
    }
    if (t_size != sizeof(T)) {
        throw std::runtime_error("deserialize<Segm<Vec>>: scalar size mismatch");
    }

    // Read items
    items.resize(count);
    for(std::size_t i = 0; i < count; ++i) {
        // start
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].start[c] = val;
        }
        // end
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].end[c] = val;
        }
    }
}

/**********************************************************************
 *  5) Serialize/Deserialize Tri
 *     Each Tri<Vec> = a (dim), b (dim), c (dim) => 3*dim
 *
 *  Note the last signature:
 *      inline void deserialize(const char** buffer, std::vector<Tri<Vec>>&)
 *  is slightly different, so we implement that exactly as requested.
 **********************************************************************/
template<typename Vec>
inline char* serialize_tri(const std::vector<bvh::Tri<Vec>> &items, size_t& bytes_written)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t type_id  = TYPE_TRI;
    constexpr std::uint32_t dim      = static_cast<std::uint32_t>(Vec::dim);
    const std::uint64_t t_size       = sizeof(T);
    const std::uint64_t count        = items.size();

    // Tri has 3*dim components
    const std::uint64_t components_per_item = 3 * dim;
    const std::uint64_t total_data_bytes    = count * components_per_item * t_size;

    // Header
    const size_t header_size = sizeof(type_id)
                             + sizeof(dim)
                             + sizeof(t_size)
                             + sizeof(count);
    const size_t total_size  = header_size + total_data_bytes;
    bytes_written = total_size;

    char* buffer = new char[total_size];
    char* p      = buffer;

    // Write header
    writePodToBuffer(p, &type_id, sizeof(type_id));
    writePodToBuffer(p, &dim,     sizeof(dim));
    writePodToBuffer(p, &t_size,  sizeof(t_size));
    writePodToBuffer(p, &count,   sizeof(count));

    // Write data: a + b + c
    for (std::size_t i = 0; i < count; ++i) {
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].a[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].b[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
        for(std::size_t c = 0; c < dim; ++c) {
            const T val = items[i].c[c];
            writePodToBuffer(p, &val, sizeof(T));
        }
    }
    return buffer;
}

template<typename Vec>
inline void deserialize_tri(const char* buffer, std::vector<bvh::Tri<Vec>> &items)
{
    using T = typename Vec::value_type;
    constexpr std::uint32_t expected_type_id = TYPE_TRI;
    constexpr std::uint32_t expected_dim     = Vec::dim;
    const char* p = buffer; // convenience reference for reading

    // Read header
    std::uint32_t type_id;
    std::uint32_t dim;
    std::uint64_t t_size;
    std::uint64_t count;

    readPodFromBuffer(p, &type_id, sizeof(type_id));
    readPodFromBuffer(p, &dim,     sizeof(dim));
    readPodFromBuffer(p, &t_size,  sizeof(t_size));
    readPodFromBuffer(p, &count,   sizeof(count));

    // Check
    if (type_id != expected_type_id) {
        throw std::runtime_error("deserialize<Tri<Vec>>: type_id mismatch");
    }
    if (dim != expected_dim) {
        throw std::runtime_error("deserialize<Tri<Vec>>: dim mismatch");
    }
    if (t_size != sizeof(T)) {
        throw std::runtime_error("deserialize<Tri<Vec>>: scalar size mismatch");
    }

    // Read items
    items.resize(count);
    for(std::size_t i = 0; i < count; ++i) {
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].a[c] = val;
        }
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].b[c] = val;
        }
        for(std::size_t c = 0; c < dim; ++c) {
            T val;
            readPodFromBuffer(p, &val, sizeof(T));
            items[i].c[c] = val;
        }
    }
}



}


#endif //SERIALIZATION_H
