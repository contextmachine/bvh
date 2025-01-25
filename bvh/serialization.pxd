# cython: language_level=3
# distutils: language = c++
cimport cython


from libcpp.vector cimport vector

cdef extern from "vec.h" namespace "bvh" nogil:
    cdef cppclass vec3d:
        double x, y, z

cdef extern from "serialization.h" namespace "bvh" nogil:
   
    # Declare bvh namespace and classes minimally

    cdef cppclass AABB[Vec]:
        Vec min, max
        void updateCentroid()


    cdef cppclass Ray[Vec]:
        Vec start, direction

    cdef cppclass Segm[Vec]:
        Vec start, end

    cdef cppclass Tri[Vec]:
        Vec a, b, c



    # Template serialize/deserialize
    cdef char* serialize_aabb[Vec](const vector[AABB[Vec]]& items, size_t& bytes_written) except +
    cdef void deserialize_aabb[Vec](const char* buffer, vector[AABB[Vec]]& items) except +


    cdef char* serialize_tri[Vec](const vector[Tri[Vec]]& items, size_t& bytes_written) except +
    cdef void deserialize_tri[Vec](const char* buffer, vector[Tri[Vec]]& items) except +

    cdef char* serialize_ray[Vec](const vector[Ray[Vec]]& items, size_t& bytes_written) except +
    cdef void deserialize_ray[Vec](const char* buffer, vector[Ray[Vec]]& items) except +

    cdef char* serialize_segm[Vec](const vector[Segm[Vec]]& items, size_t& bytes_written) except +
    cdef void deserialize_segm[Vec](const char* buffer, vector[Segm[Vec]]& items) except +

    cdef char* serialize_vec[Vec](const vector[Vec]& items, size_t& bytes_written) except +
    cdef void deserialize_vec[Vec](const char* buffer, vector[Vec]& items) except +