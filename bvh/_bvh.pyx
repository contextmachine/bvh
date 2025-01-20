# cython: language_level=3
# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8
cimport cython
from libc.stddef cimport size_t
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.math cimport modf
cimport numpy as cnp
import numpy as np
cnp.import_array()

cdef extern from "vec.h" namespace "bvh" nogil:
    cdef cppclass vec3d:
        double x
        double y
        double z

cdef extern from "aabb.h" namespace "bvh" nogil:
    cdef cppclass AABB[Vec]:
        Vec min
        Vec max
    cdef cppclass AABB3d:
        vec3d min
        vec3d max

cdef extern from "prims.h" namespace "bvh" nogil:
    cdef cppclass Segm[Vec]:
        Vec start
        Vec end
    cdef cppclass Ray[Vec]:
        Vec start
        Vec direction

    cdef cppclass Tri[Vec]:
        Vec a
        Vec b
        Vec c

cdef extern from "bvh.h" namespace "bvh" nogil:
    cdef cppclass BVH:
        BVH() except +
        bint empty() const
        void build(const vector[Tri[vec3d]] &primitives)

cdef extern from "raycast.h" namespace "bvh" nogil:
    void raycast(const Ray[vec3d] &inp, const BVH &bvh, const vector[Tri[vec3d]] &primitives, size_t &counts)
    void raycast(const vector[Ray[vec3d]] &rays,const BVH &bvh,const  vector[Tri[vec3d]] &primitives, vector[size_t] &counts)
    void raycast(const vector[Ray[vec3d]] &rays,const BVH &bvh,const  vector[Tri[vec3d]] &primitives, vector[pair[double,vec3d]] &hits, vector[size_t] &counts)

cdef class TriangleSoup:
    cdef vector[Tri[vec3d]] primitives
    cdef BVH bvh
    def __init__(self, double[:,:,:] triangles):
        

        cdef size_t i;
        cdef size_t n = triangles.shape[0]
     
        self.primitives=vector[Tri[vec3d]](n)
        self.bvh=BVH()

        for i in range(n):
            self.primitives[i].a.x=triangles[i,0,0]
            self.primitives[i].a.y=triangles[i,0,1]
            self.primitives[i].a.z=triangles[i,0,2]
            self.primitives[i].b.x=triangles[i,1,0]
            self.primitives[i].b.y=triangles[i,1,1]
            self.primitives[i].b.z=triangles[i,1,2]
            self.primitives[i].c.x=triangles[i,2,0]
            self.primitives[i].c.y=triangles[i,2,1]
            self.primitives[i].c.z=triangles[i,2,2]
        
    cdef bint _is_bvh_empty(self):
        return self.bvh.empty()
    def has_bvh(self):
        cdef bint res=self._is_bvh_empty()
        return not res

    def build_bvh(self):
        if self._is_bvh_empty():
            self.bvh.build(self.primitives)
        else:
            raise ValueError("BVH already built")
    
    cdef void raycast_counts(self, double[:,:,:] rays, size_t[:] counts) noexcept nogil:
        cdef size_t i;
        cdef size_t n = rays.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef vector[size_t] counts_cpp=vector[size_t]()
        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]
        
        raycast(rays_cpp,self.bvh,self.primitives,counts_cpp)
        for i in range(n):
            counts[i]=counts_cpp[i]

    cdef void raycast_hits(self, double[:,:,:] rays, double[:,:] hits, size_t[:] counts) noexcept nogil:
        cdef size_t i;
        cdef size_t n = rays.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]
        cdef vector[pair[double,vec3d]] hits_cpp=vector[pair[double,vec3d]]()
        cdef vector[size_t] counts_cpp=vector[size_t]()
        raycast(rays_cpp,self.bvh,self.primitives,hits_cpp,counts_cpp)
        for i in range(n):
            counts[i]=counts_cpp[i]
            hits[i,0]=hits_cpp[i].second.x
            hits[i,1]=hits_cpp[i].second.y
            hits[i,2]=hits_cpp[i].second.z
            hits[i,3]=hits_cpp[i].first
        


    def points_inside(self, double[:,:] points):
        cdef size_t i;
        cdef size_t n = points.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef bint[:] inside=np.zeros((n,),int)


        for i in range(n):
            rays_cpp[i].start.x=points[i,0]
            rays_cpp[i].start.y=points[i,1]
            rays_cpp[i].start.z=points[i,2]
            rays_cpp[i].direction.x=1.
            rays_cpp[i].direction.y=1.
            rays_cpp[i].direction.z=1.
        cdef vector[size_t] counts_cpp=vector[size_t]()
        raycast(rays_cpp,self.bvh,self.primitives,counts_cpp)

        for i in range(n):
            inside[i]=counts_cpp[i]%2
        return inside

        
     
        