# cython: language_level=3
# distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8
cimport cython
from libc.stddef cimport size_t
from cython.operator cimport dereference as deref
from libc.stdlib cimport malloc,free,realloc
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.math cimport modf
from libcpp cimport bool
from libcpp.unordered_map cimport unordered_map
from cython.operator import dereference, postincrement
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
        bool inside(const Vec &pt) const
        bool insideStrict(const Vec &pt) const
        
        
    cdef cppclass AABB3d:
        vec3d min
        vec3d max
        bool inside(const vec3d &pt) const
        bool insideStrict(const vec3d &pt) const
        double sd( const vec3d &pt) const
        void sd( const vector[vec3d] &pts, vector[double] &sdf) const


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
        AABB3d bbox() const
        void build(const vector[Tri[vec3d]] &primitives)
        void build(const vector[AABB[vec3d]] &primitives)


cdef extern from "raycast.h" namespace "bvh" nogil:
    cdef cppclass hit3d:
        double x
        double y
        double z
        double w

    cdef cppclass Hit:
        double first
        vec3d second
        size_t prim
    cdef cppclass PrimHit:
        double t
        size_t ray
        vec3d pt
    cdef cppclass pair_ulong_hash:
        size_t operator()(const pair[size_t,size_t]& p) const 

    void raycast(const vector[Ray[vec3d]] &rays,const BVH &bvh,const  vector[Tri[vec3d]] &primitives, vector[size_t] &counts)
    void raycast(const vector[Ray[vec3d]] &rays,const BVH &bvh,const  vector[Tri[vec3d]] &primitives, vector[Hit] &hits, vector[size_t] &counts)
    void raycast_single_omp(const vector[Ray[vec3d]] &rays,const BVH &bvh,vector[Tri[vec3d]] &primitives,vector[bool] &mask)
    void raycast_bvh(const vector[Ray[vec3d]] &rays,const BVH &bvh,  unordered_map[pair[size_t,size_t],hit3d, pair_ulong_hash] &prims_rays_hits )

    

cdef extern from "reflection.h" namespace "bvh" nogil:
    void reflect(const vector[Ray[vec3d]] &rays,const BVH &bvh,const  vector[Tri[vec3d]] &primitives, vector[Ray[vec3d]]&reflected, vector[bool] &mask)



cdef class BVHTree:
    
    cdef vector[AABB[vec3d]] primitives
    cdef BVH bvh
    
    def __init__(self, double[:,:,:] boxes):
        

        cdef size_t i;
        cdef size_t n = boxes.shape[0]
     
        self.primitives=vector[AABB[vec3d]](n)
        self.bvh=BVH()

        for i in range(n):
            self.primitives[i].min.x=boxes[i,0,0]
            self.primitives[i].min.y=boxes[i,0,1]
            self.primitives[i].min.z=boxes[i,0,2]
            self.primitives[i].max.x=boxes[i,1,0]
            self.primitives[i].max.y=boxes[i,1,1]
            self.primitives[i].max.z=boxes[i,1,2]
    

    def leafs_count(self):
        return self.primitives.size()
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def get_leafs(self):
        cdef size_t n = self.primitives.size()
        cdef double[:,:,:] prims=np.empty((n,2,3))
        cdef size_t i;
        for i in range(n):
            prims[i,0,0]=self.primitives[i].min.x
            prims[i,0,1]=self.primitives[i].min.y
            prims[i,0,2]=self.primitives[i].min.z
            prims[i,1,0]=self.primitives[i].max.x
            prims[i,1,1]=self.primitives[i].max.y
            prims[i,1,2]=self.primitives[i].max.z
        return prims

    


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def bbox(self):
        cdef AABB3d bb
        cdef double[:,:] bb_arr
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        bb = self.bvh.bbox()
        bb_arr=np.empty((2,3))

        bb_arr[0,0]=bb.min.x
        bb_arr[0,1]=bb.min.y
        bb_arr[0,2]=bb.min.z
        bb_arr[1,0]=bb.max.x
        bb_arr[1,1]=bb.max.y
        bb_arr[1,2]=bb.max.z
        return bb_arr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def in_root_bbox(self, double[:,:] points, bool[:] result=None):
        cdef AABB3d bb
        cdef vec3d pt
        cdef size_t i;
        cdef size_t n = points.shape[0]
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        bb = self.bvh.bbox()
        if result is None:
            result=np.empty((n,),dtype=np.bool_)
        for i in range(n):

            pt.x=points[i,0]
            pt.y=points[i,1]
            pt.z=points[i,2]
            if not bb.inside(pt):
                result[i]=False
            else:
                result[i]=True
        return result


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def has_bvh(self):
        cdef bint res= self.bvh.empty()
        return not res
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def build_bvh(self):
        if  self.bvh.empty():
            self.bvh.build(self.primitives)
        else:
            raise ValueError("BVH already built")


   

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def raycast(self, double[:,:,:] rays) :
        cdef size_t i,row, col;
        cdef size_t cols=4
        
        cdef size_t n = rays.shape[0]
        cdef dict py_prims_rays_hits=dict()
        cdef tuple pykey;
   
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef unordered_map[pair[size_t,size_t],hit3d, pair_ulong_hash]  prims_rays_hits_cpp=unordered_map[pair[size_t,size_t],hit3d, pair_ulong_hash] ()

        
        
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")

        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]

        raycast_bvh(rays_cpp,self.bvh,prims_rays_hits_cpp)

        
        cdef unordered_map[ pair[size_t,size_t],hit3d, pair_ulong_hash].iterator it = prims_rays_hits_cpp.begin()
        cdef size_t x=0;
        cdef double[:, :] bhit3d=np.empty((prims_rays_hits_cpp.size(),4))
        
        while(it != prims_rays_hits_cpp.end()):
            # let's pretend here I just want to print the key and the value
            pykey=(dereference(it).first.first,      dereference(it).first.second)
            py_prims_rays_hits[pykey]=x
            bhit3d[x,0]=dereference(it).second.x
            bhit3d[x,1]=dereference(it).second.y
            bhit3d[x,2]=dereference(it).second.z
            bhit3d[x,3]=dereference(it).second.w
            
            postincrement(it) # Increment the iterator to the net element
            
            x+=1
        
        

        return py_prims_rays_hits, bhit3d



        

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

    def leafs_count(self):
        return self.primitives.size()
     
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def get_leafs(self):
        cdef size_t n = self.primitives.size()
        cdef double[:,:,:] prims=np.empty((n,3,3))
        cdef size_t i;
        for i in range(n):
            prims[i,0,0]=self.primitives[i].a.x
            prims[i,0,1]=self.primitives[i].a.y
            prims[i,0,2]=self.primitives[i].a.z
            prims[i,1,0]=self.primitives[i].b.x
            prims[i,1,1]=self.primitives[i].b.y
            prims[i,1,2]=self.primitives[i].b.z
            prims[i,2,0]=self.primitives[i].c.x
            prims[i,2,1]=self.primitives[i].c.y
            prims[i,2,2]=self.primitives[i].c.z
        return prims

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def bbox(self):
        cdef AABB3d bb
        cdef double[:,:] bb_arr
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        bb = self.bvh.bbox()
        bb_arr=np.empty((2,3))

        bb_arr[0,0]=bb.min.x
        bb_arr[0,1]=bb.min.y
        bb_arr[0,2]=bb.min.z
        bb_arr[1,0]=bb.max.x
        bb_arr[1,1]=bb.max.y
        bb_arr[1,2]=bb.max.z
        return bb_arr
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def in_root_bbox(self, double[:,:] points, bool[:] result=None):
        cdef AABB3d bb
        cdef vec3d pt
        cdef size_t i;
        cdef size_t n = points.shape[0]
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        bb = self.bvh.bbox()
        if result is None:
            result=np.empty((n,),dtype=np.bool_)
        for i in range(n):

            pt.x=points[i,0]
            pt.y=points[i,1]
            pt.z=points[i,2]
            if not bb.inside(pt):
                result[i]=False
            else:
                result[i]=True
        return result


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.initializedcheck(False)
    @cython.nonecheck(True)
    def has_bvh(self):
        cdef bint res= self.bvh.empty()
        return not res
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def build_bvh(self):
        if  self.bvh.empty():
            self.bvh.build(self.primitives)
        else:
            raise ValueError("BVH already built")


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def raycast_single(self, double[:,:,:] rays, bint[:] mask=None):
        cdef size_t i;
        cdef size_t n = rays.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef vector[bool] mask_cpp=vector[bool]()
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        if mask is None:
            mask=np.empty((n,),dtype=np.bool_)
        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]
        raycast_single_omp(rays_cpp, self.bvh, self.primitives, mask_cpp)
        for i in range(n):
          mask[i]=mask_cpp[i]
        return mask

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def reflection(self, double[:,:,:] rays, double[:,:,:] result=None, bool[:] mask=None):
        cdef size_t i;
        cdef size_t n = rays.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef vector[Ray[vec3d]] reflected_cpp=vector[Ray[vec3d]]()
        cdef vector[bool] mask_cpp=vector[bool]()
        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]
        if mask is None:
            mask=np.empty((n,),dtype=np.bool_)
        if result is None:
            result=np.empty(rays.shape)
        reflect(rays_cpp,self.bvh,self.primitives,reflected_cpp,mask_cpp)
        for i in range(n):
            result[i,0,0]=reflected_cpp[i].start.x
            result[i,0,1]=reflected_cpp[i].start.y
            result[i,0,2]=reflected_cpp[i].start.z
            result[i,1,0]=reflected_cpp[i].direction.x
            result[i,1,1]=reflected_cpp[i].direction.y
            result[i,1,2]=reflected_cpp[i].direction.z
            mask[i]=mask_cpp[i]

        return result, mask


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def raycast_counts(self, double[:,:,:] rays, size_t[:] counts=None):
        cdef size_t i;
        cdef size_t n = rays.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef vector[size_t] counts_cpp=vector[size_t]()
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        if counts is None:
            counts=np.empty((n,),dtype=np.uint)
        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]
        
        raycast(rays_cpp, self.bvh, self.primitives, counts_cpp)
        for i in range(n):
            counts[i]=counts_cpp[i]
        return counts

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def raycast_hits(self, double[:,:,:] rays) :
        cdef size_t i,row, col;
        cdef size_t cols=4

        cdef size_t n = rays.shape[0]
        cdef size_t[:] counts=np.empty((n,),dtype=np.uint)
        cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
        cdef vector[Hit] hits_cpp=vector[Hit]()
        cdef vector[size_t] counts_cpp=vector[size_t]()
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")

        for i in range(n):
            rays_cpp[i].start.x=rays[i,0,0]
            rays_cpp[i].start.y=rays[i,0,1]
            rays_cpp[i].start.z=rays[i,0,2]
            rays_cpp[i].direction.x=rays[i,1,0]
            rays_cpp[i].direction.y=rays[i,1,1]
            rays_cpp[i].direction.z=rays[i,1,2]

        raycast(rays_cpp,self.bvh,self.primitives,hits_cpp,counts_cpp)
        cdef double[:,:] hits=np.empty((hits_cpp.size(),4))
        cdef size_t[:] prims=np.empty((hits_cpp.size(),),dtype=np.uint)
        for i in range(n):
            counts[i]=counts_cpp[i]

        for i in range(hits_cpp.size()):
            hits[i,0]=hits_cpp[i].second.x
            hits[i,1]=hits_cpp[i].second.y
            hits[i,2]=hits_cpp[i].second.z
            hits[i,3]=hits_cpp[i].first
            prims[i]=hits_cpp[i].prim
        return hits,counts,prims


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def points_inside(self, double[:,:] points, bool[:] inside=None):
        cdef vec3d pt
        cdef size_t i;
        cdef size_t n = points.shape[0]
        cdef vector[Ray[vec3d]] rays_cpp
        cdef vector[size_t] counts_cpp=vector[size_t]()
        cdef bool is_inside
        if self.bvh.empty():
            raise ValueError("BVH is empty. Build it first.")
        rays_cpp=vector[Ray[vec3d]](n)
        if inside is None:
            inside=np.zeros((n,),np.bool_)
        
        for i in range(n):

            rays_cpp[i].start.x=points[i,0]
            rays_cpp[i].start.y=points[i,1]
            rays_cpp[i].start.z=points[i,2]
            rays_cpp[i].direction.x=1.
            rays_cpp[i].direction.y=1.
            rays_cpp[i].direction.z=1.

        raycast(rays_cpp,self.bvh,self.primitives,counts_cpp)

        for i in range(n):

            inside[i]=(counts_cpp[i]%2)
        return inside
    """
    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sd(self, double[:,:] points, double[:] result=None):
        cdef vector[vec3d] pts=vector[vec3d](points.shape[0])
        cdef vector[double] result_cpp=vector[double](points.shape[0])
        cdef size_t i;
        cdef size_t n = points.shape[0]
        if result is None:
            result=np.empty((n,),dtype=np.double)
        for i in range(n):
            pts[i].x=points[i,0]
            pts[i].y=points[i,1]
            pts[i].z=points[i,2]
        
        self.bvh.sd(pts,result_cpp)
        raycast_first(rays_cpp,self.bvh,self.primitives,counts_cpp)
        memcpy(result,result_cpp.data(),n*sizeof(double))
        
        return result
    """


        
   