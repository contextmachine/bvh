
cimport cython
cimport bvh.serialization

from libcpp.string cimport string

cimport numpy as cnp
import numpy as np
cnp.import_array()


def serialize_aabb3d(double[:,:,:] boxes):

    cdef size_t bytes_written = 0
    cdef size_t n = boxes.shape[0]
    cdef size_t i;
    cdef vector[AABB[vec3d]] aabbs=vector[AABB[vec3d]](n)
    for i in range(n):
        aabbs[i].min.x = boxes[i,0,0]
        aabbs[i].min.y = boxes[i,0,1]
        aabbs[i].min.z = boxes[i,0,2]
        aabbs[i].max.x = boxes[i,1,0]
        aabbs[i].max.y = boxes[i,1,1]
        aabbs[i].max.z = boxes[i,1,2]
        aabbs[i].updateCentroid()


    cdef char* buffer = serialize_aabb[vec3d](aabbs, bytes_written)
    cdef bytes result = buffer[:bytes_written]
    return result


def deserialize_aabb3d(const unsigned char[:] buffer, double[:,:,:] boxes=None):
    cdef vector[AABB[vec3d]] aabbs;
    
    cdef const char* buf = <const char*> (&buffer[0])

    
    deserialize_aabb[vec3d](buf, aabbs)
    cdef size_t n = aabbs.size()
    cdef size_t i;
    if boxes is None:
        boxes = np.zeros((n,2,3))
    for i in range(n):
        boxes[i,0,0]=aabbs[i].min.x 
        boxes[i,0,1]=aabbs[i].min.y 
        boxes[i,0,2]=aabbs[i].min.z 
        boxes[i,1,0]=aabbs[i].max.x 
        boxes[i,1,1]=aabbs[i].max.y 
        boxes[i,1,2]=aabbs[i].max.z 


 
    return boxes

def serialize_tri3d(double[:,:,:] tris):

    cdef size_t bytes_written = 0
    cdef size_t n = tris.shape[0]
    cdef size_t i;
    cdef vector[Tri[vec3d]] tris_cpp=vector[Tri[vec3d]](n)
    for i in range(n):
        tris_cpp[i].a.x = tris[i,0,0]
        tris_cpp[i].a.y = tris[i,0,1]
        tris_cpp[i].a.z = tris[i,0,2]
        tris_cpp[i].b.x = tris[i,1,0]
        tris_cpp[i].b.y = tris[i,1,1]
        tris_cpp[i].b.z = tris[i,1,2]
        tris_cpp[i].c.x = tris[i,2,0]
        tris_cpp[i].c.y = tris[i,2,1]
        tris_cpp[i].c.z = tris[i,2,2]



    cdef char* buffer = serialize_tri[vec3d](tris_cpp, bytes_written)
    cdef bytes result = buffer[:bytes_written]
    return result


def deserialize_tri3d(const unsigned char[:] buffer, double[:,:,:] tris=None):
    cdef vector[Tri[vec3d]] tris_cpp    ;

    cdef const char* buf = <const char*> (&buffer[0])


    deserialize_tri[vec3d](buf, tris_cpp)
    cdef size_t n = tris_cpp.size()
    cdef size_t i;
    if tris is None:
        tris = np.zeros((n,3,3))
    for i in range(n):
        tris[i,0,0]=tris_cpp[i].a.x
        tris[i,0,1]=tris_cpp[i].a.y
        tris[i,0,2]=tris_cpp[i].a.z
        tris[i,1,0]=tris_cpp[i].b.x
        tris[i,1,1]=tris_cpp[i].b.y
        tris[i,1,2]=tris_cpp[i].b.z
        tris[i,2,0]=tris_cpp[i].c.x
        tris[i,2,1]=tris_cpp[i].c.y
        tris[i,2,2]=tris_cpp[i].c.z



    return tris


def serialize_ray3d(double[:,:,:] rays):

    cdef size_t bytes_written = 0
    cdef size_t n = rays.shape[0]
    cdef size_t i;
    cdef vector[Ray[vec3d]] rays_cpp=vector[Ray[vec3d]](n)
    for i in range(n):
        rays_cpp[i].start.x = rays[i,0,0]
        rays_cpp[i].start.y = rays[i,0,1]
        rays_cpp[i].start.z = rays[i,0,2]
        rays_cpp[i].direction.x = rays[i,1,0]
        rays_cpp[i].direction.y = rays[i,1,1]
        rays_cpp[i].direction.z = rays[i,1,2]
 


    cdef char* buffer = serialize_ray[vec3d](rays_cpp, bytes_written)
    cdef bytes result = buffer[:bytes_written]
    return result

def deserialize_ray3d(const unsigned char[:] buffer, double[:,:,:] rays=None):
    cdef vector[Ray[vec3d]] rays_cpp;

    cdef const char* buf = <const char*> (&buffer[0])


    deserialize_ray[vec3d](buf, rays_cpp)
    cdef size_t n = rays_cpp.size()
    cdef size_t i;
    if rays is None:
        rays = np.zeros((n,2,3))
    for i in range(n):
        rays[i,0,0]=rays_cpp[i].start.x
        rays[i,0,1]=rays_cpp[i].start.y
        rays[i,0,2]=rays_cpp[i].start.z
        rays[i,1,0]=rays_cpp[i].direction.x
        rays[i,1,1]=rays_cpp[i].direction.y
        rays[i,1,2]=rays_cpp[i].direction.z



    return rays

def serialize_seg3d(double[:,:,:] segs):
    cdef size_t bytes_written = 0
    cdef size_t n = segs.shape[0]
    cdef size_t i;
    cdef vector[Segm[vec3d]] segs_cpp=vector[Segm[vec3d]](n)
    for i in range(n):
        segs_cpp[i].start.x = segs[i,0,0]
        segs_cpp[i].start.y = segs[i,0,1]
        segs_cpp[i].start.z = segs[i,0,2]
        segs_cpp[i].end.x = segs[i,1,0]
        segs_cpp[i].end.y = segs[i,1,1]
        segs_cpp[i].end.z = segs[i,1,2]


    cdef char* buffer = serialize_segm[vec3d](segs_cpp, bytes_written)
    cdef bytes result = buffer[:bytes_written]
    return result

def deserialize_seg3d(const unsigned char[:] buffer, double[:,:,:] segs=None):
    cdef vector[Segm[vec3d]] segs_cpp;

    cdef const char* buf = <const char*> (&buffer[0])


    deserialize_segm[vec3d](buf, segs_cpp)
    cdef size_t n = segs_cpp.size()
    cdef size_t i;
    if segs is None:
        segs = np.zeros((n,2,3))
    for i in range(n):
        segs[i,0,0]=segs_cpp[i].start.x
        segs[i,0,1]=segs_cpp[i].start.y
        segs[i,0,2]=segs_cpp[i].start.z
        segs[i,1,0]=segs_cpp[i].end.x
        segs[i,1,1]=segs_cpp[i].end.y
        segs[i,1,2]=segs_cpp[i].end.z



    return segs
def serialize_vec3d(double[:,:] vecs):
    cdef size_t bytes_written = 0
    cdef size_t n = vecs.shape[0]
    cdef size_t i;
    cdef vector[vec3d] vec=vector[vec3d](n)
    for i in range(n):
        vec[i].x = vecs[i,0]
        vec[i].y = vecs[i,1]
        vec[i].z = vecs[i,2]

    cdef char* buffer = serialize_vec[vec3d](vec, bytes_written)
    cdef bytes result = buffer[:bytes_written]
    return result

def deserialize_vec3d(const unsigned char[:] buffer, double[:,:] vecs=None):
    cdef vector[vec3d] vec;

    cdef const char* buf = <const char*> (&buffer[0])


    deserialize_vec[vec3d](buf, vec)
    cdef size_t n = vec.size()
    cdef size_t i;
    if vecs is None:
        vecs = np.zeros((n,3))
    for i in range(n):
        vecs[i,0]=vec[i].x
        vecs[i,1]=vec[i].y
        vecs[i,2]=vec[i].z



    return vecs

def deserialize(const unsigned char[:] buffer, int typecode, size_t dim=3):
    if dim != 3:
        raise NotImplementedError("Only 3D arrays are supported")
    if typecode == 1:
        return deserialize_vec3d(buffer)
    elif typecode == 2:
        return deserialize_aabb3d(buffer)
    elif typecode == 3:
        return deserialize_ray3d(buffer)
    elif typecode == 4:
        return deserialize_seg3d(buffer)
    elif typecode == 5:
        return deserialize_tri3d(buffer)

def serialize(arr, int typecode, size_t dim=3):
    """
        TYPE_VEC  = 1,   // For std::vector< Vec >
        TYPE_AABB = 2,   // For std::vector< AABB<Vec> >
        TYPE_RAY  = 3,   // For std::vector< Ray<Vec> >
        TYPE_SEGM = 4,   // For std::vector< Segm<Vec> >
        TYPE_TRI  = 5    // For std::vector< Tri<Vec> >"""
    if dim != 3:
        raise NotImplementedError("Only 3D arrays are supported")
    if typecode == 1:
        return serialize_vec3d(arr)
    elif typecode == 2:
        return serialize_aabb3d(arr)
    elif typecode == 3:
        return serialize_tri3d(arr)
    elif typecode == 4:
        return serialize_ray3d(arr)
    elif typecode == 5:
        return serialize_seg3d(arr)
