# test_vec
test_vec = executable('test_vec',
               'test_vec.cpp',
               include_directories : my_inc,
               install : false
)
test('test_vec',
     test_vec)

# test_bvh
test_bvh = executable('test_bvh',
                       'test_bvh.cpp',
                       include_directories : my_inc,
                       install : false
)
test('test_bvh',
     test_bvh)

# test_aabb
test_aabb = executable('test_aabb',
                       'test_aabb.cpp',
                       include_directories : my_inc,
                       install : false
)
test('test_aabb',
     test_aabb)

# test_raycast
test_raycast = executable('test_raycast',
               'test_raycast.cpp',
               include_directories : my_inc,
               install : false,
     dependencies: [openmp_dep],
          link_args:link_args


)
test('test_raycast',
     test_raycast)

# test_raycast_eigen
#test_raycast_eigen = executable('test_raycast_eigen',
#                                'test_raycast_eigen.cpp',
#                                include_directories : my_inc,
#                                install : false,
#dependencies: [openmp_dep,eigen_dep],
#)
#test('test_raycast_eigen',
#     test_raycast_eigen)

# test_mesh
test_mesh = executable('test_mesh',
               'test_mesh.cpp',
               include_directories : my_inc,
               install : false
)
test('test_mesh',
     test_mesh)
