project('bvh', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++17','optimization=3'])
# Detect the host operating system
simd_args = ['-flto']
link_args=[ '-flto']
eigen_dep = dependency('eigen3', required: false)

if host_machine.system() == 'darwin'
    #simd_args = ['-march=armv8-a+simd']
     simd_args += ['-mcpu=apple-m1','-funroll-loops']#,#'-flto', ]
    openmp_dep = dependency('openmp', required: false)
        # Fallback for macOS where OpenMP is not part of Clang by default
        if not openmp_dep.found() and host_machine.system() == 'darwin'
          openmp_dep = cpp.find_library('omp', required: false)
          message('darwin dep omp:')
        endif
        # Ensure OpenMP is found on either platform
        if not openmp_dep.found()

          error('OpenMP not found! Install libomp on macOS or check your compiler setup on Linux.')
        endif

elif host_machine.cpu() == 'aarch64'
        simd_args += ['-mcpu=native','-march=armv8-a+simd',]


else
  simd_args += [
    '-msse', '-msse2', '-msse3', '-mssse3',
    '-msse4.1', '-msse4.2', '-mavx', '-mavx2'
  ]
endif



message('LINK ARGS:',link_args)
message('CXX ARGS',simd_args)
add_project_arguments(simd_args, language: 'cpp')
my_inc = include_directories(['.'])
#bvh = executable('bvh', 'main.cpp', install : false, include_directories : my_inc)
#test('test', bvh)

libbvh=shared_library('bvh', ['bvh/vec.h','bvh/aabb.h','bvh/bvh.h','bvh/prims.h','bvh/moller.h','bvh/raycast.h','bvh/serialization.h'] , version : '0.1.0', install : false, include_directories : my_inc)
libbvh_dep = declare_dependency(link_with : libbvh, include_directories : my_inc)
subdir('tests')
subdir('examples')
