project('bvh', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++17','optimization=3'])
# Detect the host operating system
simd_args = ['-flto','-O3']
link_args=[ '-flto']
if host_machine.system() == 'darwin'
    #simd_args = ['-march=armv8-a+simd']
     simd_args += ['-mcpu=apple-m1','-funroll-loops']#,#'-flto', ]

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
bvh = executable('bvh', 'main.cpp', install : false, include_directories : my_inc)

test('test', bvh)
subdir('tests')
