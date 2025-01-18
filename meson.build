project('bvh', 'cpp',
        version : '1.0.0',
        default_options : ['warning_level=3', 'cpp_std=c++17'])
my_inc = include_directories(['.'])
bvh = executable('bvh', 'main.cpp', install : false, include_directories : my_inc)

test('test', bvh)
subdir('tests')
