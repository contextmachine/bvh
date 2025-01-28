import dataclasses
import subprocess
from dataclasses import field
import platform
import sys
import os
import warnings
from pathlib import Path
import re
from typing import Any

import numpy as np



def parse_meson_options(file_path):
    options = {}
    # Define regex to capture option names and values
    option_regex = re.compile(
        r"option\(\s*'(?P<name>\w+)'\s*,\s*type:\s*'(?P<type>\w+)'\s*,\s*value:\s*(?P<value>[^,]+),?"
    )
    with open(file_path, 'r') as file:
        for line in file:
            match = option_regex.search(line)
            if match:
                name = match.group('name')
                value = match.group('value').strip()
                # Convert value based on type
                if value in ["true", "false"]:
                    value = value == "true"
                elif value.isdigit():
                    value = int(value)
                elif value.startswith("'") and value.endswith("'"):
                    value = value.strip("'")
                options[name] = value
    return options

sys.path.append(os.getcwd())

if sys.platform.startswith("win"):
    import pyMSVC

    environment = pyMSVC.setup_environment()
    #print(environment)
opt=Path(__file__).parent/"meson_options.txt"
if not opt.exists():
    #print(opt ,"is missing")
    meson_options={}
else:
    meson_options=parse_meson_options(opt)

    #print(meson_options)

import setuptools
import numpy

# rest of setup code here
from setuptools import Extension, Distribution
from setuptools.command.build_ext import build_ext

from Cython.Build import cythonize
from pathlib import Path
compile_args = ["-O3","-std=c++17","-DNPY_NO_DEPRECATED_API=NPY_2_0_API_VERSION"]
link_args = [ ]
include_dirs = [numpy.get_include(),os.getcwd()]
define_macros = []



#print(sys.platform)
if sys.platform == "darwin" and platform.machine()=="arm64":
    #print("Darwin")

    subprocess.run("brew")
    os.environ['ARCHFLAGS'] = "-arch arm64"
    os.environ['CFLAGS'] = "-arch arm64"
    os.environ['CXXFLAGS'] = "-arch arm64 -stdlib=libc++"
    os.environ['LDFLAGS'] = "-arch arm64 -stdlib=libc++"
    os.environ['CC'] = 'clang'
    os.environ['CXX'] = 'clang++'
    compile_args += [ "-mcpu=apple-m1", "-funroll-loops","-flto"]

    OPENMP_PATH = '/opt/homebrew/opt/llvm'
    LLVM_PATH = '/opt/homebrew/opt/libomp'
    if not Path(LLVM_PATH).exists():
        warnings.warn('LLVM is not found. Install it using homebrew:\nbrew install llvm\n')
        os.environ['HOMEBREW_NO_AUTO_UPDATE'] = '1'
        subprocess.run('brew install llvm'.split(' '))
    if not Path(OPENMP_PATH).exists():
        warnings.warn('OpenMP is not found. Install it using homebrew:\nbrew install libomp\n')
        os.environ['HOMEBREW_NO_AUTO_UPDATE'] = '1'
        subprocess.run('brew install libomp'.split(' '))

    else:
        compile_args += [
            '-fopenmp'
        ]
        link_args += [f'-L{LLVM_PATH}/bin',f'-L{LLVM_PATH}/lib',f'-L{OPENMP_PATH}/lib',f'-L{OPENMP_PATH}/lib', '-fopenmp', '-lomp']

        include_dirs += [f'-L{LLVM_PATH}/include',f'{OPENMP_PATH}/include']

    link_args += compile_args



elif sys.platform == "win32":

    compile_args[1] = "/std:c++20"
    compile_args[0] = "/Ox"

    compile_args+=['/openmp','/fp:fast','/nologo', '/EHsc', '/MD' ,'/favor:Intel64' ]
    link_args+=['/openmp']

    link_args += compile_args

elif platform.machine()=="x86_64":

    sys.path.append('/usr/include/x86_64-linux-gnu')
    sys.path.append('/usr/lib/x86_64-linux-gnu')
    compile_args += ['-lm']
    include_dirs+=['/usr/include/x86_64-linux-gnu']
    link_args += ['-L/usr/lib/x86_64-linux-gnu']
    link_args += ['-fopenmp','-lgomp']
    compile_args += [
        '-msse', '-msse2', '-msse3', '-mssse3',
        '-msse4.1', '-msse4.2', '-mavx', '-mavx2'
    ]

elif platform.machine() in ["arm64","aarch64"]:

    compile_args+=[ '-lm','-mcpu=native','-march=armv8-a+simd']
   
    compiler_args = ['-flto']
    link_args += ['-fopenmp','-lgomp']

    link_args+=    compile_args

 
else:
    raise Exception("unknown platform")
# see pyproject.toml for other metadata


print(compile_args,link_args,include_dirs,define_macros)

logo = f"""
------------------------------------------------------------------------------------------------------------------------
    contextmachine bvh
------------------------------------------------------------------------------------------------------------------------
{platform.uname()}                                       
compile_args: {compile_args}
"""

extensions = [
Extension("bvh._bvh",   ["bvh/_bvh.pyx"],
          extra_compile_args=compile_args,

        extra_link_args=link_args,
          define_macros=define_macros,
          include_dirs=include_dirs,
          language="c++"),
          Extension("bvh.serialization",   ["bvh/serialization.pyx"],
          extra_compile_args=compile_args,

        extra_link_args=link_args,
          define_macros=define_macros,
          include_dirs=include_dirs,
          language="c++"),
]



compiler_directives = dict(
    boundscheck=False,
    wraparound=False,
    cdivision=True,
    nonecheck=False,
    overflowcheck=False,
    initializedcheck=False,
    embedsignature=False,
    language_level="3str",

)

print(__name__)
if __name__ == "__main__":

        print(logo)


        try:
            ext_modules = cythonize(
                extensions,
                include_path=include_dirs,
                compiler_directives=compiler_directives,
                force=True,

            )
            dist = Distribution({"ext_modules": ext_modules})
            cmd = build_ext(dist)
            cmd.ensure_finalized()
            cmd.run()
        except Exception as e:

                homebrew = os.getenv("HOMEBREW_PREFIX")
                if homebrew is None:
                    raise e
                if not (Path(homebrew)/'opt/llvm').exists() or not (Path(homebrew)/'opt/libomp').exists() :
                    os.environ['HOMEBREW_NO_AUTO_UPDATE']='1'
                    subprocess.run('brew install llvm'.split(' '))
                    subprocess.run('brew install libomp'.split(' '))
                    include_dirs.append(Path(homebrew)/'opt/llvm/include')
                    link_args.extend( [Path(homebrew)/'opt/llvm/bin','-Wl','-rpath',Path(homebrew)/'opt/llvm/lib'])

                    ext_modules = cythonize(
                                extensions,
                                include_path=include_dirs,
                                compiler_directives=compiler_directives,
                                force=True,

                            )
                    dist = Distribution({"ext_modules": ext_modules})
                    cmd = build_ext(dist)
                    cmd.ensure_finalized()
                    cmd.run()
                else:
                    raise e

        import os, shutil

        for output in cmd.get_outputs():
            print(output)
            relative_extension = os.path.relpath(output, cmd.build_lib)
            print(relative_extension)
            shutil.copyfile(output, relative_extension)

