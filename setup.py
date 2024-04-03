import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand
from shutil import copyfile, copymode


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ['--', '-j30']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        # edit .bashrc file
        bashrcpath=os.path.join(os.path.expanduser("~"),'.bashrc')
        datapat=re.compile(r'export SCUBA_DATAPATH=(\S*)\n?')
        with open(bashrcpath,'r') as f:
            lines=f.readlines()
            lines=[datapat.sub('export SCUBA_DATAPATH={}\n'.format(ext.sourcedir+"/pySCUBA/data"),line) for line in lines]
        if not "SCUBA_DATAPATH" in os.environ.keys():
            lines.append('export SCUBA_DATAPATH={}\n'.format(ext.sourcedir+"/pySCUBA/data"))
        if not "PYSCUBA_WORKSPACE" in os.environ.keys():
            lines.append('export PYSCUBA_WORKSPACE={}\n'.format(os.path.join(os.path.expanduser("~"),'pySCUBA_workspace')))    
        with open(bashrcpath,'w') as f:
            f.writelines(lines)
        subprocess.check_call(['bash','-c','source '+bashrcpath])
modules=[
    CMakeExtension('pySCUBA'),
    ]

setup(
    name='pySCUBA',
    version='1.0',
    author='yhh',
    author_email='yuhh@mail.ustc.edu.cn',
    install_requires=[
        'numpy'
        ],
    description=' ',
    long_description='',
    packages=find_packages(''),
    # package_dir={'':'pySCUBA/pybind11/'},
    ext_modules=modules,
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)