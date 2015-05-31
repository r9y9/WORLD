APPNAME = 'world'
VERSION = '0.2.0'

from waflib import Options
import sys
import os
import re
import waflib

subdirs = [
    'src',
    'example'
]

top = '.'
out = 'build'


def options(opt):
    opt.load('compiler_cxx')


def configure(conf):
    conf.check_waf_version(mini='1.7.11')
    conf.load('compiler_cxx')
    conf.define('WORLD_VERSION', VERSION)
    conf.env['VERSION'] = VERSION
    conf.env['APPNAME'] = APPNAME

    if conf.env.COMPILER_CXX != 'msvc':
        conf.env.append_unique('CXXFLAGS', ['-O2', '-Wall'])
    conf.env.HPREFIX = conf.env.PREFIX + '/include/world'
    conf.recurse(subdirs)

    print """
world has been configured as follows:

[Build information]
Package:                 %s
build (compile on):      %s
host endian:             %s
Compiler:                %s
Compiler version:        %s
CXXFLAGS:                %s
""" % (
        APPNAME + '-' + VERSION,
        conf.env.DEST_CPU + '-' + conf.env.DEST_OS,
        sys.byteorder,
        conf.env.COMPILER_CXX,
        '.'.join(conf.env.CC_VERSION),
        ' '.join(conf.env.CXXFLAGS)
    )

    conf.write_config_header('src/world-config.h')


def build(bld):
    bld.recurse(subdirs)

    libs = []
    for tasks in bld.get_build_iterator():
        if tasks == []:
            break
        for task in tasks:
            if isinstance(task.generator, waflib.TaskGen.task_gen) and 'cxxshlib' in task.generator.features:
                libs.append(task.generator.target)
    ls = ''
    for l in set(libs):
        ls = ls + ' -l' + l

    bld(source='world.pc.in',
        prefix=bld.env['PREFIX'],
        exec_prefix='${prefix}',
        libdir=bld.env['LIBDIR'],
        libs=ls,
        includedir='${prefix}/include',
        PACKAGE=APPNAME,
        VERSION=VERSION)
