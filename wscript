#!/usr/bin/env python
# encoding: utf-8

from waflib.Tools import waf_unit_test, python
from waflib.Task import Task
from waflib import Logs
from waflib import Utils
import os
import glob

top = '.'
out = 'build'

cppfiles = [
    #'EGFRDSimulatorWrapper.cpp',
    'epdp/Logger.cpp',
    'epdp/ConsoleAppender.cpp',
    'epdp/utils.cpp',
    'epdp/findRoot.cpp',
    'epdp/funcSum.cpp',
    'epdp/GreensFunction1DAbsAbs.cpp',
    'epdp/GreensFunction1DRadAbs.cpp',
    'epdp/PairGreensFunction.cpp',
    'epdp/GreensFunction3D.cpp',
    'epdp/GreensFunction3DAbs.cpp',
    'epdp/GreensFunction3DAbsSym.cpp',
    'epdp/GreensFunction3DRadAbsBase.cpp',
    'epdp/GreensFunction3DRadAbs.cpp',
    'epdp/GreensFunction3DRadInf.cpp',
    'epdp/GreensFunction3DSym.cpp',
    'epdp/SphericalBesselGenerator.cpp',
    'epdp/BasicNetworkRulesImpl.cpp',
    'epdp/Model.cpp',
    'epdp/NetworkRules.cpp',
    'epdp/ParticleModel.cpp',
    'epdp/SpeciesType.cpp',
    'epdp/StructureType.cpp',
    ]

def options(opt):
    opt.load('compiler_cxx waf_unit_test boost')

def configure(conf):
    conf.load('compiler_cxx waf_unit_test python boost')
    conf.check_cxx(lib = 'gsl')
    conf.check_cxx(lib = 'gslcblas')
    conf.check_cxx(lib = 'm')
    conf.check_python_module('scipy')

    conf.check(
        header_name = 'unordered_map.hpp',
        define_name = "HAVE_UNORDERED_MAP_HPP", mandatory = False)
    conf.check(
        header_name = 'boost/unordered_map.hpp',
        define_name = "HAVE_BOOST_UNORDERED_MAP_HPP", mandatory = False)
    conf.check(
        header_name = 'boost/functional/hash.hpp',
        define_name = "HAVE_BOOST_FUNCTIONAL_HASH_HPP", mandatory = False)

    conf.check_cc(fragment = '''
#include <tr1/unordered_map>
std::tr1::unordered_map<int, int> a, b( a); ''',
               define_name = "HAVE_TR1_UNORDERED_MAP",
               mandatory = False)

    conf.check_cc(fragment = '''
#include <tr1/functional>
int main() {std::tr1::hash<int>(); return 0; }''',
               define_name = "HAVE_TR1_FUNCTIONAL",
               mandatory = False)

    conf.check_cxx(fragment='''
#include <math.h>
int main() { double a = INFINITY; return (int)a * 0; }
''',
        define_name='HAVE_DECL_INFINITY',
        mandatory=False)
    conf.check_cxx(fragment='''
#include <math.h>
int main() { isfinite(0.); return 0; }
''',
        define_name='HAVE_ISFINITE',
        mandatory=False)
    conf.check_cxx(fragment='''
#include <math.h>
int main() { double a, b; sincos(0., &a, &b); return 0; }
''',
        define_name='HAVE_SINCOS',
        mandatory=False)

    # conf.check_cxx(lib = 'hdf5')
    # conf.check_cxx(lib = 'hdf5_cpp')

    conf.check_boost(lib = 'regex')
    #conf.check_cxx(lib = 'ecell4-core')

    conf.write_config_header('epdp/config.h', guard='__ECELL4_W_CONFIG_H_WAF')
    #conf.recurse(subdirs)


def pre(bld):
    sjy_table = bld.path.find_or_declare('epdp/SphericalBesselTable.hpp')
    if not os.path.exists(sjy_table.abspath()):
        bld.exec_command(' '.join([
            bld.env['PYTHON'][0],
            bld.path.find_resource('epdp/make_sjy_table.py').abspath(),
            sjy_table.abspath(),
            ]))
        Utils.check_dir(bld.env['PREFIX'] + '/include/ecell4/egfrd/epdp')
        bld.exec_command(' '.join([
            'cp', sjy_table.abspath(), '${PREFIX}/include/ecell4/egfrd/epdp/.']))
    cjy_table = bld.path.find_or_declare('epdp/CylindricalBesselTable.hpp')
    if not os.path.exists(cjy_table.abspath()):
        bld.exec_command(' '.join([
            bld.env['PYTHON'][0],
            bld.path.find_resource('epdp/make_cjy_table.py').abspath(),
            cjy_table.abspath(),
            ]))
        Utils.check_dir(bld.env['PREFIX'] + '/include/ecell4/egfrd/epdp')
        bld.exec_command(' '.join([
            'cp', cjy_table.abspath(), '${PREFIX}/include/ecell4/egfrd/epdp/.']))

def build(bld):
    bld.add_pre_fun(pre)

    #bld.install_files(
    #    '${PREFIX}/include/ecell4/egfrd', hppfiles)
    bld.install_files(
        '${PREFIX}/include/egfrd_core', ['epdp/config.h'])

    bld.install_files(
        '${PREFIX}/include/egfrd_core', glob.glob('epdp/*.hpp'))
    bld.install_files(
        '${PREFIX}/include/egfrd_core', ['epdp/compat.h'])
    bld.install_files(
        '${PREFIX}/include/egfrd_core/utils', glob.glob('epdp/utils/*.hpp'))

    bld.shlib(
        source = cppfiles,
        includes = ['.', 'epdp'],
        defines = ['HAVE_CONFIG_H', 'HAVE_INLINE'],
        #lib = ['ecell4-core', 'gsl', 'gslcblas', 'm', 'hdf5', 'hdf5_cpp'],
        lib = ['gsl', 'gslcblas', 'm', 'hdf5', 'hdf5_cpp'],
        use = ['BOOST'],
        target = 'egfrd_core')

    #bld.recurse(subdirs)

    #bld.add_post_fun(summary)
    #bld.options.all_tests = True
