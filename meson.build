project('FreshNBody', 'cpp', version : '0.1.0', default_options : ['cpp_std=c++17', 'warning_level=3'])

omp_dep = dependency('openmp')

source_files = files('src/particle.cpp', 'src/leapfrog.cpp', 'src/mercurius.cpp', 'src/whfast.cpp', 'src/ias15.cpp', 'src/simulation.cpp', 'src/accel.cpp', 'src/transform.cpp')

inc = include_directories('include')

lib_freshnbody = library('freshnbody', source_files, dependencies : omp_dep, include_directories : inc, install : true)

freshnbody_dep = declare_dependency(
  include_directories : inc,
  link_with : lib_freshnbody
)

pkg = import('pkgconfig')
pkg.generate(
  lib_freshnbody,
  name : 'freshnbody',
  description: 'simulation library',
  subdirs: 'freshnbody'
)

do_unit_tests = get_option('build-tests')

install_subdir('include', install_dir : join_paths(get_option('includedir'), 'freshnbody'), strip_directory : true)

if do_unit_tests
  gtest_main_dep = dependency('gtest_main', fallback : ['gtest', 'gtest_main_dep'])

  test_files = {
    'test1': 'test1.cpp',
    'leapfrog-basic': 'leapfrog-basic.cpp',
    'leapfrog-compensated': 'leapfrog-compensated.cpp',
    'jacobitransformtest': 'jacobitransformtest.cpp',
    'whfast': 'whfast.cpp'
  }

  foreach test_name, test_file : test_files
    test_file = 'unit-tests' / test_file
    test_exe = executable(test_name,
      sources : test_file,
      include_directories : inc,
      link_with : lib_freshnbody,
      dependencies : gtest_main_dep,
      install : false,
      override_options : ['cpp_std=c++20'] 
    )

    test(test_name, test_exe)
    endforeach
endif
