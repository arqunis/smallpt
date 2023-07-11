project('smallpt',
  version : '0.1.0',
  default_options : [
    'warning_level=3',
    'cpp_std=c++17',
  ]
)

if 'cpp' in get_option('languages')
  add_languages('cpp', required : true)

  smallpt_inc = include_directories(['cpp/'])
  smallpt_src = files([
    'cpp/smallpt/constants.hpp',
    'cpp/smallpt/primitives.hpp',
    'cpp/smallpt/primitives.cpp',
    'cpp/smallpt/vec.hpp',
    'cpp/smallpt/vec.cpp',
    'cpp/main.cpp',
  ])

  executable('smallpt', smallpt_src,
    cpp_args : ['-fopenmp'],
    link_args : ['-fopenmp'],
    include_directories : [smallpt_inc])
endif

if 'd' in get_option('languages')
  add_languages('d', required : true)

  smallptd_inc = include_directories(['d/'])
  smallptd_src = files([
    'd/smallpt/primitives.d',
    'd/smallpt/rng.d',
    'd/smallpt/vec.d',
    'd/main.d',
  ])
  executable('smallptd', smallptd_src,
    include_directories : [smallptd_inc])
endif
