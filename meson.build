project('smallpt', ['cpp'],
  version : '0.1.0',
  default_options : [
    'warning_level=3',
    'cpp_std=c++17',
  ]
)

executable('smallpt', files(['src/main.cpp']),
  cpp_args : ['-fopenmp'],
  link_args : ['-fopenmp'])
