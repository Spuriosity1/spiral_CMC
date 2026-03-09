project('classical_MC', ['cpp', 'c'],
  version: '1.0.0',
  meson_version : '>= 1.5.0',
  default_options : [ 'warning_level=3', 'buildtype=release', 'cpp_std=c++2a' ]
)

if get_option('buildtype') == 'release'
  add_project_arguments('-Ofast',  language : 'cpp')
endif

json_dep = dependency('nlohmann_json', required: true)

latlib_dep = dependency('liblatindex', version : '>=1.1', required: true)

fftw_dep = dependency('fftw3', required: true)

hdf5_dep = dependency('hdf5', required: true, language: 'c')

argparse_dep = dependency('argparse', required: true)

cc=meson.get_compiler('c')

# Detect whether HDF5 was built with MPI support
if cc.compiles('''
    #include <hdf5.h>
    int main() { return 0; }
''',
    name: 'Check for parallel HDF5',
    dependencies: hdf5_dep
)

    message('Serial HDF5 detected')
else
 # suspect MPI is needed
    # Get MPI headers only (no linking)
    mpi_dep = dependency('mpi', required: true, include_type: 'system')

    # Combine deps so your target automatically gets mpi.h
    hdf5_dep = declare_dependency(
        dependencies: [hdf5_dep, mpi_dep]
    )
    message('parallel HDF5 detected')
endif   




inc=[
  include_directories(['include'])
  ]

g_deps = [
    latlib_dep,
    json_dep,
    fftw_dep,
    hdf5_dep,
    argparse_dep
    ]

subdir('src')

subdir('test')
