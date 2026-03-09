

# build the MC library
subdir('MC')

anneal_bin = executable('anneal', 
  files('anneal.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: inc,
  install: true
  )

spiral_optimise = executable('spiral_optimise',
  files('spiral_optimise.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: inc
  )

anneal_quench_bin = executable('anneal_and_quench', 
  files('anneal_and_quench.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: inc,
  install: true
  )

bzdump = executable('bzdump', 
  files('bzdump.cpp'),
  dependencies: g_deps,
  link_with: mc_lib,
  include_directories: inc,
  install: true
  )

