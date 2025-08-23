project('seismic-rna', 'c', version: '0.24.3', default_options: ['c_std=c99'])
python = import('python').find_installation(pure: false)

# Install Python package files.
py_install_dir = python.get_install_dir()
install_subdir('src/seismicrna',
               install_dir: py_install_dir)

# Build and install C extension modules.
add_project_arguments(['-DNDEBUG', '-O2'],
                      language: 'c')
python.extension_module('relate',
                        sources: ['src/seismicrna/relate/cx/relate.c'],
                        subdir: 'seismicrna/relate/cx',
                        install: true)
