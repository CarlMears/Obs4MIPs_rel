project('calc_tb_sounder', 'c', 'fortran', version: '0.1.0', default_options: ['buildtype=debug','warning_level=0','b_vscrt=static_from_buildtype'])

fortran_compiler = meson.get_compiler('fortran')

if fortran_compiler.get_id() == 'gcc'
  # Flags for gfortran (GCC)
  gfortran_flags = [
    '-ffree-form',          # Use free-form source code
    '-ffree-line-length-none', # Disable line length limits
    '-fimplicit-none',      # Require explicit declaration of all variables
    '-fbounds-check',        # Enable array bounds checking
    '-fno-range-check'     # Disable range checking (needed for NAN support)
  ]
  add_global_arguments(gfortran_flags, language : 'fortran')
  fortran_args = gfortran_flags
  package_ext = ''
else
    if get_option('buildtype') == 'debug'
        dialect = ['/nologo', '/free', '/fpp', '/debug:full']
        package_ext = '-d'
    else
        dialect = ['/nologo', '/free', '/fpp']
        package_ext = ''
    endif
    # Setup global arguments.
    bounds = '/CB'

    # make record length be in units of bytes, not 4-byte words.
    fortran_args = ['/assume:byterecl']
    add_global_arguments(dialect, language : 'fortran')
    add_global_arguments(bounds, language : 'fortran')
endif

#py = import('python').find_installation(required: true)
py_mod = import('python')
py = py_mod.find_installation('python',pure: false)
py_dep = py.dependency()
netcdf_dep = dependency('netcdf-fortran', required: false)
if not netcdf_dep.found()
  netcdf_dep = dependency('netcdff', required: true)
endif


fortran_src = [
  'src/AMSU_Constants.f90',
  'src/rtm.f90',
  'src/rtm_tables_AMSU.f90',
  'src/Calc_Tb_multiview_AMSU.f90',
]
f2py_wrapper = custom_target('f2py_wrapper',
  input : fortran_src,
  output : ['calc_tb_soundermodule.c', 'calc_tb_sounder-f2pywrappers2.f90'], # f2py often generates a second .f file
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'calc_tb_sounder', '--lower'] #,
  #'only:','compute_absorption_tables_amsu_dsb','compute_absorption_tables_amsu', ':']
) # I had to change output files to be the names of the generated files, because f2py doesn't always generate the same .f file name, 
  # and meson needs to know the exact output file names to use them in the extension module.

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true
).stdout().strip()

fortranobject_c = incdir_f2py / 'fortranobject.c'

# Add the manylinux-specific Fortran module path
extra_fortran_inc = '/usr/lib64/gfortran/modules'
# Verify it exists to avoid warnings on other systems
if not fortran_compiler.has_header('netcdf.mod', args: '-I' + extra_fortran_inc)
  # Fallback: if it's not there, it might be in the standard include
  extra_fortran_inc = '/usr/include'
endif

py.extension_module('calc_tb_sounder',
  [f2py_wrapper, fortran_src,fortranobject_c],
  include_directories: [incdir_numpy, incdir_f2py, extra_fortran_inc],
  dependencies : [py_dep, netcdf_dep],
  fortran_args : fortran_args,
  install : true,
  subdir : 'calc_tb_sounder',
)

# py.install_sources('data/abs_tables/amsu_05_abs_table_q_per_Pa.dat', subdir : 'calc_tb_sounder/data/abs_tables')
# py.install_sources('data/abs_tables/amsu_05_cld_abs_table.dat', subdir : 'calc_tb_sounder/data/abs_tables')
# py.install_sources('data/emiss_tables/amsu_05_emiss_table_sea_ice.dat', subdir : 'calc_tb_sounder/data/emiss_tables')
# py.install_sources('data/emiss_tables/amsu_05_emiss_table_W.dat', subdir : 'calc_tb_sounder/data/emiss_tables')

py.install_sources('data/abs_tables/amsu_05_abs_table_per_Pa_q.4.4.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_05_abs_table_per_Pa_q.4.5.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_07_abs_table_per_Pa_q.4.4.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_07_abs_table_per_Pa_q.4.5.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_09_abs_table_per_Pa_q.4.4.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_09_abs_table_per_Pa_q.4.5.nc', subdir : 'calc_tb_sounder/data/abs_tables')


# py.install_sources('data/abs_tables/amsu_07_abs_table_q_per_Pa.dat', subdir : 'calc_tb_sounder/data/abs_tables')
# py.install_sources('data/abs_tables/amsu_07_cld_abs_table.dat', subdir : 'calc_tb_sounder/data/abs_tables')
# py.install_sources('data/emiss_tables/amsu_07_emiss_table_sea_ice.dat', subdir : 'calc_tb_sounder/data/emiss_tables')
# py.install_sources('data/emiss_tables/amsu_07_emiss_table_W.dat', subdir : 'calc_tb_sounder/data/emiss_tables')

# py.install_sources('data/abs_tables/amsu_09_abs_table_q_per_Pa.dat', subdir : 'calc_tb_sounder/data/abs_tables')
# py.install_sources('data/abs_tables/amsu_09_cld_abs_table.dat', subdir : 'calc_tb_sounder/data/abs_tables')
# py.install_sources('data/emiss_tables/amsu_09_emiss_table_sea_ice.dat', subdir : 'calc_tb_sounder/data/emiss_tables')
# py.install_sources('data/emiss_tables/amsu_09_emiss_table_W.dat', subdir : 'calc_tb_sounder/data/emiss_tables')

py.install_sources('data/ocean_emiss_tables/ocean_emissivity_table_AMSU_channel_05.nc', subdir : 'calc_tb_sounder/data/ocean_emiss_tables')
py.install_sources('data/ocean_emiss_tables/ocean_emissivity_table_AMSU_channel_07.nc', subdir : 'calc_tb_sounder/data/ocean_emiss_tables')
py.install_sources('data/ocean_emiss_tables/ocean_emissivity_table_AMSU_channel_09.nc', subdir : 'calc_tb_sounder/data/ocean_emiss_tables')

py.install_sources('data/abs_tables/amsu_05_cld_abs_table.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_07_cld_abs_table.nc', subdir : 'calc_tb_sounder/data/abs_tables')
py.install_sources('data/abs_tables/amsu_09_cld_abs_table.nc', subdir : 'calc_tb_sounder/data/abs_tables')

py.install_sources('data/emiss_tables/amsu_05_emiss_table_sea_ice.nc', subdir : 'calc_tb_sounder/data/emiss_tables')
py.install_sources('data/emiss_tables/amsu_07_emiss_table_sea_ice.nc', subdir : 'calc_tb_sounder/data/emiss_tables')
py.install_sources('data/emiss_tables/amsu_09_emiss_table_sea_ice.nc', subdir : 'calc_tb_sounder/data/emiss_tables')