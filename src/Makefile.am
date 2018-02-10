AM_CFLAGS = -Wno-format
bin_PROGRAMS = ffluid_disc ffluid_hplane


ffluid_disc_SOURCES = ffluid_disc.c
if WITH_QUAD 
 ffluid_disc_SOURCES += mapping/discq.c memory/memoryq.c math/equationsq.c
else
 ffluid_disc_SOURCES += array/array_func.c
 ffluid_disc_SOURCES += io/io.c io/input_parameters.c io/input_data.c io/output_disc.c
 ffluid_disc_SOURCES += mapping/mapping.c mapping/disc.c 
 ffluid_disc_SOURCES += math/equations.c math/disc_surface.c
 ffluid_disc_SOURCES += memory/memory.c
 ffluid_disc_SOURCES += messages/misc.c
 ffluid_disc_SOURCES += timemarching/stepping.c timemarching/runge_kutta_4.c timemarching/runge_kutta_6.c
endif

ffluid_hplane_SOURCES = ffluid_hplane.c
if WITH_QUAD 
 ffluid_hplane_SOURCES += mapping/halfplaneq.c memory/memoryq.c math/equationsq.c
else
 ffluid_hplane_SOURCES += array/array_func.c
 ffluid_hplane_SOURCES += io/io.c io/input_parameters.c io/input_data.c io/output_halfplane.c
 ffluid_hplane_SOURCES += mapping/mapping.c mapping/halfplane.c 
 ffluid_hplane_SOURCES += math/equations.c math/halfplane_surface.c 
 ffluid_hplane_SOURCES += memory/memory.c 
 ffluid_hplane_SOURCES += messages/misc.c
 ffluid_hplane_SOURCES += timemarching/stepping.c timemarching/runge_kutta_4.c timemarching/runge_kutta_6.c
endif
