AM_CFLAGS = -Wno-format
bin_PROGRAMS = idealfluid
idealfluid_SOURCES = main.c 
if WITH_QUAD 
 idealfluid_SOURCES += mapping/discq.c mapping/halfplaneq.c memory/memoryq.c timemarching/runge_kutta_4q.c array/array_funcq.c math/equationsq.c
else
 idealfluid_SOURCES += mapping/mapping.c mapping/disc.c mapping/halfplane.c 
 idealfluid_SOURCES += memory/memory.c 
 idealfluid_SOURCES += timemarching/stepping timemarching/runge_kutta_4.c timemarching/stepping.c 
 idealfluid_SOURCES += array/array_func.c 
 idealfluid_SOURCES += math/equations.c math/halfplane_surface.c math/disc_surface.c
 idealfluid_SOURCES += io/io.c io/input_parameters.c io/input_data.c io/output.c
 idealfluid_SOURCES += messages/misc.c
endif

