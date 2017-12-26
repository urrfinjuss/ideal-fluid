bin_PROGRAMS = idealfluid
idealfluid_SOURCES = main.c 
if WITH_QUAD 
 idealfluid_SOURCES += mapping/discq.c mapping/halfplaneq.c memory/memoryq.c timemarching/runge_kutta_4q.c array/array_funcq.c math/equationsq.c
else
 idealfluid_SOURCES += mapping/disc.c mapping/halfplane.c 
 idealfluid_SOURCES += memory/memory.c 
 idealfluid_SOURCES += timemarching/runge_kutta_4.c timemarching/stepping.c 
 idealfluid_SOURCES += array/array_func.c 
 idealfluid_SOURCES += math/equations.c
 idealfluid_SOURCES += io/input.c
 idealfluid_SOURCES += messages/misc.c
endif

