bin_PROGRAMS = idealfluid
idealfluid_SOURCES = main.c 
if WITH_QUAD 
 idealfluid_SOURCES += mapping/discq.c mapping/halfplaneq.c
else
 idealfluid_SOURCES += mapping/disc.c mapping/halfplane.c
endif

