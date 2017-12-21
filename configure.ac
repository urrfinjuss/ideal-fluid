AC_INIT([ideal-fluid], [1.0], [urrfinjuss@gmail.com], [ideal-fluid], [https://github.com/urrfinjuss/ideal-fluid])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AM_SILENT_RULES([yes])
AC_PROG_CC

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
src/Makefile
Makefile
])
# enable debugging
AC_ARG_ENABLE([debug], 
  [AS_HELP_STRING([--enable-debug], [enable debug messages])],
[ AC_DEFINE([DEBUG], [1], [Enable extra debugging information]) ], 
[ AC_DEFINE([DEBUG], [0], [Disable extra debugging information]) ] 
)

# enable quad-precision
AC_ARG_ENABLE([quad-precision], 
  [AS_HELP_STRING([--enable-quad-precision], [Use __float128 as default requires quad-precision FFTW])],
  [ AC_DEFINE([USE_QUAD],   [1], [Use __float128 as default requires quad-precision FFTW library]) ], 
  [ AC_DEFINE([USE_DOUBLE], [1], [Use double as default requires standard FFTW library]) ] 
)
AC_OUTPUT
