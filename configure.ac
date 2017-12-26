AC_INIT([ideal-fluid], [1.0], [urrfinjuss@gmail.com], [ideal-fluid], [https://github.com/urrfinjuss/ideal-fluid])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AM_SILENT_RULES([yes])
AC_PROG_CC
AC_CHECK_LIB([fftw3],[fftw_malloc],[],[AC_MSG_ERROR([FFTW not installed])])
AC_CHECK_LIB([m],[pow],[],[AC_MSG_ERROR([Could not find pow])])
AC_CONFIG_HEADERS([src/headers/config.h])
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
AC_ARG_WITH([quad-precision], 
  [AS_HELP_STRING([--with-quad-precision], [Use __float128 as default requires quad-precision FFTW])],
  [], 
  [with_quad=no] 
)

# define USE_QUAD and compile __float128 sources
AS_IF([test "x$with_quad" != xno],
	[
	AC_CHECK_LIB([fftw3q],[fftwq_malloc],[],[AC_MSG_ERROR([quad precision FFTW not installed])])
	AM_CONDITIONAL(WITH_QUAD, true)
	AC_DEFINE([USE_QUAD],   [], [Use __float128 as default requires quad-precision FFTW library])
	], [
	AM_CONDITIONAL(WITH_QUAD, false)
	])

AC_OUTPUT
