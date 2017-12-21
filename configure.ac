AC_INIT([ideal-fluid], [1.0], [urrfinjuss@gmail.com])
AC_ARG_ENABLE(
  [quad-precision], 
  [AS_HELP_STRING([--enable-quad-precision], [enable the use of __float128 data type instead of double data type. Requires FFTW with quad-precision])],
  [AC_CHECK_LIB([fftw3q],[fftwq_malloc], [] , AC_MSG_ERROR([ Dependency not satisfied: FFTW with quad-precision is required ])) ],
  [AC_CHECK_LIB([fftw3], [fftw_malloc], [], AC_MSG_ERROR([Dependency not satisfied: FFTW library not found])) ]  
)
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AM_SILENT_RULES([yes])
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
src/Makefile
Makefile
])
AC_OUTPUT
