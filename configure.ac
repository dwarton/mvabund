
AC_INIT([mvabund], 1.0-1)


AC_CHECK_LIB(gsl, gsl_matrix_set,
             [LIBS="${LIBS} -lgsl -lgslcblas"], 
             [AC_MSG_ERROR([Cannot find Gnu Scientific Library.])],
             [-lgslcblas])
AC_CHECK_HEADER([gsl/gsl_version.h], ,
             [AC_MSG_ERROR([Cannot find GSL headers.])])
AC_MSG_CHECKING([if GSL version >= 1.10])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit (minor < 10); /* GNU GSL 1.10 tested for mvabund */
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.10])
else
  AC_MSG_RESULT([yes])
fi

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1) 
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
    # Use gsl-config for header and linker arguments
    GSL_CFLAGS=`${GSL_CONFIG} --cflags`
    GSL_LIBS=`${GSL_CONFIG} --libs`
else 
    AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

## Use Rscript to query Rcpp for compiler and linker flags
# Rcpp compile flag providing header directory containing Rcpp.h
RCPP_CXXFLAGS=`"${R_HOME}/bin/Rscript" -e 'Rcpp:::CxxFlags()'`
## link flag providing libary as well as path to library, and optionally rpath
RCPP_LDFLAGS=`"${R_HOME}/bin/Rscript" -e 'Rcpp:::LdFlags()'`

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_SUBST(RCPP_CXXFLAGS)
AC_SUBST(RCPP_LDFLAGS)

AC_OUTPUT(src/Makevars)
