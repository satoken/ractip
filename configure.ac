#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.15])
AC_INIT([ractip], [0.0.1], [satoken@k.u-tokyo.ac.jp])
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/ractip.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.


# gurobi optimizer
AC_ARG_WITH([gurobi], AS_HELP_STRING([--with-gurobi], [support Gurobi optimizer]),
  if test "$with_gurobi" != "yes" && test "$with_gurobi" != "no"; then
     CPPFLAGS="-I${with_gurobi}/include $CPPFLAGS"
     LDFLAGS="-L${with_gurobi}/lib $LDFLAGS"
  fi
)
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_HEADERS([gurobi_c++.h])
AC_CHECK_LIB([gurobi20], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi20 $LIBS")
#AC_CHECK_LIB([gurobi_c++], [GRBEnv])

# Vienna RNA package
AC_ARG_WITH([vienna-rna], AS_HELP_STRING([--with-vienna-rna], [support Vienna RNA package]),
  if test "$with_vienna_rna" != "yes" && test "$with_vienna_rna" != "no" ; then
     CPPFLAGS="-I${with_vienna_rna}/include $CPPFLAGS"
     LDFLAGS="-L${with_vienna_rna}/lib $LDFLAGS"
  fi
)
AC_CHECK_LIB(RNA, init_pf_fold)


AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
