## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: fft.m4 10978 2013-07-11 15:28:46Z micael $
##

AC_DEFUN([ACX_FFT],
[
  acx_fft_ok=no

  if test "x${SINGLE_PRECISION}" != x; then
    fft_func="sfftw_plan_dft_1d"
    fft_libs="fftw3f"
  else
    fft_func="dfftw_plan_dft_1d"
    fft_libs="fftw3"
  fi

  dnl Check if the library was given in the command line
  if test $acx_fft_ok = no; then
    AC_ARG_WITH(fft-lib, [AS_HELP_STRING([--with-fft-lib=<lib>], [use FFT library <lib>])])
    case $with_fft_lib in
      yes | "") ;;
      -* | */* | *.a | *.so | *.so.* | *.o) LIBS_FFT="$with_fft_lib" ;;
      *) LIBS_FFT="-l$with_fft_lib" ;;
    esac
  fi

  acx_fft_save_LIBS="$LIBS"
  LIBS="$LIBS_FFT $LIBS $FLIBS"

  dnl First, check LIBS_FFT environment variable
  if test $acx_fft_ok = no; then
    if test "x$LIBS_FFT" != x; then
      AC_MSG_CHECKING([for $fft_func in $LIBS_FFT])
      AC_LINK_IFELSE([AC_LANG_CALL([], [$fft_func])], [acx_fft_ok=yes], [])
      if test $acx_fft_ok = no; then
        AC_MSG_RESULT([$acx_fft_ok])
      else
        AC_MSG_RESULT([$acx_fft_ok ($LIBS_FFT)])
      fi
    fi
  fi

  dnl FFTW linked to by default?
  if test $acx_fft_ok = no; then
    AC_CHECK_FUNC($fft_func, [acx_fft_ok=yes])
  fi

  dnl search libraries
  for fftl in $fft_libs; do
    if test $acx_fft_ok = no; then
      AC_CHECK_LIB($fftl, $fft_func,
        [acx_fft_ok=yes; LIBS_FFT="$LIBS_FFT -l$fftl"], [], [$FLIBS])
    fi
  done

  AC_MSG_CHECKING([for dfftw_init_threads in $LIBS_FFT])
  AC_LINK_IFELSE([AC_LANG_CALL([], [dfftw_init_threads])], [acx_fft_threads_ok=yes], [acx_fft_threads_ok=no])
  AC_MSG_RESULT([$acx_fft_threads_ok])
  if test x"$acx_fft_threads_ok" != xyes; then
    LIBS="$LIBS -lfftw3_threads"
    AC_MSG_CHECKING([for dfftw_init_threads in $LIBS_FFT -lfftw3_threads])
    AC_LINK_IFELSE([AC_LANG_CALL([], [dfftw_init_threads])],
      [acx_fft_threads_ok=yes; LIBS_FFT="$LIBS_FFT -lfftw3_threads"], [acx_fft_threads_ok=no])
    AC_MSG_RESULT([$acx_fft_threads_ok])
  fi
  if test x"$acx_fft_threads_ok" = xyes; then  
    AC_DEFINE(HAVE_FFTW3_THREADS, 1,[Define if the threaded version of FFTW3 is available.])
  fi

  AC_SUBST(LIBS_FFT)
  LIBS="$acx_fft_save_LIBS"

  # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
  if test x"$acx_fft_ok" != xyes; then
    if test $acx_fft_ok != disable; then
      AC_MSG_ERROR([Could not find required FFT library.])
    fi
    $2
  fi
])
