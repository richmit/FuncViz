#!/usr/bin/env -S bash
# -*- Mode:Shell-script; Coding:us-ascii-unix; fill-column:158 -*-
#########################################################################################################################################################.H.S.##
##
# @file      configure.sh
# @author    Mitch Richling http://www.mitchr.me/
# @date      2024-07-31
# @brief     Just a little helper for people accustomed to GNU autotools.@EOL
# @std       bash
# @copyright 
#  @parblock
#  Copyright (c) 2024, Mitchell Jay Richling <http://www.mitchr.me/> All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
#  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  @endparblock
#########################################################################################################################################################.H.E.##

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
if [[ "${@}" == *'-h'* ]]; then
  cat <<EOF

  Run this script from the 'build' directory.

  If you don't have a 'build' directory yet, then create one!

  Use: configure.sh [configure options] [cmake arguments]

    Configure Options
     * -C Clean the build directory before running cmake

    Common Cmake Arguments:
     * Target -- leave it off to get the default
       - -G 'MSYS Makefiles'           <-- Default on MSYS2
       - -G 'Visual Studio 17 2022'
       - -G 'Unix Makefiles'           <-- Default on Linux ('Linux' means 'Not MSYS2')
       - -G Ninja
     * Compiler -- leave it off to get the default
       - -DCMAKE_CXX_COMPILER=clang++
       - -DCMAKE_CXX_COMPILER=g++      <-- Default for 'MSYS Makefiles'
       - -DCMAKE_CXX_COMPILER=g++-##   <-- Default for 'Unix Makefiles' if /usr/bin/g++-[0-9][0-9] exists
                                           in which case the highest numbered version is selected.
                                           This code base needs at least GCC-14.
       - -DCMAKE_CXX_COMPILER=g++      <-- Default for 'Unix Makefiles' if /usr/bin/g++-[0-9][0-9] missing
       -                               <-- Default for 'Visual Studio 17 2022'
     * Optional features -- leave them off to enable everything
       - -DO_DOXYGEN=[YES|NO]  -- Doxygen (to build documentation)
       - -DO_BTEST=[YES|NO]    -- BOOT unit tests (to run unit tests)
       - -DO_MRASTER=[YES|NO]  -- MRaster (used for some examples)
EOF
exit
fi

if [ "$1" == "-C" ]; then
  shift
  CLEAN_DIR='YES'
else
  CLEAN_DIR='NO'
fi


#---------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ -e ../CMakeLists.txt ]; then
  if [ "$(basename $(pwd))" == "build" ]; then
    #
    # If we need to clean, then CLEAN!!
    if [ "$CLEAN_DIR" == 'YES' ]; then
      ls
      rm -rI *
    fi
    #
    # Figure out target
    CMAKE_TARGET=''
    CMAKE_TARGET_SRC='AUTO'
    for arg in "$@"; do
      if [ -z "$CMAKE_TARGET" -a "$CMAKE_TARGET_SRC" == 'COMMAND_LINE' ]; then
        CMAKE_TARGET="$arg"
      fi
      if [ "$arg" == "-G" ]; then
        CMAKE_TARGET_SRC='COMMAND_LINE'
      fi
    done
    if [ "$CMAKE_TARGET" == 'LOOKING' ]; then
      echo "ERROR: Found -G but no following target argument"
    fi
    # No -G, figure out default!
    if [ -z "$CMAKE_TARGET" ]; then
      if [ -n "$MSYSTEM" ]; then
        CMAKE_TARGET='MSYS Makefiles'
      else
        CMAKE_TARGET='Unix Makefiles'
      fi
    fi
    #
    # Figure out the -DCMAKE_CXX_COMPILER argument
    CMAKE_CARG=''
    if [[ "$@" != *'-DCMAKE_CXX_COMPILER'* ]]; then
      if [ "$CMAKE_TARGET" == 'MSYS Makefiles' ]; then
        CMAKE_CARG='-DCMAKE_CXX_COMPILER=g++.exe'
      fi
      if [ "$CMAKE_TARGET" == 'Unix Makefiles' ]; then
        HIGCC=$(ls /usr/bin/gcc-[0-9][0-9] 2>/dev/null | sort | tail -n 1)
        if [ -x "$HIGCC" ]; then
          CMAKE_CARG='-DCMAKE_CXX_COMPILER=g++-14'
        else
          CMAKE_CARG='-DCMAKE_CXX_COMPILER=g++'
        fi
      fi
    fi
    if [ "$DEBUG" == 'DEBUG' ]; then
      echo $CMAKE_TARGET
      echo $CMAKE_TARGET_SRC
      echo $CMAKE_CARG
    fi
    #
    # Run cmake
    if [ "$CMAKE_TARGET_SRC" == "COMMAND_LINE" ]; then
      if [ -z "$CMAKE_CARG" ]; then
        echo AA cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' "$@" ../
        cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' "$@" ../
      else
        echo BB cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' "$CMAKE_CARG" "$@" ../
        cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' "$CMAKE_CARG" "$@" ../
      fi
    else
      if [ -z "$CMAKE_CARG" ]; then
        echo CC cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' -G "$CMAKE_TARGET" "$@" ../
        cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' -G "$CMAKE_TARGET" "$@" ../
      else
        echo DD cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' -G "$CMAKE_TARGET" "$CMAKE_CARG" "$@" ../
        cmake '-DCMAKE_EXPORT_COMPILE_COMMANDS=1' -G "$CMAKE_TARGET" "$CMAKE_CARG" "$@" ../
      fi
    fi
  else
    echo "ERROR: Must run from build directory"
  fi
else
  if [ -e ./CMakeLists.txt ]; then
    echo "ERROR: It looks like you are running from the base of the repo"
    if [ -d ./build ]; then
      echo "ERROR: cd into the build directory to run this script"
    else
      echo "ERROR: Create a directory called 'bulid'."
      echo "ERROR: cd into the build directory to run this script"
    fi
  else
    echo "ERROR: Missing ../CMakeLists.txt"
  fi
fi
