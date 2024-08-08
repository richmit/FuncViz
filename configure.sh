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

  Use: configure.sh [cmake arguments]

    Common Arguments:
     * Target -- leave it off to get the default
       - -G 'MSYS Makefiles'           <-- Default on MSYS2
       - -G 'Visual Studio 17 2022'
       - -G 'Unix Makefiles'           <-- Default on Linux ('Linux' means 'Not MSYS2')
       - -G Ninja
     * Compiler -- leave it off to get the default
       - -DCMAKE_CXX_COMPILER=clang++
       - -DCMAKE_CXX_COMPILER=g++      <-- Default on MSYS2
       - -DCMAKE_CXX_COMPILER=g++-14   <-- Default on Linux if /usr/bin/g++-14 exists
       - -DCMAKE_CXX_COMPILER=g++      <-- Default on Linux if g++-14 wasn't found
     * Optional features -- leave them off to enable everything
       - -DO_DOXYGEN=[YES|NO]  -- For documentation
       - -DO_BTEST=[YES|NO]    -- Used for BOOT unit tests
EOF
exit
fi

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
if [ -e ../CMakeLists.txt ]; then
  if [ "$(basename $(pwd))" == "build" ]; then
    CMAKE_G=''
    if [[ "$CMD_LINE_ARGS" != *'-G'* ]]; then
      if [ -n "$MSYSTEM" ]; then
        CMAKE_G='MSYS Makefiles'
      else
        CMAKE_G='Unix Makefiles'
      fi
    fi
    if [ -n "$1" ]; then
      echo cmake -G "$CMAKE_G" "$@" ../
      cmake -G "$CMAKE_G" "$@" ../
    else
      CMAKE_C=''
      if [ -n "$MSYSTEM" ]; then
        CMAKE_C='-DCMAKE_CXX_COMPILER=g++.exe'
      else
        if [ -x '/usr/bin/g++-14' ]; then
          CMAKE_C='-DCMAKE_CXX_COMPILER=g++-14'
        else
          CMAKE_C='-DCMAKE_CXX_COMPILER=g++'
        fi
      fi
      echo cmake -G "$CMAKE_G" "$CMAKE_C" ..
      cmake -G "$CMAKE_G" "$CMAKE_C" ..
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
