#!/bin/bash

# NOTE: This is a bash script that works within MSYS2 on WINDOWS!

# NOTE: Run from the build directory!

# NOTE: This is also a horrible way to do functional tests with floating point output.

if [ "$(basename $(pwd))" == "build" -a -e '../CMakeLists.txt' ]; then
  echo " "
  echo " "
  echo "Building Tests Now"
  echo " "
  for f in ../ftests/*.vtu; do
    t=$(echo $(basename $f) | sed 's/\.vtu$//')
    make "$t"
  done
  echo " "
  echo " "
  echo "Running Tests Now"
  echo " "
  for f in ../ftests/*.vtu; do
    t=$(echo $(basename $f) | sed 's/\.vtu$//')
    b="$t"
    if [ -x "$b"'.exe' ]; then
      b="$t"'.exe'
    fi
    if [ -x "$b" ]; then
      ./"$b" > /dev/null
      if ruby ../ftests/float_diff.rb "$t".vtu ../ftests/"$t".vtu >/dev/null; then
        echo "PASS: $t"
      else
        echo "FAIL: $t (invalid output)"
      fi
    else
      echo "FAIL: $t (no executable)"
    fi
  done
  echo " "
else
  echo "ERROR(run_tests.sh): Run this script from the build directory!"
fi
