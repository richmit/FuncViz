#!/bin/bash

# Find files too big for github that are not listed in .gitignore

if [ -e paraview -a -e docs/func-viz ]; then
  find paraview/     -size +50M -a -not -iname '*lossless.webm' -a -not -iname '*crf01.mp4' -exec sh -c "grep -q '^/{}$' .gitignore || echo '{}'" \;
  find docs/func-viz/ -size +50M
else
  echo "ERROR: Run from base of FuncViz repository"
fi
