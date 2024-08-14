#!/bin/bash

grep -EB 10000 '(!|#|\*|=|_|-|%|;|C|/){70,}\.H\.E\.' $1 | grep -EA 10000 '@filedetails' | head -n -2 | grep -v '@filedetails' | \
  sed -e 's/ @f\$/ \\(/g; s/@f\$\([ .,!:]\)/\\)\1/g;     # Rewrite inline math'           \
      -e 's/@f\[/\\[/g; s/@f]/\\]/g;                     # Rewrite displaied maty'        \
      -e 's/^ *\\verbatim/#+BEGIN_EXAMPLE/;              # verbatim -> BEGIN_EXAMPLE'     \
      -e 's/^ *\\endverbatim/#+END_EXAMPLE/;             # endverbatim -> END_EXAMPLE'    \
      -e 's/ \([0-9a-zA-Z_-]*\).cpp/ [[#\1][\1.cpp]]/g   # links from source file names'
