#!/bin/bash

grep -EA 10000 '(!|#|\*|=|_|-|%|;|C|/){70,}\.H\.E\.' $1 | grep -vE '(!|#|\*|=|_|-|%|;|C|/){70,}' | grep -vE '@(cond|endcond)' | perl -0 -p -e 's/\n{3,}/\n\n/g' | sed '1{/^$/d}'
