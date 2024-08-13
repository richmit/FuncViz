#!/bin/bash

# Quick little script to report on lines of code.

echo ""
echo "Non-C++"
cloc --force-lang=text,org --match-f='(tex|org|sh|txt|cmake)$' .

echo ""
echo "C++"
cloc --force-lang=text,org --match-f='(cpp|hpp)$' .
