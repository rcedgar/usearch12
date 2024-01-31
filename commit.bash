#!/bin/bash -e

git add .
git commit -m "WIP"
echo cpp= `ls *.cpp | wc -l` h= `ls *.h | wc -l`
