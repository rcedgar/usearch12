#!/bin/bash -e

msg="$1"

if [ x"$msg" == x"" ] ; then
	msg=WIP
fi

git add .
git commit -m "$msg"
echo cpp= `ls *.cpp | wc -l` h= `ls *.h | wc -l`
