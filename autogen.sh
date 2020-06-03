#! /bin/sh

test -e README || ln -s README.md README
aclocal
autoheader 
autoconf
automake --add-missing
