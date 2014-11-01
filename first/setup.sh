#!/bin/bash
CHIAPETPATH=`pwd`

echo 'Compiling C codes...'
cd $CHIAPETPATH/src/c
make

echo 'Compiling 3rd-party codes...'
cd $CHIAPETPATH/3rd/batman/src
make

echo 'Compiling Java codes...'
cd $CHIAPETPATH/src/java
sh build.sh

