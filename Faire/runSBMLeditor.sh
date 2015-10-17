#!/bin/sh

cd `dirname $0`

## libsbml configuration
## comment these lines or modify them if libsbml is installed somewhere else.

export OS=`uname`
export OS_ARCH=`uname -m`

# usual location
LIBSBML_HOME=./lib/libsbml-libxml2-4.0.0-1.i386
LD_LIBRARY_PATH=${LIBSBML_HOME}/lib

if [ "$OS_ARCH" == "x86_64" ] ; then
    LIBSBML_HOME=./lib/libsbml-libxml2-4.0.0-1.x86_64
    LD_LIBRARY_PATH=${LIBSBML_HOME}/lib64
fi

if [ "$OS" == "Darwin" ] ; then
    LIBSBML_HOME=./lib/libsbml-libxml2-4.0.0.macosx
    DYLD_LIBRARY_PATH=${LIBSBML_HOME}/lib
    export DYLD_LIBRARY_PATH
fi

export LD_LIBRARY_PATH

# custom location LIBSBML_HOME=/usr/lib/libsbml-4.0.0

LIBSBML_VERSION=4.0.0


SBML_CLASSPATH=./lib/libsbmlj-${LIBSBML_VERSION}.jar:./lib/libsbml-${LIBSBML_VERSION}-utils.jar

## end of libsbml configuration



JCOMPNEUR_CLASSPATH=./lib/compneur-2.1.jar
CLASSPATH=./lib/jsbml-1.0-a1-20120628-with-dependencies.jar:$SBML_CLASSPATH:./lib/SBWCore.jar:$JCOMPNEUR_CLASSPATH:lib/SBMLeditor.jar:./lib/log4j-1.2.15.jar:./lib/libsbmlj-5.5.0.jar

## add this parameter between java and -Xms512m if you want to have different properties files for each user
## -DSBMLeditor.properties.dir="%HOMEDRIVE%%HOMEPATH%" under windows
## -DSBMLeditor.properties.dir=$HOME/.SBML under Linux/Mac


java -DSBMLeditor.properties.dir=$HOME/.SBML -Xmx512m -classpath ${CLASSPATH} uk.ac.ebi.sbml.editor.Main
