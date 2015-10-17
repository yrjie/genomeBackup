#!/bin/sh

cd `dirname $0`

JCOMPNEUR_CLASSPATH=./lib/compneur-2.1.jar

SBW_CLASSPATH=./lib/SBWCore.jar

COMMAND=java

XML_CLASSPATH=./lib/jsbml-1.0-a1-20120628-with-dependencies.jar:./lib/libsbmlj-5.5.0.jar

## libsbml configuration
## comment these lines or modify them if libsbml is installed somewhere else.

# usual location
LIBSBML_HOME=/usr/local

LD_LIBRARY_PATH=${LIBSBML_HOME}/lib:$LD_LIBRARY_PATH

$COMMAND -Xmx256m -classpath $XML_CLASSPATH:$SBW_CLASSPATH:$JCOMPNEUR_CLASSPATH:lib/SBMLeditor.jar uk.ac.ebi.sbml.editor.Main -sbwregister

