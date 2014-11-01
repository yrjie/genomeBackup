if [ $# -lt 1 ]
then
    echo 'Usage: -p iffFold'
    exit
fi
JAVAOPTS="-Xms100M -Xmx16000M"
#JAVAOPTS="-Xmx8000M"

#outdir="test"
#outdir="K562"
java $JAVAOPS -cp Fpeak.jar:lib/commons-cli-1.1.jar fpeak.TestBG $*
