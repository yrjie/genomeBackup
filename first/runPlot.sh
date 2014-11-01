if [ $# -lt 1 ]
then
    echo "Usage: sh runPlot.sh cell"
    exit
fi
cell=$1
#paste output/$cell'FP.tab' output/$cell'M.tab' output/$cell'DF.tab' output/$cell'Fseq.tab' >temp.dat
#paste output/$cell'FP.tab' output/$cell'M.tab' output/$cell'DF.tab' output/$cell'Fseq.tab' output/$cell'Z.tab' >temp.dat
paste output/$cell'FP.tab' output/$cell'M.tab' output/$cell'DF.tab' output/$cell'Fseq.tab' output/$cell'Z.tab' output/$cell'NoSt.tab' >temp.dat
echo "call \"plotRoc.gp\" \"$cell\"" |gnuplot
