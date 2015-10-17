args=commandArgs()
if (length(args)<7){
print("Usage: R --no-save --slave --args infile outputdir prefix <runZinba.R")
q()
}

infile=args[5]
out=args[6]
prefix=args[7]

library(zinba)

# The alignability files are /home/ruijie/Faire/zinba/chr*.wig
#generateAlignability(
#  mapdir='/home/ruijie/data/hg19/map50/',
#  outdir=out,
#  athresh=4,
#  extension=134,
#  twoBitFile='/home/ruijie/data/hg19/hg19.2bit'
#)

bcfile=paste(out,'/temp.bc',sep="")
print(bcfile)
basealigncount(
  inputfile=infile,
  outputfile=bcfile,
  extension=134,
  filetype='bed',
  twoBitFile='/home/ruijie/data/hg19/hg19.2bit'
)

run.zinba(
  seq=infile,
  input="none",
  filetype="bed",
  twoBit="/home/ruijie/data/hg19/hg19.2bit",
  winSize=300,
  offset=75,
  extension=134,
  basecountfile=bcfile,
  align=out,
  selectmodel=F,
  formula= exp_count ~ gcPerc + align_perc + exp_cnvwin_log,
  formulaE= exp_count ~ gcPerc + align_perc + exp_cnvwin_log,
  formulaZ= exp_count ~ gcPerc + align_perc + exp_cnvwin_log,
  threshold=0.1,
  refinepeaks=1, 
  numProc=4,
  winGap=0,   
  FDR=TRUE,
  cleanup=TRUE,
  outfile=paste(out,"/",prefix,sep=""),
  printFullOut=1,
  method="mixture"
)
