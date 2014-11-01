args=commandArgs()
if (length(args)<6){
print("Usage: R --no-save --slave --args file1 file2 <getCorr.R")
q()
}
a=as.matrix(read.table(args[5]))
b=as.matrix(read.table(args[6]))
print(cor(a,b))
