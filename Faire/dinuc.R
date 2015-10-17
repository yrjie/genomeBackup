dinuc = function(output, twoBit="hg18.2bit", win=200,chrm="chr21"){
	if(!file.exists(twoBit)) stop("twobit file does not exist")
	.C("gcwin",as.character(output),as.character(twoBit), as.integer(win), as.character(chrm), package="zinba")
}

