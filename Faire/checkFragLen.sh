ls Fpeak/test/AtacGm12878_*Ex.bed|xargs -n1 sh -c 'windowBed -a $0 -b Fpeak/test/AtacGm12878_300+NoStL100Ex.bed -u -w 0|wc -l'
