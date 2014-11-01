cat data/FaireK562_c3.bed data/FaireK562_c4.bed data/FaireK562_c5.bed >temp.bed
windowBed -a ~/DNase/Faire/wgEncodeOpenChromFaireK562Pk.narrowPeak -b temp.bed -u -w 0 >Cneg.bed

cat data/FaireH1hesc_c2.bed data/FaireH1hesc_c3.bed data/FaireH1hesc_c4.bed >temp.bed
windowBed -a ~/DNase/Faire/wgEncodeOpenChromFaireH1hescPk.narrowPeak -b temp.bed -u -w 0 >Cneg.bed

cat data/FaireGm12878_c2.bed data/FaireGm12878_c3.bed data/FaireGm12878_c4.bed >temp.bed
windowBed -a ~/DNase/Faire/wgEncodeOpenChromFaireGm12878Pk.narrowPeak -b temp.bed -u -w 0 >Cneg.bed
