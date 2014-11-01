root dir
/data1/BASIC/basic

mapping
vim engine/fileformat_defs.yaml
sudo /usr/sbin/httpd -k restart

data
/data3/shaojiang/ChIA-PET

create table
./python console/table_util.py create -l 'Huvec-chr_mark' -i /data3/shaojiang/ChIA-PET/wgEncodeBroadHistoneHuvecH3k27acStdPk.broadPeak -p bed12 hg19 'wgEncodeBroadHistoneHuvecH3k27acStdPk'

create track
./python console/track_util.py new -n 'wgEncodeBroadHistoneHuvecH3k27acStdPk' 189 trfac

install info
http://biogpu.ddns.comp.nus.edu.sg/wiki/index.php/BASIC_Miscellaneous_Topics


south Africa
musa
