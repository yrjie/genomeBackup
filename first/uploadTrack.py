import sys,os
import getpass
import json
from urllib import quote

if len(sys.argv)<2:
    print 'Usage: python uploadTrack.py <cfg file> [create]'
    print '\nformat of cfg (tab-delimited):\n\tfile format asm lib_name shortLabel longLabel tags BCS'
    print 'example (tab-delimited):\n\t/tmp/test.bw bigWig hg19 TestLib_0 test1S test1L {"assembly":"hg19"} http://biogpu.d1.comp.nus.edu.sg:8003/'
    print '\nNote:'
    print '1. the file should be accessible to BCS'
    print '2. supported formats: bed3-12, bigWig, bigBed, vcf, psl, bam, encodeNarrowPeak, encodeBroadPeak, chiapet.itx, pair.cluster'
    print '3. multiple lines are allowed'
    print '4. create is the option for creating tracks\n\t(values: [0 - NO track creation, 1 - create the tracks], default: 1)'
    print '5. all lines starting with "#" will be ignored'
    print ''
    exit(1)

user=raw_input('username: ')
pwd = getpass.getpass('password: ')

if len(sys.argv)>2:
    create=sys.argv[2]
else:
    create='1'

fi=open(sys.argv[1])
for line in fi:
    line=line.strip()
    if len(line)<1 or line[0]=='#':
    	continue
    temp=line.split('\t')
    file=temp[0]
    trkName=os.path.basename(temp[0]).replace('.','_')
    trkName=quote(trkName)
    format=temp[1]
    asm=temp[2]
    lib=temp[3]
    shortL=temp[4]
    longL=temp[5]
    tags=temp[6]
    bcs=temp[7].strip('/')
    if (format == 'bed'): format='bed12'
    if 'bed' in format or 'vcf' in format or 'psl' in format or 'encode' in format or 'chiapet.itx' in format or 'pair.cluster' in format:
    	driver='mysql_table'
    elif 'bigWig' in format:
    	driver='bigWig'
    elif 'bam' in format:
    	driver='bam'
    elif 'bigBed' in format:
    	driver='bigBed'
    else:
    	print "Unknown format"
	print "Supported formats are:"
	print "  bed3-12, bigWig, bigBed, vcf, psl, bam, encodeNarrowPeak, encodeBroadPeak, chiapet.itx, pair.cluster"
    	sys.exit(1)
    t2 = json.loads(tags)
    t2['name'] = os.path.basename(file)
    tags = json.dumps(t2,separators=(',', ':'))
    hook='http://biogpu.d1.comp.nus.edu.sg/basic2/reqbin/?asm='+asm+'&lib='+lib+'&track='+trkName+'&shortLabel='+shortL+'&longLabel='+longL+'&user='+user+'&bcs='+bcs+'&create='+create
    #hook='http://biogpu.d1.comp.nus.edu.sg:8001/reqbin/?asm='+asm+'&lib='+lib+'&track='+trkName+'&shortLabel='+shortL+'&longLabel='+longL+'&user='+user+'&bcs='+bcs+'&create='+create
    if file.startswith('http') or file.startswith('ftp'):
	if driver!='bigWig' and driver!='bigBed':
	    raise Exception('only bigWig and bigBed are supported for url')
	cmd='curl --user '+user+':'+pwd+' -F \'tags='+tags+'\' -F \'format='+format+'\' -F \'input='+file+'\' -F \'importmode=remote\' -F \'hook='+hook+'\' -X POST '+bcs+'/create/'+driver+'/'
    else:
    	cmd='curl --user '+user+':'+pwd+' -F \'tags='+tags+'\' -F \'format='+format+'\' -F \'input=upload://xxxfile\' -F xxxfile=@\"' + file + '\" -F \'hook='+hook+'\' -X POST '+bcs+'/create/'+driver+'/'
    	#cmd='curl --user '+user+':'+pwd+' -F \'tags='+tags+'\' -F \'format='+format+'\' -F \'input=file://'+file+'\' -F \'hook='+hook+'\' -X POST '+bcs+'/create/'+driver+'/'
#    curL -F 'tags={"assembly":"hg19"}' -F 'format=bigWig' -F 'importmode=remote' -F 'input=http://www.somewhere.com/path/file.bigWig' -X POST http://biogpu.ddns.comp.nus.edu.sg:8002/create/bigWig/
#    print driver
    print cmd
    os.system(cmd)
fi.close()
