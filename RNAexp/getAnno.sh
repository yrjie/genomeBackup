if [ $# -lt 1 ]
then
    echo 'Auto script'
    exit
fi

cut -f4,7,10 GPL15931_Agilent_GEO_annotation_file.txt |awk 'BEGIN{FS="\t"} NF==3{ print}' >anno_a.tab
cut -f1,3,5 GPL15932_annotation.txt|awk -F"\t" 'NF==3{print}' >anno_x.tab

