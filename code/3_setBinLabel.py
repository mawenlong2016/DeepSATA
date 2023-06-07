# The script used for setting bins label from feature information
# Example: python 3_setBinLabel.py 2_conBin.bed 3_Feature.bed
# Format of 2_conBin.bed just as follows:
# chr1    4001    4200
# chr1    4201    4400
# chr1    4601    4800
# Format of 3_Feature.bed just as follows:
# chr10   219886  219951  D8_beiji_duluoke
# chr10   369633  369704  D8_beiji_duluoke
# chr10   381624  381841  D8_beiji_duluoke

from subprocess import *
import sys
import os
check_call(["bedtools intersect -a " + sys.argv[1]+" -b "+sys.argv[2]+" -wo > 3_conBin_tmp.bed"],shell=True)
check_call(["awk '{print $1"+'"\t"$2'+'"\t"$3'+'"\t"$7' \
		+'"\t"$8+1 }' +"' 3_conBin_tmp.bed" \
		+" | awk '{a[$1"+'"\t"$2'+'"\t"$3' \
		+'"\t"$4]+=$5}END{for (i in a) print i'+'"\t"a[i]}' \
		+"' > 3_conBin.bed"],shell=True)
# assign 1 if half of 200bp overlap with feature signal else 0
check_call(["awk '$5> 100 {print $1"+'"\t"'+"$2"+'"\t"'+"$3"+'"\t"'+"$4"+'"\t"'+"$5"+'"\t"'+"1}' 3_conBin.bed > 3_binLabel.bed"],shell=True)
check_call(["awk  '$5<= 100 {print $1"+'"\t"'+"$2"+'"\t"'+"$3"+'"\t"'+"$4"+'"\t"'+"$5"+'"\t"'+"0}' 3_conBin.bed >> 3_binLabel.bed"],shell=True)
check_call(["rm 3_conBin_tmp.bed"],shell=True)
check_call(["sort -k5nr 3_binLabel.bed | sort -V -k1 -k2n > 3_binLabel_tmp.bed"],shell=True)
check_call(["mv 3_binLabel_tmp.bed 3_binLabel.bed"],shell=True)
#check_call(["perl -p -i -e 's/ /\t/g' 3_binLabel.bed"],shell=True)
check_call(["awk '!a[$4]++ {print $4}' 3_binLabel.bed > label_name.txt"],shell=True)
