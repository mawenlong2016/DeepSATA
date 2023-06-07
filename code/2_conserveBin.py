# The script used for conserving bins with ATAC signals
# Example: python 2_conserveBin.py 1_ChrPos.txt.bin.bed 2_TFBind.bed
# Format of 1_pigChrPos.txt.bin.bed just as follows:
# chr1    401     600
# chr1    601     800
# chr1    801     1000
# Format of 2_TFBind.bed just as follows:
# chr10   219886  219951  D8_beiji_duluoke
# chr10   369633  369704  D8_beiji_duluoke
# chr10   381624  381841  D8_beiji_duluoke

from subprocess import *
import sys
import subprocess
import os
import numpy as np

check_call(["bedtools intersect -a "+sys.argv[1]+" -b "+sys.argv[2]+" -wa > 2_conBin_tmp.bed"],shell=True)
check_call(["awk '!a[$1,$2,$3]++' 2_conBin_tmp.bed > 2_conBin.bed"],shell=True)
check_call(["rm 2_conBin_tmp.bed "],shell=True)
check_call(["perl -p -i -e 's/ /\t/g' 2_conBin.bed"],shell=True)
# check_call(["shuf 2_conBin.bed > 2_conBin-ori.bed"],shell=True)

# if not os.path.exi	sts("2_conBin/"):
#         os.makedirs("2_conBin/")

# num = int(subprocess.check_output("wc -l 2_conBin-ori.bed",shell=True).split()[0])
# split = 10
# internal = int(num/split) + 1
# tep_position = np.arange(1,num,internal).tolist()
# tep_position.append(num)
# pre = 1
# flag = 1
# for i in tep_position[1:len(tep_position)]:  
#     check_call(["sed -n "+str(pre)+","+str(i)+"p 2_conBin-ori.bed > 2_conBin/2_conBin_"+str(flag)+".bed"],shell=True)
#     pre = i + 1
#     flag = flag + 1
# check_call(["cp 2_conBin/2_conBin_1.bed 2_conBin.bed"],shell=True)
