# The script used for predicting TF and TFBS of given sequence
# Example: python 7_predictTF.py valid 4_TF_meme.txt 30
# Format of 4_TF_meme.txt just as follows:
# MOTIF MA0006.1 Ahr::Arnt
# letter-probability matrix: alength= 4 w= 6 nsites= 24 E= 0
# 0.125000  0.333333  0.083333  0.458333
# 0.000000  0.000000  0.958333  0.041667
# 0.000000  0.958333  0.000000  0.041667
# 0.000000  0.000000  0.958333  0.041667
# 0.000000  0.000000  0.000000  1.000000
# 0.000000  0.000000  1.000000  0.000000
# URL http://jaspar.genereg.net/matrix/MA0006.1
# Error: cannot allocate vector of size 1219.8 Gb

import sys
from subprocess import *
import multiprocessing as mp
import numpy as np
import os
import subprocess
import datetime

infilename=sys.argv[1]

if not os.path.exists(sys.argv[1]):
        os.makedirs(sys.argv[1])

def predictTF(pre, current, output, TF_meme, Fasta):
    output = output+"/temp_"+str(pre)+"_"+str(current)
    if not os.path.exists(output):
            os.makedirs(output)
    check_call(["sed -n "+str(pre)+","+str(current)+"p "+Fasta+" > "+output+"/temp.fasta"],shell=True)
    check_call(["fimo -oc  " + output + " " + TF_meme + " "+output+"/temp.fasta"],shell=True)

start1=datetime.datetime.now()
pool = mp.Pool(int(sys.argv[3]))

if 'vcf' in sys.argv[1]:
    totalSeq = int(subprocess.check_output("wc -l "+sys.argv[1] + ".fasta.ref.bed",shell=True).split()[0])
else:
    totalSeq = int(subprocess.check_output("wc -l 4_"+sys.argv[1] + "Bin.bed",shell=True).split()[0])
internal = int(totalSeq/int(sys.argv[3])) + 1

tep_position = np.arange(1,totalSeq,internal).tolist()
tep_position.append(totalSeq)
if len(tep_position)==1:
    tep_position.append(totalSeq)

if 'vcf' in sys.argv[1]:
    pre = 1
    for i in tep_position[1:len(tep_position)]:  
        i = i*14  
        pool.apply_async(predictTF, args=(pre, i, sys.argv[1], sys.argv[2], sys.argv[1] + ".fasta.ref.fasta",))
        pre = i + 1
else:
    pre = 1
    for i in tep_position[1:len(tep_position)]:  
        i = i*14  
        pool.apply_async(predictTF, args=(pre, i, sys.argv[1], sys.argv[2], "4_"+sys.argv[1] + "Bin.bed.wt1000.fasta",))
        pre = i + 1
   
pool.close()
pool.join()



if 'vcf' in sys.argv[1]:
    pre = 1
    for current in tep_position[1:len(tep_position)]: 
        current = current * 14
        num = int(subprocess.check_output("wc -l "+sys.argv[1]+"/temp_"+str(pre)+"_"+str(current) + "/fimo.tsv",shell=True).split()[0])-4
        check_call(["sed -n 2"+","+str(num)+"p  "+sys.argv[1]+"/temp_"+str(pre)+"_"+str(current) + "/fimo.tsv"+" >> "+sys.argv[1]+"/fimo_back.tsv"],shell=True)
        pre = current + 1
else:
    pre = 1
    for current in tep_position[1:len(tep_position)]: 
        current = current * 14
        num = int(subprocess.check_output("wc -l "+sys.argv[1]+"/temp_"+str(pre)+"_"+str(current) + "/fimo.tsv",shell=True).split()[0])-4
        check_call(["sed -n 2"+","+str(num)+"p  "+sys.argv[1]+"/temp_"+str(pre)+"_"+str(current) + "/fimo.tsv"+" >> "+sys.argv[1]+"/fimo_back.tsv"],shell=True)
        pre = current + 1

check_call(["bash ../../code/7_filterTF.sh " + sys.argv[1] +" " +sys.argv[2]+" " +sys.argv[3]],shell=True)

end1=datetime.datetime.now()    
print "Finished preparation for predict TF, time elasped: %s", (end1 - start1)





