# The script used for converting sequence into numeric value of hdf5 format
# Example: python 7_seq2num.py valid 4_validBin.bed.wt1000.fasta valid/fimo.tsv
# Format of 4_validBin.bed.wt1000.fasta just as follows:
# >__chr7_30508801_30509800
# AGCCAAGAGGTCTGGCTCCAGAATCCTCACTCATTAAATATGTATTTGCCTTGGCCTCCTTTCTTTCTTTTCTCATTTTT
# TCCTCGAGGTATAGTGGACACACAGCTCTGTATCAGTTTAGGGTGTGCAGCGCAATGATTTGACTTATCTACATCATGAA

from Bio import SeqIO
import numpy as np
import sys
import h5py
from subprocess import *
import multiprocessing as mp
from itertools import islice 
import datetime
import os
import subprocess
    

start1=datetime.datetime.now()
inputwindow=1000
mutpos=inputwindow/2-1
fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
np.random.seed(1)
check_call(["grep '>' "+sys.argv[2]+" | cut -c 2- | cut -d ' ' -f -1 > 4_"+sys.argv[1]+"seqBin.names"],shell=True)

if sys.argv[1]=="train":
    if os.path.exists("4_motif.names"):
        check_call(["rm 4_motif.names "],shell=True)
        
if not os.path.exists("4_motif.names"):
    check_call(["awk '!a[$1]++ {print $1}' "+sys.argv[1]+"/fimo.tsv"+" > 4_motif.names "],shell=True)
    # check_call(["awk '!a[$1]++ {print $1}' "+sys.argv[1]+"/fimo.tsv"+" | awk 'NR!=1' | grep -w '\w' > 4_motif.names "],shell=True)

tfdict={}
flag_motif=1
with open("4_motif.names", "r") as f:
    for motif in f.readlines():
        motif = motif.strip('\n')  
        tfdict[motif]=flag_motif
        flag_motif=flag_motif+1

tfPFM={}
motif_matrix = []
motif_flag = []
flag = 0
#with open("test_meme.txt", "r") as f:
with open(sys.argv[3], "r") as f:
    for motif in f.readlines():
        motif = motif.strip('\n')  
        motif = motif.split()
        if len(motif):
            motif_matrix.append(motif)
            if 'MOTIF' in motif:
                motif_flag.append(flag)
            flag = flag + 1

tfPFM={}
for flag in range(0,len(motif_flag)-1):
    motif_name = motif_matrix[motif_flag[flag]][1]
    temp_matrix = []
    for i in range(motif_flag[flag] + 2, motif_flag[flag+1]-1):
        temp_matrix.append(motif_matrix[i])
    tfPFM[motif_name] = temp_matrix

motif_name = motif_matrix[motif_flag[len(motif_flag)-1]][1]
temp_matrix = []
for i in range(motif_flag[len(motif_flag)-1] + 2, len(motif_matrix)-1):
    temp_matrix.append(motif_matrix[i])
tfPFM[motif_name] = temp_matrix

item_dict={'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3}
item_dict_a={'A':0,'C':2,'G':1,'T':3,'a':0,'c':2,'g':1,'t':3}
tfseq={}
seqBin=[]
flag2=0
with open("4_"+sys.argv[1]+"seqBin.names", "r") as f:
    for line in f.readlines():
        line = line.strip('\n')
        seqBin.append(line)
        tfseq[line]=flag2
        flag2=flag2+1


#write mat / h5 file
def writeh5(seqs,filename,offset_values=None):
    seqsnp=np.zeros((len(seqs),4,1000,flag_motif),np.bool_)
    seqsnp=seqsnp.astype(np.float64)
    mydict={'A':np.asarray([1,0,0,0]),'G':np.asarray([0,1,0,0]),'C':np.asarray([0,0,1,0]),'T':np.asarray([0,0,0,1]),'N':np.asarray([0,0,0,0]),'H':np.asarray([0,0,0,0]),'a':np.asarray([1,0,0,0]),'g':np.asarray([0,1,0,0]),'c':np.asarray([0,0,1,0]),'t':np.asarray([0,0,0,1]),'n':np.asarray([0,0,0,0])}
    n=0
    offset_values = np.zeros(len(seqs))
    for line,o in zip(seqs,offset_values):  
        if len(line)<1000:
            print len(line)
            continue
            raise Exception("Each fasta sequence has to be at least 1000bp.")
        #if the user specified region/sequence is longer than 1000bp, use the center 1000bp
        cline = line[((int(o) + (len(line)-1000)/2)):(int(o)+(len(line)+1000)/2)]
        for c,i in zip(cline,range(len(cline))):
            seqsnp[n,:,i,0]=mydict[c] 
        n=n+1

    if 'wt1000' in  sys.argv[2]:
        with open(sys.argv[1]+"/fimo.tsv","r") as f:
            for line in f.readlines(): #islice(f,1,int(subprocess.check_output("wc -l "+sys.argv[1]+"/fimo.tsv",shell=True).split()[0])-4):           
                tf_name,_,seq_name,start,end,_,_,_,_,motif_seq = line.split('\t')
                motif_seq = motif_seq.strip('\n')
                temp_matrix=np.array([0,0,0,0]).reshape(4,1)
                flag=0
                for item in motif_seq:
                    if item in ['A','C','G','T','a','c','g','t']:
                        item_temp = [0,0,0,0]
                        item_temp[item_dict_a[item]] = tfPFM[tf_name][flag][item_dict[item]]
                        item_temp = np.array(item_temp).reshape(4,1)
                        temp_matrix=np.concatenate([temp_matrix, item_temp], axis=1)
                        flag=flag+1   
                    else:
                        temp_matrix = np.array([-1,-1,-1,-1]).reshape(4,1)
                        flag=flag+1      
                seqsnp[tfseq[seq_name],:,(int(start)-1):int(end),tfdict[tf_name]]=temp_matrix[:,1:temp_matrix.shape[1]]
 
    elif 'mut1000' in  sys.argv[2]:
        with open(sys.argv[1]+"/fimo.tsv","r") as f:
            for line in f.readlines(): #islice(f,1,int(subprocess.check_output("wc -l "+sys.argv[1]+"/fimo.tsv",shell=True).split()[0])-4):           
                tf_name,_,seq_name,start,end,_,_,_,_,motif_seq = line.split('\t')
                motif_seq = motif_seq.strip('\n')
                temp_matrix=np.array([0,0,0,0]).reshape(4,1)
                flag=0
                for item in motif_seq:
                    if item in ['A','C','G','T','a','c','g','t']:
                        item_temp = [0,0,0,0]
                        item_temp[item_dict_a[item]] = tfPFM[tf_name][flag][item_dict[item]]
                        item_temp = np.array(item_temp).reshape(4,1)
                        temp_matrix=np.concatenate([temp_matrix, item_temp], axis=1)
                        flag=flag+1   
                    else:
                        temp_matrix = np.array([-1,-1,-1,-1]).reshape(4,1)
                        flag=flag+1      
                seqsnp[tfseq[seq_name],:,(int(start)-1):int(end),tfdict[tf_name]]=temp_matrix[:,1:temp_matrix.shape[1]]
 
    else:
        with open(sys.argv[1]+"/fimo.tsv","r") as f:
            for line in f.readlines(): #islice(f,1,int(subprocess.check_output("wc -l "+sys.argv[1]+"/fimo.tsv",shell=True).split()[0])-4):           
                tf_name,_,seq_name,start,end,_,_,_,_,motif_seq = line.split('\t')
                motif_seq = motif_seq.strip('\n')
                temp_matrix=np.array([0,0,0,0]).reshape(4,1)
                flag=0
                for item in motif_seq:
                    if item in ['A','C','G','T','a','c','g','t']:
                        item_temp = [0,0,0,0]
                        item_temp[item_dict_a[item]] = tfPFM[tf_name][flag][item_dict[item]]
                        item_temp = np.array(item_temp).reshape(4,1)
                        temp_matrix=np.concatenate([temp_matrix, item_temp], axis=1)
                        flag=flag+1   
                    else:
                        temp_matrix = np.array([-1,-1,-1,-1]).reshape(4,1)
                        flag=flag+1      
                seqsnp[tfseq[seq_name],:,(int(start)-1):int(end),tfdict[tf_name]]=temp_matrix[:,1:temp_matrix.shape[1]]
    
    dataflip=seqsnp[:,::-1,::-1,::-1]
    seqsnp=np.concatenate([seqsnp, dataflip],axis=0)
    f=h5py.File(filename,'w')
    f.create_dataset('testxdata', data= seqsnp,compression="gzip")
    f.close()    
    
seqs=[str(fasta.seq) for fasta in fasta_sequences]
if True:
    writeh5(seqs,'7_'+sys.argv[1]+'.fasta.ref.h5')
    end1=datetime.datetime.now()    
    print "Finished preparation for input matrix, time elasped: %s", (end1 - start1)




