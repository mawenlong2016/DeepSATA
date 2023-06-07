# The script used for cross species predictions
# Example: python ../../code/11_crossPredictDeepSATA --Filename example.vcf --Reference pig --MEME example_4_TF_meme.txt --Case pig --CPU 5 --Flag example.vcf
# The output files will be in outdir.

from subprocess import *
from tempfile import mkdtemp
import sys
import os 

import argparse
parser = argparse.ArgumentParser(description="Whole argparse")
parser.add_argument('-m','--MEME', default='4_TF_meme.txt',help='The transcriptional factor affinity file')
parser.add_argument('-u','--CPU', default='20',help='The cpu number')
parser.add_argument('-i','--Filename', default='example.bed',help='The input filename')
parser.add_argument('-r','--Reference', default='pig',help='The reference genome')
parser.add_argument('-s','--Case', default='pig',help='The case genome')
parser.add_argument('-g','--Flag', default='pig',help='The flag')

args = parser.parse_args()
check_call(["cp ../../data/motif.enrich.name-"+args.Case+" ./motif.enrich.name"],shell=True)  
check_call(["cp ../../data/4_motif.names-"+args.Case+" ./4_motif.names"],shell=True)  
check_call(["cp ../../data/7_train.label.name-"+args.Case+" ./7_train.label.name"],shell=True) 


infilename=args.Filename
#for fasta input 
if infilename.endswith('fasta'):    
    check_call(["grep '>' "+args.Filename+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".wt1000.fasta.ref.bed"],shell=True)  
    check_call(["cp "+args.Filename+".wt1000.fasta.ref.bed 4_"+args.Flag+"Bin.bed"],shell=True)   
    check_call(["cp "+args.Filename+".wt1000.fasta 4_"+args.Flag+"Bin.bed.wt1000.fasta"],shell=True)   
    check_call(["python ../../code/7_predictTF.py "+  args.Flag +" ../../data/" + args.MEME+"-"+args.Case + " "+ args.CPU],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+args.Flag+" "+args.Filename+".wt1000.fasta" + " ../../data/" + args.MEME+"-"+args.Case],shell=True)
    check_call(["grep '>' "+args.Filename+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".wt1000.fasta.ref.bed"],shell=True)
    print "Successfully converted to input format"    
    check_call(["luajit  ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+args.Flag+".fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath ../../data/bestmodel.net"+"-"+args.Case],shell=True)
    print "Finished running DeepSATA. Now prepare output files..."
    check_call(["python  ../../code/deepsata_predict/3_h5ToOutput.py "+args.Filename+" 7_"+args.Flag+".fasta.ref.h5.pred.h5"],shell=True)
    print "Everything done."
#for bed input 
elif infilename.endswith('bed'):    
    check_call(["Rscript  ../../code/6_extractSeq.R "+args.Filename+" "+args.Reference],shell=True)
    check_call(["python ../../code/7_predictTF.py "+  args.Flag +" ../../data/" + args.MEME+"-"+args.Case + " "+ args.CPU],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+args.Flag+" "+args.Filename+".wt1000.fasta" + " ../../data/" + args.MEME+"-"+args.Case],shell=True)
    check_call(["grep '>' "+args.Filename+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".wt1000.fasta.ref.bed"],shell=True)
    print "Successfully converted to input format"    
    check_call(["luajit  ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+args.Flag+".fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath ../../data/bestmodel.net"+"-"+args.Case],shell=True)
    print "Finished running DeepSATA. Now prepare output files..."
    check_call(["python  ../../code/deepsata_predict/3_h5ToOutput.py "+args.Filename+" 7_"+args.Flag+".fasta.ref.h5.pred.h5"],shell=True)
    print "Everything done."
#for vcf input 
elif infilename.endswith('vcf'):
    check_call(["Rscript  ../../code/6_extractSeq.R "+args.Filename+" "+args.Reference],shell=True)        
    # retrieve 1100bp instead of 1000bp for supporting deletion variants (<100bp) 
    check_call(["python ../../code/deepsata_predict/1_fasta2input.py  "+args.Filename+".wt1000.fasta 1000"],shell=True) 
    check_call(["grep '>' "+args.Filename+".wt1000.fasta.ref.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".wt1000.fasta.ref.bed"],shell=True)
    check_call(["grep '>' "+args.Filename+".mut1000.fasta.ref.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".mut1000.fasta.ref.bed"],shell=True)
    check_call(["python ../../code/7_predictTF.py "+  args.Flag +".wt1000 ../../data/" + args.MEME+"-"+args.Case + " "+ args.CPU],shell=True)
    check_call(["python ../../code/7_predictTF.py "+  args.Flag +".mut1000 ../../data/" + args.MEME+"-"+args.Case + " "+ args.CPU],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+args.Flag+".wt1000 "+args.Filename+".wt1000.fasta.ref.fasta" + " ../../data/" + args.MEME+"-"+args.Case],shell=True) 
    check_call(["python  ../../code/7_seq2num.py "+args.Flag+".mut1000 "+args.Filename+".mut1000.fasta.ref.fasta" + " ../../data/" + args.MEME+"-"+args.Case],shell=True)
    print "Successfully converted to input format"


    check_call(["luajit ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+args.Filename+".wt1000.fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath ../../data/bestmodel.net"+"-"+args.Case],shell=True)
    check_call(["luajit ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+args.Filename+".mut1000.fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath ../../data/bestmodel.net"+"-"+args.Case],shell=True)
    print "Finished running DeepSEA. Now prepare output files..."
    check_call(["python ../../code/deepsata_predict/3_h5ToOutput.py "+args.Filename+" 7_"+args.Filename+".wt1000.fasta.ref.h5.pred.h5  7_"+args.Filename+".mut1000.fasta.ref.h5.pred.h5"],shell=True) 
    print "Everything done."


