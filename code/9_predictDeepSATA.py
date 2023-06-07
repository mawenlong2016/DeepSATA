# The script used for test set prediction
# Example: python ./code/9_predictDeepSEA.py 4_testBin.bed ./predict_results test 4_TF_meme.txt 30
# The output files will be in outdir.

from subprocess import *
from tempfile import mkdtemp
import sys
import os 

infilename=sys.argv[1]
#for fasta input 
if infilename.endswith('fasta'):    
    check_call(["grep '>' "+sys.argv[1]+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+sys.argv[1]+".wt1000.fasta.ref.bed"],shell=True)  
    check_call(["cp "+sys.argv[1]+".wt1000.fasta.ref.bed 4_"+sys.argv[3]+"Bin.bed"],shell=True)   
    check_call(["cp "+sys.argv[1]+".wt1000.fasta 4_"+sys.argv[3]+"Bin.bed.wt1000.fasta"],shell=True)   
    check_call(["python ../../code/7_predictTF.py "+  sys.argv[3] +" " + sys.argv[4] + " "+ sys.argv[5]],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+sys.argv[3]+" "+sys.argv[1]+".wt1000.fasta" + " " + sys.argv[4]],shell=True)
    check_call(["grep '>' "+sys.argv[1]+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+sys.argv[1]+".wt1000.fasta.ref.bed"],shell=True)
    print "Successfully converted to input format"    
    check_call(["luajit  ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+sys.argv[3]+".fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath ./train_results/bestmodel.net "],shell=True)
    print "Finished running DeepSATA. Now prepare output files..."
    check_call(["python  ../../code/deepsata_predict/3_h5ToOutput.py "+sys.argv[1]+" 7_"+sys.argv[3]+".fasta.ref.h5.pred.h5"],shell=True)
    print "Everything done."
#for bed input 
elif infilename.endswith('bed'):    
    check_call(["Rscript  ../../code/6_extractSeq.R "+sys.argv[1]+" "+sys.argv[2]],shell=True)
    check_call(["python ../../code/7_predictTF.py "+  sys.argv[3] +" " + sys.argv[4] + " "+ sys.argv[5]],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+sys.argv[3]+" "+sys.argv[1]+".wt1000.fasta" + " " + sys.argv[4]],shell=True)
    check_call(["grep '>' "+sys.argv[1]+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+sys.argv[1]+".wt1000.fasta.ref.bed"],shell=True)
    print "Successfully converted to input format"    
    check_call(["luajit  ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+sys.argv[3]+".fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath ./train_results/bestmodel.net "],shell=True)
    print "Finished running DeepSATA. Now prepare output files..."
    check_call(["python  ../../code/deepsata_predict/3_h5ToOutput.py "+sys.argv[1]+" 7_"+sys.argv[3]+".fasta.ref.h5.pred.h5"],shell=True)
    print "Everything done."
#for vcf input 
elif infilename.endswith('vcf'):
    check_call(["Rscript  ../../code/6_extractSeq.R "+sys.argv[1]+" "+sys.argv[2]],shell=True)        
    # retrieve 1100bp instead of 1000bp for supporting deletion variants (<100bp) 
    check_call(["python ../../code/deepsata_predict/1_fasta2input.py  "+sys.argv[1]+".wt1000.fasta 1000"],shell=True) 
    check_call(["grep '>' "+sys.argv[1]+".wt1000.fasta.ref.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+sys.argv[1]+".wt1000.fasta.ref.bed"],shell=True)
    check_call(["grep '>' "+sys.argv[1]+".mut1000.fasta.ref.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+sys.argv[1]+".mut1000.fasta.ref.bed"],shell=True)
    check_call(["python ../../code/7_predictTF.py "+  sys.argv[3] +".wt1000 " + sys.argv[4] + " "+ sys.argv[5]],shell=True)
    check_call(["python ../../code/7_predictTF.py "+  sys.argv[3] +".mut1000 " + sys.argv[4] + " "+ sys.argv[5]],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+sys.argv[3]+".wt1000 "+sys.argv[1]+".wt1000.fasta.ref.fasta" + " " + sys.argv[4]],shell=True) 
    check_call(["python  ../../code/7_seq2num.py "+sys.argv[3]+".mut1000 "+sys.argv[1]+".mut1000.fasta.ref.fasta" + " " + sys.argv[4]],shell=True)
    print "Successfully converted to input format"


    check_call(["luajit ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+sys.argv[1]+".wt1000.fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath './train_results/bestmodel.net' "],shell=True)
    check_call(["luajit ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+sys.argv[1]+".mut1000.fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath './train_results/bestmodel.net' "],shell=True)
    print "Finished running DeepSEA. Now prepare output files..."
    check_call(["python ../../code/deepsata_predict/3_h5ToOutput.py "+sys.argv[1]+" 7_"+sys.argv[1]+".wt1000.fasta.ref.h5.pred.h5  7_"+sys.argv[1]+".mut1000.fasta.ref.h5.pred.h5"],shell=True) 
    print "Everything done."


