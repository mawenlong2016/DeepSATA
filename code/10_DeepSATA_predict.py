# The script used for DeepSATA prediction
# The input file should be in: VCF, BED, or Fasta format
# VCF:
# chr5  66103958    .   G   A
# BED:
# chr5  66103459    66104458
# Fasta:
# >G_A_chr5_66103459_66104458_.
# GGACCCCTGTGCCCTGTAGCCCCAGAACTGGGGCTCACGGCGTGGGGTGGGGAGTGGCTGAACATGGCAGCCAAGAGGAG
# CAGAGGCTGTGGAAGCCAGGCCAGCCCCCGGCCACTGCCCAGCTCTCCTTGGGGTTGGTCCAAAGGAGGTGACTTGAGTG
# GTCCAAAGCCAGGTCCAGCCTGGGGGACCCAGGGCAAGTCACTTTCCCTTTGACCGCACCCCCCCCCCCAGGCTTCCGCC
# ACTGCCTGGATCCCAGGTGTCCAAGAGGAGGGTGGCCAACTAACTATCAACTAACTAACTAGTTTTTCTGTCTTCCACTT
# GCAGCTCAGCACCAGGTCATGCATGGGTCGTGACCAGCACTTGTCCTTGTTGTTGCTGAACTGAACAAAACAGAATAGAA
# AGTATCAGGAAGTCATCAGCGAAGCGGAGGAGTGAATACTCTATTGAGCAACTTTTGTATCGAACTGAGGTACCCTGGCA
# GCTCACTGGAGTTGTGCTAGAAAGTCTGTTTCTTCCTTTACAGCCTAGGTCCAAAGGCTGAAATCCCAGAAGCCTTTCAA
# CTGGCTTTCACTTAAGGCTCTGGAGTAAGAAGGACCCCCAAGAAAAAGAGTGAGGAAGAGACACCAGAACATAACACCGC
# TTCCCTCCACCAGGCAGCAAAGTCACCAGGCATCTCTCTCCATCAGGCTGCTCCACACAGAGGCCCCTGGGGCCGCCCTC
# CGAAACCCCCACCCCCCCGACTGGGAGAGGAACTGTCGAGAAACACTGGAGGAAAACATCTGAGGGCTGGCGTCCTTGTA
# CTCGCTCTCAGTAGCGAACTTCCACGCTCACCCGCCGAAAGTTCCACCCAGCCTGTGCAGCCCCTTCTGGAAAGCTTACA
# CTTTGGCTGTGGAACTGGACTGAGAAATAGTCACCATATCTACTTGGTCGAGGGGGGTTTAAAAGAAAAAAATCAGACTT
# CATCAGTGTGTAACCATCATCAAGATCACGAGCATCCAGG

# Example usage:
# python DeepSATA_predict.py --Filename example.bed --Reference pig --Case pig --Flag predict --MEME 4_TF_meme.txt 

python ../../code/9_predictDeepSATA.py --Filename 4_testBin.bed --Reference pig --Flag test --MEME 4_TF_meme.txt --Case  pig --CPU 5

# The output files will be in outdir:

import argparse
from subprocess import *
from tempfile import mkdtemp
import sys
import os 

parser = argparse.ArgumentParser(description="Whole argparse")
parser.add_argument('-f','--Filename', default='./data/example.bed',help='the input filename')
parser.add_argument('-r','--Reference', default='pig',help='the reference genome')
parser.add_argument('-c','--Case', default='pig',help='the case genome')
parser.add_argument('-u','--CPU', default='pig',help='the cpu number')
parser.add_argument('-g','--Flag', default='pig',help='the flag')
parser.add_argument('-m','--MEME', default='4_TF_meme.txt',help='The transcriptional factor affinity file')


infilename=args.Filename

python ../../code/9_predictDeepSEAT.py 4_testBin.bed ./predict_results test 4_TF_meme.txt 20

#for fasta input 
if infilename.endswith('fasta'):
    try:
        tempdir = mkdtemp()
        check_call(['cp',infilename,tempdir+'/infile.fasta'])
        print "Successfully copied input to working directory."
        check_call(["grep '>'  "+tempdir+'/infile.fasta '+">"+ tempdir+'/infile.fasta.name'],shell=True)
        try:
            check_call(["python 1_fasta2input.nomut.py  "+tempdir+"/infile.fasta"],shell=True)
        except:
            raise Exception("Fasta format error.")
        print "Successfully converted to input format"
        #check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.fasta.ref.h5"],shell=True)
        check_call(["luajit 2_DeepSEA.lua -test_file_h5 "+tempdir+"/infile.fasta.ref.h5"+"-threads 10 -type 'cuda' -netPath './bestmodel.net' "],shell=True)
        print "Finished running DeepSEA. Now prepare output files..."
        check_call(["python 3_h5ToOutput.py "+tempdir+"/infile.fasta "+tempdir+"/infile.fasta.ref.h5.pred.h5"],shell=True) 
        if cpoutdir:
            check_call(["cp", tempdir+"/infile.fasta.out" ,outdir ])
    finally:
        print "Finished creating output file. Now clean up..."
        call(['rm',tempdir,'-r'])
        print "Everything done."
#for bed input 
elif infilename.endswith('bed'):
    check_call(["Rscript  ../../code/6_extractSeq.R "+infilename+" "+args.Reference],shell=True)
    check_call(["python ../../code/7_predictTF.py "+args.Flag +" "+ args.MEME+"-"+args.Case + " "+ args.CPU],shell=True)
    check_call(["python  ../../code/7_seq2num.py "+args.Case+" "+args.Filename+".wt1000.fasta"],shell=True) 
    check_call(["grep '>' "+args.Filename+".wt1000.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".wt1000.fasta.ref.bed"],shell=True)
    print "Successfully converted to input format"    
    check_call(["luajit  ../../code/deepsata_predict/2_DeepSATA.lua -test_file_h5 7_"+args.Flag+".fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath './train_results/*/bestmodel.net' "],shell=True)
    print "Finished running DeepSATA. Now prepare output files..."
    check_call(["python  code/deepsata_predict/3_h5ToOutput.py "+args.Filename+" 7_"+args.Flag+".fasta.ref.h5.pred.h5"],shell=True)      
    print "Everything done."
    
#for vcf input 
elif infilename.endswith('vcf'):
    try:
        check_call(["Rscript  code/6_extractSeq.R "+args.Filename+" "+args.Reference],shell=True)        
    except:
        raise Exception('Vcf format error.')
    #retrieve 1000bp instead of 1000bp for supporting deletion variants (<100bp) 
    check_call(["python code/deepsea_predict/1_fasta2input.py  "+args.Filename+".wt1000.fasta 1000"],shell=True) 
    check_call(["grep '>' "+args.Filename+".wt1000.fasta.ref.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".wt1000.fasta.ref.bed"],shell=True)
    check_call(["grep '>' "+args.Filename+".mut1000.fasta.ref.fasta |sed 's/>__//g'|sed 's/_/\t/g' > "+args.Filename+".mut1000.fasta.ref.bed"],shell=True)
    check_call(["python code/7_predictTF.py "+  args.Case +".wt1000 " + "4_TF_meme_"+args.Case+".txt" + " "+ args.CPU],shell=True)
    check_call(["python code/7_predictTF.py "+  args.Case +".mut1000 " + "4_TF_meme_"+args.Case+".txt" + " "+ args.CPU],shell=True)
    check_call(["python  code/7_seq2num.py "+args.Case+".wt1000 "+args.Filename+".wt1000.fasta.ref.fasta"],shell=True) 
    check_call(["python  code/7_seq2num.py "+args.Case+".mut1000 "+args.Filename+".mut1000.fasta.ref.fasta"],shell=True)    

    print "Successfully converted to input format"
    #check_call(["luajit code/deepsea_predict/2_DeepSEA.lua -test_file_h5 "+args.Filename+".wt1000.fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath './train_results/model,L1Sparsity=1e-08,save=train_results/bestmodel.net' "],shell=True)
    #check_call(["luajit code/deepsea_predict/2_DeepSEA.lua -test_file_h5 "+args.Filename+".mut1000.fasta.ref.h5"+" -threads 10 -type 'cuda' -netPath './train_results/model,L1Sparsity=1e-08,save=train_results/bestmodel.net' "],shell=True)
    #print "Finished running DeepSEA. Now prepare output files..."
    #check_call(["python code/deepsea_predict/3_h5ToOutput.py "+args.Filename+" "+args.Filename+".wt1000.fasta.ref.h5.pred.h5  "+args.Filename+".mut1000.fasta.ref.h5.pred.h5"],shell=True) 
    print "Everything done."

