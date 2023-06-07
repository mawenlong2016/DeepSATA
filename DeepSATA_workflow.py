# The step-by-step workflow of DeepSATA pipeline
# please make directory under DeepSATA/species/XXX
# then cd DeepSATA/species/XXX
# the thorough DeepSATA pipeline of building model for species XXX
# just run:
# python DeepSATA_workflow.py --ChrPos 1_ChrPos.txt --TFBind 2_TFBind.bed --Feature 3_Feature.bed --Reference pig --MEME 4_TF_meme.txt --Case pig --CPU 5
import argparse
from subprocess import *
import datetime
parser = argparse.ArgumentParser(description="Whole argparse")
parser.add_argument('-c','--ChrPos', default='1_ChrPos.txt',help='The chromosome position file')
parser.add_argument('-t','--TFBind', default='2_TFBind.bed',help='The open chromatin region file')
parser.add_argument('-f','--Feature', default='3_Feature.bed',help='The chromatin feature file')
parser.add_argument('-m','--MEME', default='4_TF_meme.txt',help='The transcriptional factor affinity file')
parser.add_argument('-u','--CPU', default='20',help='The cpu number')
parser.add_argument('-i','--Filename', default='example.bed',help='The input filename')
parser.add_argument('-r','--Reference', default='pig',help='The reference genome')
parser.add_argument('-s','--Case', default='pig',help='The case genome')
parser.add_argument('-g','--Flag', default='pig',help='The flag')

args = parser.parse_args()
print "Preparing input matrix begins, please wait..."
start1=datetime.datetime.now()
check_call(["Rscript  ../../code/1_setBin.R " + args.ChrPos],shell=True)
check_call(["python  ../../code/2_conserveBin.py " + args.ChrPos+".bin.bed " + args.TFBind],shell=True)
check_call(["python  ../../code/3_setBinLabel.py 2_conBin.bed " + args.Feature],shell=True)
check_call(["python  ../../code/4_splitBin.py"],shell=True)
check_call(["python   ../../code/5_splitLabel.py"],shell=True)
check_call(["Rscript  ../../code/6_extractSeq.R 4_trainBin.bed "+args.Reference],shell=True)
check_call(["Rscript  ../../code/6_extractSeq.R 4_validBin.bed "+args.Reference],shell=True)
check_call(["python ../../code/7_transH5.py train "+ args.MEME +"-"+args.Case + " " + args.CPU],shell=True)
check_call(["python ../../code/7_transH5.py valid "+ args.MEME +"-"+args.Case+ " " + args.CPU],shell=True)
end1=datetime.datetime.now()    
print "Finished preparation for input matrix, time elasped: %s", (end1 - start1)
print "Model training begins..."
start2=datetime.datetime.now()
check_call(["sh ../../code/8_trainDeepSATA.sh"],shell=True)
end2=datetime.datetime.now()    
print "Finished model training, time elasped", (end2 - start2)
print "Test set prediction begins..."
check_call(["python ../../code/9_predictDeepSATA.py 4_testBin.bed "+args.Reference+" test "+ args.MEME + "-" + args.Case + " " + args.CPU],shell=True)
print "Finished test prediction..."
print "Test set evaluation begins..."
check_call(["Rscript ../../code/10_evaluation.R"],shell=True)
print "Finished test evaluation..."




