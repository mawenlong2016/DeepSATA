# The script used for transforming input information into numeric feature matrix of H5 format
# Example: python 7_transH5.py valid JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt
# It contains three steps:
# First, predicting TF and TFBS of given sequence
# Second, converting sequence into numeric value of hdf5 format
# Third, combining feature matrix and label into a single H5 file

import sys
from subprocess import *
import datetime

start1=datetime.datetime.now()
check_call(["python ../../code/7_predictTF.py " + sys.argv[1] + " " + sys.argv[2] + " " + sys.argv[3]],shell=True)
check_call(["python ../../code/7_seq2num.py " + sys.argv[1] + " 4_" + sys.argv[1] + "Bin.bed.wt1000.fasta"+ " " + sys.argv[2]],shell=True)
check_call(["Rscript ../../code/7_labelH5.R " + sys.argv[1] + " 7_" + sys.argv[1] + '.fasta.ref.h5' + " 5_"+sys.argv[1]+"Label.bed"],shell=True)
end1=datetime.datetime.now()    
print "Finished preparation for input matrix, time elasped: %s", (end1 - start1)
