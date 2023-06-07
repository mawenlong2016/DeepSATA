# The script used for splitting bins into training, validation, test sets
# Example: python 4_splitBin.py

from subprocess import *
# sequence defined as surrounding bins with size of 1000bp(+-400bp)
# chr1 as validation set, chr2 and chr3 as test set, the others as training set
check_call(["grep -w 'chr1' 2_conBin.bed | awk '{print $1"+'"\t"'+"$2-400"+'"\t"'+"$3+400}'  > 4_validBin.bed"],shell=True)
check_call(["grep -w 'chr2\|chr3' 2_conBin.bed | awk '{print $1"+'"\t"'+"$2-400"+'"\t"'+"$3+400}' > 4_testBin.bed"],shell=True)
check_call(["grep -w -v 'chr1\|chr2\|chr3' 2_conBin.bed | awk '{print $1"+'"\t"'+"$2-400"+'"\t"'+"$3+400}' > 4_trainBin.bed"],shell=True)

