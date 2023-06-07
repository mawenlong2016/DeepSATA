# The script used for splitting label of bins into training, validation, test sets
# Example: python 5_splitLabel.py
from subprocess import *

check_call(["grep -w 'chr1' 3_binLabel.bed | awk '{print $1,$2-400,$3+400,$4,$6}' | awk '{print $1"+'"_"'+"$2"+'"_"'+"$3"+'"\t"'+"$4"+'"\t"'+"$5}' " + ' > 5_validLabel.bed'],shell=True)
check_call(["grep -w 'chr2\|chr3' 3_binLabel.bed | awk '{print $1,$2-400,$3+400,$4,$6}' | awk '{print $1"+'"_"'+"$2"+'"_"'+"$3"+'"\t"'+"$4"+'"\t"'+"$5}' " + ' > 5_testLabel.bed'],shell=True)
check_call(["grep -w -v 'chr1\|chr2\|chr3' 3_binLabel.bed | awk '{print $1,$2-400,$3+400,$4,$6}' | awk '{print $1" + '"_"' +"$2"+'"_"'+"$3"+'"\t"'+"$4"+'"\t"'+"$5}' " + ' > 5_trainLabel.bed'],shell=True)
