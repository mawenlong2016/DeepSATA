# ___DeepSATA: A Deep Learning-Based Sequence Analyzer Incorporating Transcription Factor Affinity for Cross-Species Deciphering of Non-Coding Variant Effects___ <br>
The 'DeepSATA' pipeline is writen in python language, which can be used to leverage epigenomic data to predict the regulatory activity of genomic sequences, which is a promising tool for prioritizing functional single-nucleotide polymorphisms (SNPs) based on predicted chromatin effect signals. The effectiveness of same-species and cross-species predictions of DeepSATA has been demonstrated in predicting chromatin features on serval tissues of pigs (_Sus scrofa_), chickens (_Gallus gallus_), cattle (_Bos taurus_), and mice (_Mus musculus_). Furthermore, we showcased its effectiveness in analyzing pig genetic variants associated with economic traits. Overall, DeepSATA represents a valuable tool that broadens our knowledge of the potential influence of non-coding variants on complex diseases or traits across different species.
<br>
## DeepSATA dependencies <br>
DeepSATA denpends on the following software environments：<br>
1. [R](https://www.r-project.org) - The R (version 4.1) is needed. <br>
2. [Python](https://www.python.org) - The python (version 2.7) is needed. <br>
3. [Torch](http://torch.ch) - The torch (version 7) is needed. <br>
4. GPU should be supported. <br>
## DeepSATA usage <br>
### DeepSATA pipeline description: <br>
1. The thorough step-by-step pipeline of building a DeepSATA model has been described in **DeepSATA_workflow.py**:
```python
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

```
2. The thorough step-by-step pipeline of using DeepSATA model to make same-species or cross-species predicions has been described in **crossPredictDeepSATA.py**:
```python
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
```
### Example usage: <br>
1. Firstly, download the DeepSATA file
2. Secondly, prepare the following data just as the same format in **DeepSATA/species/example**
3. For training a new DeepSATA model

```
cd /DeepSATA/species/example
# training for pig
# example_1_ChrPos.txt represents the size of chromosome:
# chr1	1234567
# example_2_TFBind.bed represents the open chromatin regions:
# chr1	12345	13345	Duroc_Muscle
# example_3_Feature.bed represents the chromatin features:
# chr1	12345	13345	Duroc_Fat
# example_4_TF_meme.txt represents the transcription factor affinity
# letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
# 0.200000  0.800000  0.000000  0.000000
# 0.950000  0.000000  0.050000  0.000000
# 0.000000  1.000000  0.000000  0.000000
# 0.000000  0.000000  1.000000  0.000000
# 0.000000  0.000000  0.000000  1.000000
# 0.000000  0.000000  1.000000  0.000000

# the Reference presents the background species, here we focus on pig
# the Case presents the interested species, here we focus on pig
python ../../DeepSATA_workflow.py --ChrPos example_1_ChrPos.txt --TFBind example_2_TFBind.bed --Feature example_3_Feature.bed --Reference pig --MEME example_4_TF_meme.txt --Case pig --CPU 5
```
**Note**: The **4_TF_meme.txt** can be download from [JASPAR](https://jaspar.genereg.net). <br>
**Note**: The **example_2_TFBind.bed** and **example_3_Feature.bed** in **DeepSATA/species/example** have been compressed due to the file size limitation. <br>
**Note**: For exemplification, we set the **early stop epoch** as 30, the **maximum training epoch** as 500. These parameters can be modified at line 174 and 191 [**/DeepSATA/code/deepsata_train/4_train.lua**](https://github.com/mawenlong2016/DeepSATA/blob/main/code/deepsata_train/4_train.lua). <br>

4. For making same-species or cross-species predictions using DeepSATA:
```
cd /DeepSATA/species/cross
# for same-species predicion
# for BED file
python ../../crossPredictDeepSATA.py --Filename 4_example.bedBin.bed --Reference pig --MEME 4_TF_meme.txt --Case pig --CPU 5 --Flag example.bed
# for VCF file
python ../../crossPredictDeepSATA.py --Filename example.vcf --Reference pig --MEME 4_TF_meme.txt --Case pig --CPU 5 --Flag example.vcf
# for Fasta file
python ../../crossPredictDeepSATA.py --Filename 4_example.fasta --Reference pig --MEME 4_TF_meme.txt --Case pig --CPU 5 --Flag example.fasta

# for cross-species predicion
# for BED file
python ../../crossPredictDeepSATA.py --Filename 4_example.bedBin.bed --Reference pig --MEME 4_TF_meme.txt --Case mouse --CPU 5 --Flag example.bed
# for VCF file
python ../../crossPredictDeepSATA.py --Filename example.vcf --Reference pig --MEME 4_TF_meme.txt --Case mouse --CPU 5 --Flag example.vcf
# for FASTA file
python ../../crossPredictDeepSATA.py --Filename 4_example.fasta --Reference pig --MEME 4_TF_meme.txt --Case mouse --CPU 5 --Flag example.fasta
```
**Note**: Currently, only **pig**, **mouse**, **cattle**, and **chicken** are support. More species will be supported in the next updated DeepSATA. <br>
**Note**: Due to the file size limitation, we cannot upload the pre-trained best model for **pig**, **mouse**, **cattle**, and **chicken**, which are necessary in **cross-species prediction** of **DeepSATA/species/data**. Users can train the model by themselves following the instruction of **DeepSATA/DeepSATA_workflow.py**, or request these models via email: **mawenlong_nwsuaf@163.com**.<br>

## Ask questions
Please use [**DeepSATA/issues**](https://github.com/mawenlong2016/DeepSATA/issues) for how to use DeepSATA and reporting bugs.
