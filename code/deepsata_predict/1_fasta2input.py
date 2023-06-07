#Example usage: python 1_fasta2input.py infile.fasta window_size

#Convert sequences to hdf5 format. Used for processing variant format.
#It only supports fasta file produced by ../../deepsea_predict/6_extractSeq.R, which includes ref and alt
#allele information in sequence names. 

from Bio import SeqIO
import numpy as np
import sys
import h5py
import math
from os.path import basename,dirname
writeFasta=True
#writeh5=True
writeh5=False
#writebed=True
#writevcf=True
writebed=False
writevcf=False

#if vcf_original_allele_check is False, always use the ref and alt alleles user specified. 
#Otherwise ref allele must match the reference genome
vcf_original_allele_check=False#True

if len(sys.argv)==3:
    inputwindow = int(sys.argv[2])
else:
    inputwindow = 1000
mutpos=inputwindow/2-1
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
np.random.seed(1)




#write h5 file

def writeh5(seqs,filename,offset_values=None):
    seqsnp=np.zeros((len(seqs),4,1000),np.bool_)

    mydict={'A':np.asarray([1,0,0,0]),'G':np.asarray([0,1,0,0]),'C':np.asarray([0,0,1,0]),'T':np.asarray([0,0,0,1]),'N':np.asarray([0,0,0,0]),'H':np.asarray([0,0,0,0]),'a':np.asarray([1,0,0,0]),'g':np.asarray([0,1,0,0]),'c':np.asarray([0,0,1,0]),'t':np.asarray([0,0,0,1]),'n':np.asarray([0,0,0,0])}
    n=0
    for line in seqs:
        cline = line[ int(math.floor(( (len(line)-1000)/2.0))):int(math.floor(len(line)-(len(line)-1000)/2.0))]
        for c,i in zip(cline,range(len(cline))):
            seqsnp[n,:,i]=mydict[c]
        n=n+1
    
    #get the complementary sequences
    dataflip=seqsnp[:,::-1,::-1];
    seqsnp=np.concatenate([seqsnp, dataflip],axis=0)
    seqsnp = seqsnp.astype(np.uint8)
    f=h5py.File(filename,'w')
    f.create_dataset('testxdata', data= seqsnp,compression="gzip")
    f.close()

    
    
seqs=[str(fasta.seq) for fasta in fasta_sequences]

oris=[]
muts=[]
chrs=[]
poss=[]
annos=[]
names=[]
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
for fasta in fasta_sequences:
    anno = fasta.name.split('_')
    annos.append(fasta.name)
    oris.append(anno[0])
    muts.append(anno[1])
    chrs.append(anno[2])
    poss.append(int(anno[3])+mutpos)
    if len(anno)>5:
        names.append('_'.join(anno[5:]))

oris=np.asarray(oris)
muts=np.asarray(muts)
chrs=np.asarray(chrs)
seqs=np.asarray(seqs)
poss=np.asarray(poss)
annos=np.asarray(annos)
names=np.asarray(names)

seqsmut=[]
inds=[]
for i in range(len(seqs)):
    if vcf_original_allele_check:
        if type(oris[i]) is not np.string_ or seqs[i][mutpos:(mutpos+len(oris[i]))]!=oris[i]:
            continue
    else:
        if type(oris[i]) is not np.string_:
            continue
    inds.append(i)
    seqsmut.append(seqs[i][:mutpos]+ muts[i]+seqs[i][(mutpos+len(oris[i])):])
    seqs[i]=seqs[i][:mutpos]+ oris[i]+seqs[i][(mutpos+len(oris[i])):]

print "Number of valid variants:"
print(len(inds))
print "Number of input variants:"
print(len(seqs))
inds = np.asarray(inds)
seqsmut = np.asarray(seqsmut)
seqs = seqs[inds]
oris = oris[inds]
muts = muts[inds]
chrs = chrs[inds]
poss = poss[inds]
annos = annos[inds]
if len(names)>0:
    names = names[inds]

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
if writeFasta:
    allrecwt=[]
    allrecmut=[]
    for i in range(len(seqsmut)):
        allrecmut.append(SeqRecord(Seq(seqsmut[i]),id=annos[i]))
        allrecwt.append(SeqRecord(Seq(seqs[i]),id=annos[i]))
    SeqIO.write(allrecwt,open(sys.argv[1]+'.ref.fasta','w'),'fasta')
    SeqIO.write(allrecmut,open(sys.argv[1].replace('wt','mut')+'.ref.fasta','w'),'fasta')

if writeh5:
    writeh5(seqs,sys.argv[1]+'.ref.h5')
    writeh5(seqsmut,basename(sys.argv[1]).replace('wt','mut')+'.ref.h5')
    #writeh5(seqsmut,dirname(sys.argv[1])+'/'+basename(sys.argv[1]).replace('wt','mut')+'.ref.h5')

if writebed:
    myfile=open(sys.argv[1]+'.ref.bed','w')
    for i in range(len(seqs)):
        myfile.write(chrs[i]+'\t'+str(poss[i]-mutpos)+'\t'+str(poss[i]+mutpos+1)+'\t.\t0\t*\n')
    myfile.close()

if writevcf:
    myfile=open(sys.argv[1]+'.ref.vcf','w')
    if len(names)>0:
        for i in range(len(seqs)):
            myfile.write(chrs[i]+'\t'+str(poss[i])+'\t'+names[i]+'\t'+oris[i]+'\t'+muts[i]+'\n')
    else:
        for i in range(len(seqs)):
            myfile.write(chrs[i]+'\t'+str(poss[i])+'\t1\t'+oris[i]+'\t'+muts[i]+'\n')
    myfile.close()
    
        
    
