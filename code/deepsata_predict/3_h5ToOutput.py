import numpy as np
import h5py
import pandas as pd
import gzip
import sys
import os
from subprocess import *
#from statsmodels.distributions import ECDF
import joblib



evoinfo_exist=os.path.isfile('./resources/Gerp/gerp_elements.tsv.gz') and os.path.isfile('./resources/Gerp/gerp_scores.tsv.gz') and  os.path.isfile('./resources/phastCons/primates_nohuman.tsv.gz') and os.path.isfile('./resources/phyloP/primates_nohuman.tsv.gz') 

#header = np.loadtxt('./resources/predictor.names',dtype=np.str)
header = np.loadtxt('7_train.label.name',dtype=np.str)

#load  column headers
if sys.argv[1].endswith('.vcf'):
    #coordata = pd.read_csv(sys.argv[1],header=None,delimiter='\t').iloc[:,np.asarray([0,1,2,3,4])]
    #coordata = pd.read_csv(sys.argv[1]+'.wt1000.fasta.ref.vcf',header=None,delimiter='\t').iloc[:,np.asarray([0,1,2,3,4])]
    coordata = pd.read_csv(sys.argv[1]+'.wt1000.fasta.ref.bed',header=None,delimiter='\t').iloc[:,np.asarray([2,4,5,0,1])]
    coordata[0]=coordata[0].str[1:]
    coordata[4]=coordata[4]-500
    coordata.columns = ['chr','pos','name','ref','alt']
    coordata.pos=coordata.pos.astype(int)
    wfile1 = open(sys.argv[1]+'.out.ref','w')
    wfile2 = open(sys.argv[1]+'.out.alt','w')
    wfile3 = open(sys.argv[1]+'.out.logfoldchange','w')
    wfile4 = open(sys.argv[1]+'.out.diff','w')
    # wfile6 = open(sys.argv[1]+'.out.evalue','w')    
    if evoinfo_exist:
        wfile5 = open(sys.argv[1]+'.out.snpclass','w')
        wfile7 = open(sys.argv[1]+'.out.summary','w')
        wfile8 = open(sys.argv[1]+'.out.funsig','w')

if sys.argv[1].endswith('.bed') or sys.argv[1].endswith('.fasta'):
    coordata = pd.read_csv(sys.argv[1]+'.wt1000.fasta.ref.bed',header=None,delimiter='\t').iloc[:,:3]
    coordata.columns = ['chr','start','end']
    coordata.start=coordata.start.astype(int)
    coordata.end=coordata.end.astype(int)
    wfile = open(sys.argv[1]+'.out','w')

# if sys.argv[1].endswith('.fasta'):
#     coordata = pd.read_csv(sys.argv[1]+'.name',header=None,delimiter='\t')
#     coordata.columns = ['name']
#     wfile = open(sys.argv[1]+'.out','w')


#process output files
if sys.argv[1].endswith('.bed') or sys.argv[1].endswith('.fasta'):
    #write header
    #write data
    data=np.asarray(h5py.File(sys.argv[2],'r')['pred'])
    data=data[:(data.shape[0]/2),:]/2.0+data[(data.shape[0]/2):,:]/2.0
    data=pd.DataFrame(data)
    data.columns = header
    data=pd.concat([coordata,data],axis=1)
    data.to_csv(wfile, float_format='%.4e')
    wfile.close()
elif sys.argv[1].endswith('.vcf'):
    #write header
    #write data
    data1=np.asarray(h5py.File(sys.argv[2],'r')['pred'])
    data2=np.asarray(h5py.File(sys.argv[3],'r')['pred'])
    #compute relative diffrence and absolute difference of chromatin feature prediciton between reference and alternative alleles
    data=np.hstack([np.log2(data2/(1-data2+1e-12))-np.log2(data1/(1-data1+1e-12)),data2-data1])
    data=data[:(data.shape[0]/2),:]/2.0+data[(data.shape[0]/2):,:]/2.0
    
    SORTIND=np.arange(data.shape[0])   
    #write reference allele prediction, alternative allele prediction, relative difference and absolution difference files
    data1=data1[:(data1.shape[0]/2),:]/2.0+data1[(data1.shape[0]/2):,:]/2.0
    data2=data2[:(data2.shape[0]/2),:]/2.0+data2[(data2.shape[0]/2):,:]/2.0
    temp = pd.DataFrame(data1)
    temp.columns = header
    data1=pd.concat([coordata,temp],axis=1)
    data1=data1.iloc[SORTIND,:]
    data1.to_csv(wfile1, float_format='%.4e')
    temp = pd.DataFrame(data2)
    temp.columns = header
    data2=pd.concat([coordata,temp],axis=1)
    data2=data2.iloc[SORTIND,:]
    data2.to_csv(wfile2, float_format='%.4e')

    temp=pd.DataFrame(data[:,:len(header)])
    temp.columns = header
    data3=pd.concat([coordata,temp],axis=1)
    data3=data3.iloc[SORTIND,:]
    data3.to_csv(wfile3, float_format='%.4e')
    temp=pd.DataFrame(data[:,len(header):])
    temp.columns = header
    data4=pd.concat([coordata,temp],axis=1)
    data4=data4.iloc[SORTIND,:]
    data4.to_csv(wfile4, float_format='%.4e')

    wfile1.close()
    wfile2.close()
    wfile3.close()
    wfile4.close()  



    # #compute E-values for chromatin effects
    # ecdfs=joblib.load('./resources/ecdfs/ecdf_difflogdiff.rescaled.pkl')
    # datae=np.ones((data.shape[0],len(header)))
    # for i in range(len(header)):
    #     datae[:,i]=1-ecdfs[i](np.abs(data[:,i+len(header)]*data[:,i]))
    # datae[datae==0]=1e-6

    # if evoinfo_exist:
    #     #compute E-values for evolutionary features and Functional Significance scores
    #     check_call(['sh retrieve_gerpneutral.production.sh '+sys.argv[1]+'.wt1100.fasta.ref.vcf' ],shell=True) #create .evo1 .evo2 .evo3
    #     try:
    #         check_call(['python evoevalues.production.py '+sys.argv[1]+'.wt1100.fasta.ref.vcf' ],shell=True) #create .evo.evalues
    #         dataevoe=pd.read_csv(sys.argv[1]+'.wt1100.fasta.ref.vcf.evo.evalues',delimiter=',',header=None)
    #         dataevoe[0] = 'chr' + dataevoe[0].astype(str)
    #         matchedinds=pd.match(np.asarray(coordata['chr'].astype(str)+coordata['pos'].astype(str)),np.asarray(dataevoe[0].astype(str)+dataevoe[1].astype(str)))
    #         dataevoe = np.asarray(dataevoe.iloc[:,-4:])
    #         dataevoe = dataevoe[matchedinds,:]
    #         #impute evolutionary feature E-values for rare cases that evolutionary features are not available
    #         dataevoe[matchedinds==-1,:]=np.asarray([1,1,1,1])[np.newaxis,:]
    #         datadeepsea=np.exp(np.mean(np.log(datae),axis=1)+np.mean(np.log(dataevoe),axis=1))
    #     except:
    #         datadeepsea=np.exp(np.mean(np.log(datae),axis=1))


    #     temp=pd.DataFrame(datadeepsea[:,np.newaxis])
    #     temp.columns = ['Functional significance score']
    #     datadeepsea=pd.concat([coordata,temp],axis=1)

    #     SORTIND = np.argsort(np.asarray(datadeepsea.iloc[:,-1]))
    #     datadeepsea=datadeepsea.iloc[SORTIND,:]
    #     datadeepsea.to_csv(wfile8, float_format='%.4e')

    #     wfile7.write('Number of variants with functional significance score<0.05: '+str(np.sum(datadeepsea.iloc[:,-1]<0.05))+'\n')
    #     wfile7.write('Number of variants with functional significance score<0.01: '+str(np.sum(datadeepsea.iloc[:,-1]<0.01))+'\n')


    #     #compute probability output from HGMD/eQTL/GWAS variant classifiers
    #     inds=np.arange(data.shape[0])
    #     try:
    #         dataevo = pd.read_csv(sys.argv[1]+'.wt1100.fasta.ref.vcf.evoall',delimiter=',',header=None)
    #         dataevo[0] = 'chr' + dataevo[0].astype(str)
    #         matchedinds=pd.match(np.asarray(coordata['chr'].astype(str)+coordata['pos'].astype(str)),np.asarray(dataevo[0].astype(str)+dataevo[1].astype(str)))
    #         dataevo = np.asarray(dataevo.iloc[:,-4:])
    #         dataevo = dataevo[matchedinds,:]
    #         #impute evolutionary feature values for rare cases that evolutionary features are not available
    #         dataevo[matchedinds==-1,:]=np.asarray([0.115,-0.033,1.909,-0.200])[np.newaxis,:]
    #     except:
    #         dataevo = np.tile(np.asarray([0.115,-0.033,1.909,-0.200])[np.newaxis,:],[data.shape[0],1])

    #     model=np.loadtxt('./logistic/final_eqtl.model.txt')
    #     scalermean=np.loadtxt('./logistic/final_eqtl_scaler.mean.txt')
    #     scalerstd=np.loadtxt('./logistic/final_eqtl_scaler.std.txt')
    #     temp=pd.DataFrame(1/(1+np.exp(-(((np.abs(np.hstack([data[inds,:],dataevo[inds,:]]))-scalermean[np.newaxis,:])/scalerstd[np.newaxis,:]).dot(model[1:,np.newaxis])+model[0]))))
    #     temp.columns=['eQTL-probability']
    #     datapred=pd.concat([coordata.iloc[inds,:].reset_index(drop=True),temp.reset_index(drop=True)],axis=1)

    #     model=np.loadtxt('./logistic/final_gwas.model.txt')
    #     scalermean=np.loadtxt('./logistic/final_gwas_scaler.mean.txt')
    #     scalerstd=np.loadtxt('./logistic/final_gwas_scaler.std.txt')
    #     temp=pd.DataFrame(1/(1+np.exp(-(((np.abs(np.hstack([data[inds,:],dataevo[inds,:]]))-scalermean[np.newaxis,:])/scalerstd[np.newaxis,:]).dot(model[1:,np.newaxis])+model[0]))))
    #     temp.columns=['GWAS-probability']
    #     datapred=pd.concat([datapred.reset_index(drop=True),temp.reset_index(drop=True)],axis=1)

    #     model=np.loadtxt('./logistic/final_hgmd.model.txt')
    #     scalermean=np.loadtxt('./logistic/final_hgmd_scaler.mean.txt')
    #     scalerstd=np.loadtxt('./logistic/final_hgmd_scaler.std.txt')
    #     temp=pd.DataFrame(1/(1+np.exp(-(((np.abs(np.hstack([data[inds,:],dataevo[inds,:]]))-scalermean[np.newaxis,:])/scalerstd[np.newaxis,:]).dot(model[1:,np.newaxis])+model[0]))))
    #     temp.columns=['HGMD-probability']
    #     datapred=pd.concat([datapred.reset_index(drop=True),temp.reset_index(drop=True)],axis=1)
    #     datapred.to_csv(wfile5, float_format='%.4e')
    # else:
    #     SORTIND=np.arange(data.shape[0])


    # #write E-values for chromatin effects
    # temp=pd.DataFrame(datae[:,:len(header)])
    # temp.columns = header
    # datae=pd.concat([coordata,temp],axis=1)
    # datae=datae.iloc[SORTIND,:]
    # datae.to_csv(wfile6,float_format='%.4e')
    # #write reference allele prediction, alternative allele prediction, relative difference and absolution difference files
    # data1=data1[:(data1.shape[0]/2),:]/2.0+data1[(data1.shape[0]/2):,:]/2.0
    # data2=data2[:(data2.shape[0]/2),:]/2.0+data2[(data2.shape[0]/2):,:]/2.0
    # temp = pd.DataFrame(data1)
    # temp.columns = header
    # data1=pd.concat([coordata,temp],axis=1)
    # data1=data1.iloc[SORTIND,:]
    # data1.to_csv(wfile1, float_format='%.4e')
    # temp = pd.DataFrame(data2)
    # temp.columns = header
    # data2=pd.concat([coordata,temp],axis=1)
    # data2=data2.iloc[SORTIND,:]
    # data2.to_csv(wfile2, float_format='%.4e')

    # temp=pd.DataFrame(data[:,:len(header)])
    # temp.columns = header
    # data3=pd.concat([coordata,temp],axis=1)
    # data3=data3.iloc[SORTIND,:]
    # data3.to_csv(wfile3, float_format='%.4e')
    # temp=pd.DataFrame(data[:,len(header):])
    # temp.columns = header
    # data4=pd.concat([coordata,temp],axis=1)
    # data4=data4.iloc[SORTIND,:]
    # data4.to_csv(wfile4, float_format='%.4e')

    # wfile1.close()
    # wfile2.close()
    # wfile3.close()
    # wfile4.close()
    # wfile6.close()
    # if evoinfo_exist:
    #     wfile5.close()
    #     wfile7.close()
    #     wfile8.close()



