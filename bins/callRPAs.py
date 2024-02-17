from collections import defaultdict
import scipy
import glob
import scipy.stats as stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import sys

targetFile=open('%s/%sRNAProteinInteractions.csv'%(sys.argv[1],sys.argv[5]),'a')

def getIntCount(filePath):
    dicIntCount_positive=defaultdict(int)
    dicProteinCount_positive=defaultdict(int)
    dicRNACount_positive=defaultdict(int)
    fileList=glob.glob(filePath)
    for filePath in fileList:
        with open(filePath,'r') as f:
            next(f)
            for line in f:
                splitLine=line.strip().split(',')
                gene1,gene2,strand1,strand2=splitLine[4],splitLine[10],splitLine[6],splitLine[12]
                if strand1=='+':
                    rnaGene,proteinGene=gene1,gene2
                else:
                    rnaGene,proteinGene=gene2,gene1
                pair=';'.join([rnaGene,proteinGene])
                dicProteinCount_positive[proteinGene]+=1
                dicRNACount_positive[rnaGene]+=1
                dicIntCount_positive[pair]+=1        
    #sort
    sorted_x_positive = sorted(dicIntCount_positive.items(), key=lambda kv: kv[1],reverse=True)
    return dicIntCount_positive,dicProteinCount_positive,dicRNACount_positive,sorted_x_positive

def identifyRPAs_chimericAdj(
    sorted_x_positive,dicIntCount_positive,dicProteinCount_positive,dicRNACount_positive,coEff,pCutOff,orCutoff):
    factor=sum([x[1] for x in sorted_x_positive])/len(sorted_x_positive)
    chimTotal=sum([x[1] for x in sorted_x_positive])
    pvalueList=[]
    sorted_x_1_select=[]
    selectList_1=[]
    posRCList_1=[]
    orList_1=[]
    chiList=[]
    for ha in sorted_x_positive:
        [rna,protein]=ha[0].split(';')
        a=dicIntCount_positive[ha[0]]
        b=dicRNACount_positive[rna]/2-a
        c=dicProteinCount_positive[protein]/2-a
        d=chimTotal-a-b-c
        b,c=max(0,b),max(0,c)
        oddsRatio=(a+1)*(d+1)/(b+1)/(c+1)
        chi2, p, dof, ex=stats.chi2_contingency([[a+1,b+1],[c+1,d+1]])
        orList_1.append(oddsRatio)
        pvalueList.append(p)
        sorted_x_1_select.append(ha)
        selectList_1.append(ha[0])
        posRCList_1.append(a)
        chiList.append(chi2)


    stats1 = importr('stats')
    pvalueList_adj_1=stats1.p_adjust(FloatVector(pvalueList), method = 'BH')

    lolCount=0
    count=0
    list1=[]
    rcList1=[]
    pvalueSig_1=[]
    orSig_1=[]
    chiSig=[]
    for i in range(len(selectList_1)):
        ha=selectList_1[i]
        count+=1
        gene1,gene2=ha.split(';')
        pAdj=pvalueList_adj_1[i]
        rcc=posRCList_1[i]
        orr=orList_1[i]
        chichi=chiList[i]
        #targetFile.write('%s,%s,%d,%d,%f,%f\n'%(gene1,gene2,counts,counts2,cpmaFC,pAdj))
        if pAdj<=pCutOff and rcc>coEff*factor and 'MTRNR' not in gene1 and 'MTRNR' not in gene2 and orr>orCutoff:
        #if pAdj<=pCutOff and rcc>coEff*factor and orr>1:
            lolCount+=1
            list1.append(ha)
            pvalueSig_1.append(pAdj)
            orSig_1.append(orList_1[i])
            rcList1.append(rcc)
            chiSig.append(chichi)
    return list1,rcList1,orSig_1,chiSig,pvalueSig_1


dicIntCount_positive,dicProCount_positive,dicRNACount_positive,sorted_positive=getIntCount(
    '%s/%schimericReadPairs.csv'%(sys.argv[1],sys.argv[5]))


list_RPA,rcList_RPA,orList_RPA,chiList_RPA,pvalueList_RPA=identifyRPAs_chimericAdj(
    sorted_positive,dicIntCount_positive,dicProCount_positive,dicRNACount_positive,float(sys.argv[4]),float(sys.argv[2]),float(sys.argv[3]))

#write into the file
targetFile.write('RNA,Protein,ReadCount,FDR,oddsRatio,chiSquareStat\n')
for i in range(len(list_RPA)):
    ha=list_RPA[i]
    [gene1,gene2]=ha.split(';')
    oddsRatio=str(orList_RPA[i])
    chichi=str(chiList_RPA[i])
    pp=str(pvalueList_RPA[i])
    rc=str(dicIntCount_positive[ha])
    infoList=','.join([gene1,gene2,rc,pp,oddsRatio,chichi])
    targetFile.write(infoList)
    targetFile.write('\n')
    
targetFile.close()


targetFile=open('%s/%ssummary.csv'%(sys.argv[1],sys.argv[5]),'a')
targetFile.write('#RNA-protein_associations,%d\n'%(len(list_RPA)))
targetFile.close()