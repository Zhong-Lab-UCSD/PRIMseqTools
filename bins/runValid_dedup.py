import sys
from collections import defaultdict
import glob


#read in refseq dic
dicIdGeneName={}
dicIdGeneType={}
with open('%s'%(sys.argv[3]),'r') as f:
    for line in f:
        splitLine=line.strip().split(',')
        dicIdGeneName[splitLine[0]]=splitLine[1]
        dicIdGeneType[splitLine[0]]=splitLine[2]


#a read is dsRNA if its maped strdness diff from its gene strdness
#a read is fRNA if its mapped strdness the same as its gene strdness
validIdList=set()
targetFile=open('%s/%schimericReadPairs.csv'%(sys.argv[1],sys.argv[2]),'w')
targetFile.write('readId,R1Tx,R1start,R1end,R1Gene,R1Cigar,R1Strand,R2Tx,R2start,R2end,R2Gene,R2Cigar,R2Strand,R1geneType,R2geneType\n')
dicMapInfo_count=defaultdict(int)
fileList=glob.glob('%s/%sintermediateFiles/chimericReadPairs_all_bwa.csv_*'%(sys.argv[1],sys.argv[2]))
for file in fileList:
    with open(file,'r') as f:
        next(f)
        for line in f:
            splitLine=line.strip().split(',')
            readId=splitLine[0]
            tx1,tx2=splitLine[1],splitLine[7]
            [tx1,start1,end1]=splitLine[1:4]
            [tx2,start2,end2]=splitLine[7:10]
            gene1,gene2=splitLine[4],splitLine[10]
            strand1,strand2=splitLine[6],splitLine[12]
            if tx1>tx2:
                mapInfo=','.join([tx1,start1,end1,strand1,tx2,start2,end2,strand2])
            else:
                mapInfo=','.join([tx2,start2,strand2,tx1,start1,end1,strand1])
            tp1,tp2=dicIdGeneType[tx1],dicIdGeneType[tx2]
            splitLine.append(tp1)
            splitLine.append(tp2)
            #read 1 is dsRNA
            if strand1!=strand2 and strand1=='-' and tp1 == 'protein_coding' and dicMapInfo_count[mapInfo]==0:
                infoLine=','.join(splitLine)
                targetFile.write(infoLine)
                targetFile.write('\n')
                validIdList.add(readId)
            #read 2 is dsRNA
            if strand1!=strand2 and strand2=='-' and tp2=='protein_coding'and dicMapInfo_count[mapInfo]==0:
                infoLine=','.join(splitLine)
                targetFile.write(infoLine)
                targetFile.write('\n')
                validIdList.add(readId)
            dicMapInfo_count[mapInfo]+=1
targetFile.close()

chimNum=len(validIdList)
mapSum=0
fileList=glob.glob('%s/%sintermediateFiles/mappedStats_*.txt'%(sys.argv[1],sys.argv[2]))
for file in fileList:
    with open(file,'r') as f:
        for line in f:
            splitLine=line.strip().split(',')
            mapSum+=int(splitLine[0])

targetFile=open('%s/%ssummary.csv'%(sys.argv[1],sys.argv[2]),'a')
targetFile.write('#mapped_read_pairs,%d\n'%(mapSum))
targetFile.write('#chimeric_read_pairs,%d\n'%(chimNum))
targetFile.close()
