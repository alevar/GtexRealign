#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import sys
import subprocess
from pybedtools import BedTool
import re
import argparse

#first need to extract splice junctions
def extractFlagBits(data):
    data["paired"]=data["FLAG"]               &1 #template having multiple segments in sequencing
    data["aligned2Mates"]=data["FLAG"]        &2 #each segment properly aligned according to the aligner
    data["unmappedCurr"]=data["FLAG"]         &4 #segment unmapped
    data["unmappedMate"]=data["FLAG"]         &8 #next segment in the template unmapped
    data["reversedCurr"]=data["FLAG"]         &16 #SEQ being reverse complemented
    data["reversedMate"]=data["FLAG"]         &32 #SEQ of the next segment in the template being reverse complemented
    data["firstRead"]=data["FLAG"]            &64 #the first segment in the template
    data["lastRead"]=data["FLAG"]             &128 #the last segment in the template
    data["secondaryAlignment"]=data["FLAG"]   &256 #secondary alignment
    data["noPassFilter"]=data["FLAG"]         &512 #not passing filters, such as platform/vendor quality controls
    data["PCRdup"]=data["FLAG"]               &1024 #PCR or optical duplicate
    data["suppAl"]=data["FLAG"]               &2048 #supplementary alignment

def se(row):
    cigar=row["CIGAR"]
    chars=re.findall(r"[\D']+", cigar)
    ints=[int(x) for x in re.findall(r"[\d']+",cigar)]
    readLen=0
    pre=0
    post=0
    n=0
    m_pre_tem=0
    m_pre_ref=0
    m_post_tem=0
    m_post_ref=0
    indexN=0
    di=0
    blockCount=1
    blockSizes=[0]
    tStarts=[0]
    if "N" in chars:
        indexN=chars.index("N")
        blockCount=len(cigar.split("N"))
    for i in range(len(chars)):
        if i==0 and chars[i] in "SH":
            pre=ints[i]
            readLen=readLen+ints[i]
#             tStarts[0]=tStarts[0]+ints[i]
        if i==len(chars)-1 and chars[i] in "SH":
            post=ints[i]
            readLen=readLen+ints[i]
        if chars[i]=="N":
            n=n+ints[i]
            tStarts.append(tStarts[-1]+blockSizes[-1]+ints[i])
            blockSizes.append(0)
        if i<indexN and chars[i]=="M":
            m_pre_tem=m_pre_tem+ints[i]
            m_pre_ref=m_pre_ref+ints[i]
            readLen=readLen+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
        if i>=indexN and chars[i]=="M":
            m_post_tem=m_post_tem+ints[i]
            m_post_ref=m_post_ref+ints[i]
            readLen=readLen+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
        if i<indexN and chars[i]=="D":
            m_pre_ref=m_pre_ref+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
            di=di+ints[i]
        if i>=indexN and chars[i]=="D":
            m_post_ref=m_post_ref+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
            di=di+ints[i]
        if i<indexN and chars[i]=="I":
            readLen=readLen+ints[i]
            m_pre_tem=m_pre_tem+ints[i]
            di=di+ints[i]
        if i>=indexN and chars[i]=="I":
            readLen=readLen+ints[i]
            m_post_tem=m_post_tem+ints[i]
            di=di+ints[i]
    return pd.Series([pre,post,m_pre_ref,m_pre_tem,m_post_ref,m_post_tem,n,readLen,di,blockCount,",".join([str(x) for x in blockSizes]),",".join([str(x) for x in tStarts])])

def parseCIGAR(data):
    data["CIGAR"].replace("*",np.nan,inplace=True)
    data.dropna(axis=0,inplace=True)
    data.reset_index(drop=True,inplace=True)

#     data["READ_LEN"]=data.SEQ.str.len()
    data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$",expand=False).replace(np.nan,0).astype(int)
    data["END"]=data.READ_LEN-data.CIGAR_POST
    data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[SH]",expand=False).replace(np.nan,0).astype(int)

    data16=data[data["reversedCurr"]==16].reset_index(drop=True)
    data0=data[data["reversedCurr"]==0].reset_index(drop=True)
    data16["Template_start"]=data16.READ_LEN-data16.END
    data16["Template_end"]=data16.READ_LEN-data16.CIGAR_PRE
    data0["Template_start"]=data0.CIGAR_PRE
    data0["Template_end"]=data0.END

    data16["Reference_start"]=data16.READ_LEN-data16.END+data16.POS-data16.Template_start
    data16["Reference_end"]=data16.READ_LEN-data16.CIGAR_PRE-1+data16.POS-data16.Template_start+data16.N
    data0["Reference_start"]=data0.POS
    data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE+data0.N 
    
    data=pd.concat([data16,data0]).reset_index(drop=True)
    data.drop(["CIGAR_POST","CIGAR_PRE"],axis=1,inplace=True)
    return data

def tStarts(row):
    re=row["Reference_start"]
    tStarts=[]
    qStarts=row.qStarts.split(",")
    qs2=[int(x)-int(qStarts[0]) for x in qStarts]
    blockSizes=row.blockSizes.split(",")
    for i in range(row["blockCount"]):
        tStarts.append(re+int(qStarts[i]))
    return pd.Series([",".join([str(x) for x in tStarts])])

def getSS(row):
    ns=row["N"]
    bs=[int(x) for x in row["blockSizes"].split(",")]
    ts=[int(x) for x in row["tStarts"].split(",")]
    res=[]
    for i in range(len(bs)-1):
        res.append(str(ts[i]+bs[i]-1)+"-"+str(ts[i+1]-1))
    return ",".join(res)

def findBlocks(row): # returns the minimum block size for a given junction
    sjs=row["fullSS"].split(",") # list of spliceJunctions
    cursj=row["ss"]
    fbs=row["blockSizes"].split(",") # list of blocks
    fbIDX=sjs.index(cursj) # index of the first junction
    curBlocks=[int(x) for x in fbs[fbIDX:fbIDX+2]]
    return min(curBlocks)

# main part of the tool
def sam2splice(args):
    samCols=['QNAME',
            'FLAG',
            'RNAME',
            'POS',
            'MAPQ',
            'CIGAR',
            'RNEXT',
            'PNEXT',
            'TLEN',
            'SEQ',
            'QUAL',
            'strand']
    df=pd.DataFrame([],columns=samCols)

    with open(os.path.abspath(args.input),"r") as inFP:
        for line in inFP.readlines():
            if line[0]=="#" or line[0]=="@":
                continue
            lineCols=line.split("\t")

            if "N" in lineCols[5]: # cigar string contains a splice junctions
                # there has to be an XS:A: flag
                strandCols=line.split("XS:A:")
                if not len(strandCols)==2:
                    continue
                strand=strandCols[1].split("\t")[0]
                df=df.append(pd.DataFrame([lineCols[:11]+[strand]],columns=samCols))

    inFP.close()

    df["FLAG"]=df.FLAG.astype(int)
    df["POS"]=df["POS"].astype(int)
    df["MAPQ"]=df["MAPQ"].astype(int)
    df["PNEXT"]=df["PNEXT"].astype(int)
    df["TLEN"]=df["TLEN"].astype(int)

    df=df[df["CIGAR"].str.contains("N")].reset_index(drop=True)

    extractFlagBits(df)

    df["PRE"]=np.nan
    df["POST"]=np.nan
    df["MPRER"]=np.nan
    df["MPRET"]=np.nan
    df["MPOSTR"]=np.nan
    df["MPOSTT"]=np.nan
    df["N"]=np.nan
    df["blockCount"]=0
    df["blockSizes"]=""
    df["qStarts"]=""
    df["tStarts"]=""
    df[["PRE","POST","MPRER","MPRET","MPOSTR","MPOSTT","N","READ_LEN","DI","blockCount","blockSizes","qStarts"]]=pd.DataFrame(df.apply(lambda row: se(row),axis=1))
    df=parseCIGAR(df)
    df["tStarts"]=df.apply(lambda row: tStarts(row),axis=1)
    df.drop(["FLAG"],axis=1,inplace=True)

    df["ss"]=df.apply(lambda row: getSS(row),axis=1)
    df.reset_index(drop=True,inplace=True)

    df["dup"]=df.duplicated("QNAME").astype(int)
    df["QNAME"]=df["QNAME"]+"_"+df["dup"].astype(str)

    tdf=pd.concat([pd.Series(row['QNAME'], row['ss'].split(',')) for _, row in df.iterrows()]).reset_index()
    tdf.columns=["ss","QNAME"]

    tdf=tdf.merge(df[["QNAME","RNAME","reversedCurr","MAPQ","SEQ","CIGAR","blockSizes","ss","strand"]],on="QNAME",how="left")

    tdf["type"]="intron"
    tdf["source"]="hisat"
    tdf["start"]=tdf["ss_x"].str.split("-",expand=True)[0]
    tdf["end"]=tdf["ss_x"].str.split("-",expand=True)[1].astype(int)+1
    tdf["end"]=tdf["end"].astype(str)
    assert len(tdf[tdf["start"]<tdf["end"]])
    tdf["x"]="."
    tdf["y"]="."
    tdf=tdf[["RNAME","source","type","start","end","x","strand","y"]]
    tdf=pd.DataFrame(tdf.groupby(by=["RNAME","source","type","start","end","x","strand","y"]).apply(len)).reset_index()
    tdf["attributes"]="weight="+tdf[0].astype(str)
    tdf[["RNAME","source","type","start","end","x","strand","y","attributes"]].to_csv(args.output,sep="\t",index=False,header=False)


def main(argv):
    parser=argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('-i',
                        '--input',
                        required=True,
                        type=str,
                        help="splice alignment")
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help="output file")
    parser.add_argument('-f',
                        '--filter',
                        action="store_true",
                        help="filter LTR splice junctions")
    
    parser.set_defaults(func=sam2splice)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])
