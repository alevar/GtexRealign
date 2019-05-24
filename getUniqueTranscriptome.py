#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys
import csv
import os

def wrapper(args):
    gff3Cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
    readLen=76 # the target read length
    # # first load the annotation
    # df=pd.read_csv(args.input,sep="\t",names=gff3Cols,comment="#")
    
    # df=df[df['type']=='exon'].reset_index(drop=True)
    # df['p']=df.attributes.str.split('Parent=',expand=True)[1].str.split(";",expand=True)[0]
    
    # df=df[['seqid','start','end','strand','p']]
    # df.reset_index(drop=True,inplace=True)
    # df['gene']="CHS."+df.p.str.split(".",expand=True)[1]
    
    # # now to actually find things of interest
    # # we need to specify intervals and find all intervals that overlap
    # df.sort_values(by=["seqid",'strand','start','end','gene','p'],ascending=True,inplace=True)
    # df.reset_index(drop=True,inplace=True)
    # df["start"]=df.start.astype(int)
    # df["end"]=df.end.astype(int)
    
    # dfg=df.groupby(by="p")
    # setIntervals=set()
    # results=[]
    # count=0
    # for name, group in dfg:
    #     count+=1
    #     if count%1000==0:
    #         print(count)
    #     assert len(set(group["seqid"].tolist()))==1,"wrong chromosomes in : "+name
    #     assert len(set(group["seqid"].tolist()))==1,"wrong strand in : "+name
    #     chrid=group.seqid.min()
    #     strand=group.strand.min()
    #     intervals=[]
    #     for index, row in group.iterrows():
    #         intervals.extend(list(range(row["start"],row["end"]+1)))
    #     curTmp=[] # current temporarry interval for extension
    #     g=""
    #     tmp=[]
    #     extending=False
    #     for i in range(len(intervals)+readLen-1):
    #         tmp=intervals[i:i+readLen]
    #         if len(tmp)==readLen:
    #             g=chrid+":"+strand+"@"+",".join([str(x[0])+"-"+str(x[-1]) for x in np.split(tmp,np.where(np.diff(tmp)!=1)[0]+1)])
    #             if not g in setIntervals:
    #                 if extending:
    #                     curTmp.append(tmp[-1])
    #                 else:
    #                     curTmp=tmp
    #                     extending=True
    #                 setIntervals.add(g)
    #             else:
    # #                 print("hi")
    #                 if extending:
    #                     for x in np.split(curTmp,np.where(np.diff(curTmp)!=1)[0]+1):
    #                         results.append((name,)+(x[0],x[-1]))
    # #                 else:
    # #                     for x in np.split(tmp,np.where(np.diff(tmp)!=1)[0]+1):
    # #                         results.append((name,)+(x[0],x[-1]))
    #                 curTmp=[]
    #                 extending=False
    #     if extending:
    #         for x in np.split(curTmp,np.where(np.diff(curTmp)!=1)[0]+1):
    #             results.append((name,)+(x[0],x[-1]))
    #     else:
    #         if len(tmp)==readLen and not g in setIntervals:
    #             for x in np.split(tmp,np.where(np.diff(tmp)!=1)[0]+1):
    #                 results.append((name,)+(x[0],x[-1]))

    # dfn=pd.DataFrame.from_records(results,columns=['p','start','end'])
    # dfn=dfn.merge(df[['p','seqid','strand']],how='left',on='p')
    # dfn["attributes"]='Parent='+dfn['p']
    # dfn["score"]="."
    # dfn["phase"]="."
    # dfn["type"]="exon"
    # dfn["source"]="test"
    # dfn.reset_index(drop=True)
    # dfn.to_csv(args.output,index=False)
    # dfn[gff3Cols].to_csv(args.output+".gff",sep="\t",index=False,quoting=csv.QUOTE_NONE,header=False)

    df=pd.read_csv(args.output+".gff",names=gff3Cols,sep='\t')
    df.drop_duplicates(inplace=True)
    df.sort_values(by=["seqid","strand","start","end"],inplace=True)
    df.reset_index(drop=True,inplace=True)
    df.reset_index(inplace=True)
    df['attributes']="Parent="+df['index'].astype(str)+";old"+df["attributes"]
    dfT=df.copy(deep=True)
    dfT["type"]="transcript"
    dfT["id"]=dfT.attributes.str.split("oldParent=",expand=True)[1]
    dfT["attributes"]="ID="+df["index"].astype(str)+";oldName="+dfT["id"]
    dfT.drop(["id"],axis=1,inplace=True)
    df=pd.concat([df,dfT],axis=0)
    df["type"]=pd.Categorical(df["type"],categories=["transcript","exon"],ordered=True)
    df.sort_values(by=["index","type"],ascending=True,inplace=True)
    df.drop(["index"],axis=1,inplace=True)
    df.reset_index(drop=True,inplace=True)
    df[gff3Cols].to_csv(args.output+".final.gtf",sep="\t",index=False,quoting=csv.QUOTE_NONE,header=False)

def setup(args):
    parser=argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-i",
                       "--input",
                       type=str,
                       required=True,
                       help="input gff")
    parser.add_argument("-o",
                       "--output",
                       type=str,
                       required=True,
                       help="output csv")
    parser.set_defaults(func=wrapper)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    setup(sys.argv[1:])
