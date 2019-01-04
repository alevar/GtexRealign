#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import re

def mergeGFF(args):
	gff3Cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
	chrMapCols=["name","role","molecule","type","genbank","rel","refseq","unit","seqLen","ucsc"]

	# Here we shall develop a method for merging two GFF files
	# Main idea - need to convert any arbitrary features into gene-transcript-exon features and discard anything unnecessary (such as CDS)
	#    - first we shall try to guess parent-child hierarchy
	#    - based on the hierarchy we can decide which entries to keep
	#    - alternatively, we can rely on exons, and assign "transcript" to their features and "gene" to the respective grandparents

	# load the assembly stats to be used in chromosome name conversion
	chrMap=pd.read_csv(args.stats,sep="\t",names=chrMapCols,comment="#")

	df=pd.read_csv(args.ann1,sep="\t",skiprows=7,names=gff3Cols)

	# df=df[df["seqid"]=="chr11"].reset_index(drop=True)

	df["id"]=df.attributes.str.split("ID=",expand=True)[1].str.split(";",expand=True)[0]
	dfE=df[df["type"]=="exon"].reset_index(drop=True)
	dfE["transID"]=dfE.attributes.str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
	transcriptSet=set(dfE["transID"])
	dfT=df[df["id"].isin(transcriptSet)].reset_index(drop=True)
	assert len(set(dfT["id"]))==len(transcriptSet)
	dfT["geneID"]=dfT.attributes.str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
	dfT["transID"]=dfT["id"]

	# here's the time to deal with inconsistencies
	# first remove any enties that do not have a parent - those do not conform to the expected structure that we define
	# however, this way we remove PSEUDOGENES, which may not be the desired effect - need to think more/discuss
	dfT=dfT[~(dfT["geneID"].isnull())].reset_index(drop=True)
	dfE=dfE[dfE["transID"].isin(dfT["transID"])].reset_index(drop=True)

	# next it is time to iteratively, bring the parents to down
	# for instance, a feature that is defined as a transcript at this level (eg. miRNA) may have aother feature at this level as a parent
	while True:
	    dfI=dfT[["transID","geneID"]].merge(dfT[["transID","geneID"]],how="inner",left_on="geneID",right_on="transID")
	    dfI=dfI[["transID_x","geneID_y"]]
	    dfT=dfT.merge(dfI,how="left",left_on="transID",right_on="transID_x")
	    dfT["geneID"]=np.where(dfT["geneID_y"].isnull(),dfT["geneID"],dfT["geneID_y"])
	    if len(set(dfI["geneID_y"]).intersection(set(dfT["transID"])))==0:
	        break

	dfT["attributes"]=dfT["attributes"]+";old_type="+dfT["type"]
	dfT["type"]="transcript"
	geneSet=set(dfT["geneID"])
	dfG=df[df["id"].isin(geneSet)].reset_index(drop=True)
	dfG["attributes"]=dfG["attributes"]+";old_type="+dfG["type"]
	dfG["type"]="gene"
	assert len(set(dfG["id"]))==len(geneSet)
	dfG["geneID"]=dfG["id"]
	dfG["transID"]=np.inf

	dfE=dfE.merge(dfT[["transID","geneID"]],how="left",left_on="transID",right_on="transID")
	assert len(dfE[dfE["geneID"].isnull()])==0,"incompatible features"

	# now need to concatenate all these three dataframes
	dfM=pd.concat([dfG[["seqid","source","type","start","end","score","strand","phase","attributes","transID","geneID"]],dfT[["seqid","source","type","start","end","score","strand","phase","attributes","transID","geneID"]],dfE[["seqid","source","type","start","end","score","strand","phase","attributes","transID","geneID"]]],axis=0).reset_index(drop=True)
	dfM["type"]=pd.Categorical(dfM["type"],categories=["gene","transcript","exon"],ordered=True)
	dfM=dfM.sort_values(by=["geneID","transID","type"]).reset_index(drop=True)
	del dfE
	del dfT
	del dfG
	# now need to convert the chromosome names to a standard (genbank)
	chrCurMap={}
	setNames=set(dfM["seqid"])
	for i in setNames:
	    for j in ["name","molecule","genbank","refseq","ucsc"]:
	        if len(chrMap[chrMap[j]==i])>0:
	            chrCurMap[i]=chrMap[chrMap[j]==i]["genbank"].iloc[0]
	dfM.replace(chrCurMap,inplace=True)
	dfM["start"]=dfM["start"].astype(int)
	dfM["end"]=dfM["end"].astype(int)
	# some pseudogenes are still present

	dfM_g=dfM[dfM["type"].isin(["exon"])][["seqid","strand","type","start","end","transID"]].reset_index(drop=True)
	dfM_g=dfM_g.sort_values(by=["transID","start"]).reset_index(drop=True)
	dfM_g["ns"]=dfM_g.start.astype(str).shift(-1)
	dfM_g["nid"]=dfM_g.transID.shift(-1)
	dfM_g["ns"]=np.where(dfM_g["transID"]==dfM_g["nid"],dfM_g["ns"],"***")
	dfM_g["ns"]=dfM_g.ns.astype(str)+"+"
	dfM_g["ends"]=dfM_g.end.astype(str)+"+"
	dfM_g=dfM_g.groupby(by="transID").agg({'start':'min','ns':'sum','ends':'sum','end':'max','strand':'min','seqid':'min'})
	dfM_g.reset_index(inplace=True)
	dfM_g.columns=['transID','start_min',"ns",'end','end_max','strand','seqid']
	dfM_g["ns"]=dfM_g["ns"].str.rstrip("***+")
	dfM_g["end"]=dfM_g["end"].str.rstrip("***+")
	dfM_g["chain"]=dfM_g['seqid']+":"+dfM_g['strand']+'@'+dfM_g['end']+'-'+dfM_g['ns']
	dfM_g=dfM_g[['transID','start_min','end_max','chain']]

	del setNames
	del chrCurMap
	del geneSet
	del transcriptSet
	dfM=dfM[gff3Cols+["transID"]]

	# now to remove duplicate intron chains
	dfM_g.drop_duplicates("chain",keep='first',inplace=True)

	df=pd.read_csv(args.ann2,sep="\t",skiprows=7,names=gff3Cols)

	df["id"]=df.attributes.str.split("ID=",expand=True)[1].str.split(";",expand=True)[0]
	dfE=df[df["type"]=="exon"].reset_index(drop=True)
	dfE["transID"]=dfE.attributes.str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
	transcriptSet=set(dfE["transID"])
	dfT=df[df["id"].isin(transcriptSet)].reset_index(drop=True)
	assert len(set(dfT["id"]))==len(transcriptSet)
	dfT["geneID"]=dfT.attributes.str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
	dfT["transID"]=dfT["id"]

	# here's the time to deal with inconsistencies
	# first remove any enties that do not have a parent - those do not conform to the expected structure that we define
	# however, this way we remove PSEUDOGENES, which may not be the desired effect - need to think more/discuss
	dfT=dfT[~(dfT["geneID"].isnull())].reset_index(drop=True)
	dfE=dfE[dfE["transID"].isin(dfT["transID"])].reset_index(drop=True)

	# next it is time to iteratively, bring the parents to down
	# for instance, a feature that is defined as a transcript at this level (eg. miRNA) may have aother feature at this level as a parent
	while True:
	    dfI=dfT[["transID","geneID"]].merge(dfT[["transID","geneID"]],how="inner",left_on="geneID",right_on="transID")
	    dfI=dfI[["transID_x","geneID_y"]]
	    dfT=dfT.merge(dfI,how="left",left_on="transID",right_on="transID_x")
	    dfT["geneID"]=np.where(dfT["geneID_y"].isnull(),dfT["geneID"],dfT["geneID_y"])
	    if len(set(dfI["geneID_y"]).intersection(set(dfT["transID"])))==0:
	        break

	dfT["attributes"]=dfT["attributes"]+";old_type="+dfT["type"]
	dfT["type"]="transcript"
	geneSet=set(dfT["geneID"])
	dfG=df[df["id"].isin(geneSet)].reset_index(drop=True)
	del df
	dfG["attributes"]=dfG["attributes"]+";old_type="+dfG["type"]
	dfG["type"]="gene"
	assert len(set(dfG["id"]))==len(geneSet)
	dfG["geneID"]=dfG["id"]
	dfG["transID"]=np.inf

	dfE=dfE.merge(dfT[["transID","geneID"]],how="left",left_on="transID",right_on="transID")
	assert len(dfE[dfE["geneID"].isnull()])==0,"incompatible features"

	# now need to concatenate all these three dataframes
	dfM_2=pd.concat([dfG[["seqid","source","type","start","end","score","strand","phase","attributes","transID","geneID"]],dfT[["seqid","source","type","start","end","score","strand","phase","attributes","transID","geneID"]],dfE[["seqid","source","type","start","end","score","strand","phase","attributes","transID","geneID"]]],axis=0).reset_index(drop=True)
	dfM_2["type"]=pd.Categorical(dfM_2["type"],categories=["gene","transcript","exon"],ordered=True)
	dfM_2=dfM_2.sort_values(by=["geneID","transID","type"]).reset_index(drop=True)
	del dfT
	del dfG
	# now need to convert the chromosome names to a standard (genbank)
	chrCurMap={}
	setNames=set(dfM_2["seqid"])
	for i in setNames:
	    for j in ["name","molecule","genbank","refseq","ucsc"]:
	        if len(chrMap[chrMap[j]==i])>0:
	            chrCurMap[i]=chrMap[chrMap[j]==i]["genbank"].iloc[0]
	dfM_2.replace(chrCurMap,inplace=True)
	dfM_2["start"]=dfM_2["start"].astype(int)
	dfM_2["end"]=dfM_2["end"].astype(int)
	# some pseudogenes are still present

	dfM_2_g=dfM_2[dfM_2["type"].isin(["exon"])][["seqid","strand","type","start","end","transID"]].reset_index(drop=True)
	dfM_2_g=dfM_2_g.sort_values(by=["transID","start"]).reset_index(drop=True)
	dfM_2_g["ns"]=dfM_2_g.start.astype(str).shift(-1)
	dfM_2_g["nid"]=dfM_2_g.transID.shift(-1)
	dfM_2_g["ns"]=np.where(dfM_2_g["transID"]==dfM_2_g["nid"],dfM_2_g["ns"],"***")
	dfM_2_g["ns"]=dfM_2_g.ns.astype(str)+"+"
	dfM_2_g["ends"]=dfM_2_g.end.astype(str)+"+"
	dfM_2_g=dfM_2_g.groupby(by="transID").agg({'start':'min','ns':'sum','ends':'sum','end':'max','strand':'min','seqid':'min'})
	dfM_2_g.reset_index(inplace=True)
	dfM_2_g.columns=['transID','start_min',"ns",'end','end_max','strand','seqid']
	dfM_2_g["ns"]=dfM_2_g["ns"].str.rstrip("***+")
	dfM_2_g["end"]=dfM_2_g["end"].str.rstrip("***+")
	dfM_2_g["chain"]=dfM_2_g['seqid']+":"+dfM_2_g['strand']+'@'+dfM_2_g['end']+'-'+dfM_2_g['ns']
	dfM_2_g=dfM_2_g[['transID','start_min','end_max','chain']]

	del setNames
	del chrCurMap
	del geneSet
	del transcriptSet
	dfM_2=dfM_2[gff3Cols+["transID"]]

	# now to remove duplicate intron chains
	dfM_2_g.drop_duplicates("chain",keep='first',inplace=True)

	# now we can merge two annotations

	# in doing so, we need to first identify all which overlap completely (that is having the sme intron chain,start and end)
	shared_chain=dfM_g.merge(dfM_2_g,how="inner",on=["chain"])#[["transID_x",'transID_y']]
	del dfM_g
	del dfM_2_g
	shared_all=shared_chain[(shared_chain['start_min_x']==shared_chain['start_min_y']) & (shared_chain['end_max_x']==shared_chain['end_max_y'])]
	# now we need to deal with those that share intron chain but have different start end
	# in these cases, the start or end or both exons need to be updated as well as the start or end or both of the transcripts
	sc=shared_chain[shared_chain["start_min_x"]>shared_chain["start_min_y"]]
	dfM=dfM.merge(sc,how="left",left_on="transID",right_on="transID_x")
	dfM['start']=np.where(dfM['start']==dfM['start_min_x'],dfM['start_min_y'],dfM['start'])
	dfM=dfM[gff3Cols+["transID"]]
	# second case
	sc=shared_chain[shared_chain["end_max_x"]<shared_chain["end_max_y"]]
	dfM=dfM.merge(sc,how="left",left_on="transID",right_on="transID_x")

	dfM['end']=np.where(dfM['end']==dfM['end_max_x'],dfM['end_max_y'],dfM['end'])
	# now we only need to merge the rest of transcripts

	# right now there is no need to create new gene entries and standardize the parents
	# this is unnecessary since, we are using this annotation only to create a transcriptome fasta reference file
	# which by definition will only extract transcripts, and we only need to make sure there are no intron-chin duplicates
	dfM_2=dfM_2[~(dfM_2['transID'].isin(set(shared_chain["transID_x"])))].reset_index(drop=True)
	# lastly, we shall build the final GFF by merging two annotations
	dfF=pd.concat([dfM[gff3Cols],dfM_2[gff3Cols]],axis=0)
	dfF.to_csv(args.outFP,sep='\t',index=False,header=False)


def main(argv):
    parser = argparse.ArgumentParser(description='''Help Page\nMerge two annotations (GFF/GTF) without redundancy in transcriptome''')

#===========================================
#===================BUILD===================
#===========================================
    parser.add_argument('--ann1',
                        required=True,
                        type=str,
                        help="first annotation of the genome in gff/gtf format")
    parser.add_argument('--ann2',
                        required=False,
                        type=str,
                        help="second annotation of the genome in gff/gtf format")
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help="output file")
    parser.add_argument('--stats',
                        required=True,
                        type=str,
                        help="assembly report for the release of the genome assembly")

    parser.set_defaults(func=mergeGFF)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])