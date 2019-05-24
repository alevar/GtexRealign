#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Author: Ales Varabyou
"""

# this script generates a GFF3-formatted annotation of all unique introns
# from a given genome annotation provided in a GFF3/GTF format

import pandas as pd
import numpy as np
import argparse
import sys
import os

def gff2intron(args):
    gff3cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
    df=pd.read_csv(args.input,sep="\t",names=gff3cols,comment="#")
    df=df[df["type"]=="exon"].reset_index(drop=True)
    df["parent"]=df["attributes"].str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
    df.sort_values(by=["parent","seqid","strand","start","end"],ascending=True,inplace=True)
    df["next_parent"]=df["parent"].shift(-1)
    df["next_start"]=df["start"].shift(-1)
    df.dropna(inplace=True)
    df["next_start"]=df.next_start.astype(int)
    df=df[df["parent"]==df["next_parent"]].reset_index(drop=True)
    
    df["start"]=df["end"]
    df["end"]=df["next_start"]
    df["type"]="intron"
    df.drop_duplicates(["seqid","strand","start","end"],inplace=True)
    df.drop(["next_start","parent","next_parent"],axis=1,inplace=True)
    df.to_csv(args.output,sep="\t",index=False,header=False)

def main(args):
    parser=argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-i",
                       "--input",
                       required=True,
                       help="annotation in a GFF/GTF format")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="output file for the intronic annotation")
    parser.set_defaults(func=gff2intron)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])
