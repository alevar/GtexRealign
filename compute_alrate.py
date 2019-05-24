#!/usr/bin/env python
import argparse
import pysam
import sys
import os

def run_count(args):
    inSAM=pysam.AlignmentFile(os.path.abspath(args.input),"rb")
    setPaired=dict()
    setSingle=dict()
    for read in inSAM:
        if not read.is_unmapped and read.is_proper_pair and not read.mate_is_unmapped:
            setPaired.setdefault(read.qname,1)
            setPaired[read.qname]+=1
        if not read.is_unmapped and (not read.is_proper_pair or read.mate_is_unmapped):
            setSingle.setdefault(read.qname,1)
            setSingle[read.qname]+=1
    inSAM.close()

    uniquePaired=len(setPaired)*2
    uniqueSingle=len(set(setSingle).difference(set(setPaired)))
    print("unique\t"+str(uniquePaired+uniqueSingle))
    secondaryPaired=sum(setPaired.values())-uniquePaired
    secondarySingle=sum(setSingle.values())-uniqueSingle
    print("secondary\t"+str(secondaryPaired+secondarySingle))

def main(args):
    parser=argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('-i',
                        '--input',
                        required=True,
                        type=str,
                        help="alignment path")
    
    parser.set_defaults(func=run_count)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])