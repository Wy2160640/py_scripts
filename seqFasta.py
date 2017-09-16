#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: WangYaheng
Date: Sep-16-2017 13:20
"""

import argparse, sys
from pyfasta import Fasta

def get_opt():
	'''Check and parsing the opts'''
	parser=argparse.ArgumentParser(
				prog='seqFasta', 
				description='seqFasta: A program to get fasta file from BED file ', 
				usage='%(prog)s.py [options] -g genome -b bedfile'
			)
	parser.add_argument('-g','--genome',nargs='?',type=str,default=sys.stdin, help="The path of fasta genome file, the suffix must be '.fa' or '.fasta'. [Required]",required=True)
	parser.add_argument('-b','--bedfile',nargs='?',type=argparse.FileType('r'),default=sys.stdin, help="BED file. detail, http://genome.ucsc.edu/FAQ/FAQformat.html#format1 .[Required]",required=True)
	parser.add_argument('-s','--seqname',nargs='?',type=str,default='query',help="The prefix of fasta seq name",required=False)
	parser.add_argument('-f','--flank',nargs='?',type=int,default=0,help="The number of the upstream and downstream slip bases.",required=False)
	parser.add_argument('-o','--outfile',nargs='?',type=argparse.FileType('w'),default=sys.stdout,help="Output file name for storing the results, default is stdout.",required=False)
	args=parser.parse_args()
	return args


def run(args):
	genome=Fasta(args.genome)
	bed=filter(lambda x: x.strip(), args.bedfile.readlines())
	bed_list=map(lambda x:x.strip().split(), bed)
	result=map(lambda i: '>{0}_{1}\n{2}'.format(args.seqname, i+1, genome.sequence({'chr':bed_list[i][0],'start':int(bed_list[i][1])-args.flank,'stop':int(bed_list[i][2])+args.flank,'strand':bed_list[i][3]})).upper(), range(len(bed_list)))
	if args.outfile:
		args.outfile.write('\n'.join(result))
	else:
		print ''.join(result)

if __name__=='__main__':
	args=get_opt()
	run(args)