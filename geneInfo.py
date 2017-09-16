#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: WangYaheng
Date: Sep-15-2017 19:20
"""


import argparse, sys
import MySQLdb

def get_opt():
	'''Check and parsing the opts'''
	parser=argparse.ArgumentParser(
				prog='geneInfo', 
				description='geneInfo: A program to get gene BED file from UCSC databse.', 
				usage='%(prog)s.py [options] -i geneSymbol'
			)
	parser.add_argument('-i','--genesymbol',nargs='?',type=str,default=sys.stdin, help="geneSymbol for searching in UCSC databases. [Required]",required=True)
	parser.add_argument('-d','--database',nargs='?',type=str,default='hg38', help="Database name for specificity species. Human:'hg38',mouse:'mm10'.default is 'hg38'.",required=False)
	parser.add_argument('-o','--outfile',nargs='?',type=argparse.FileType('w'),default=sys.stdout,help="Output file name for storing the results, default is stdout.",required=False)
	args=parser.parse_args()
	return args
	
def conn_ucsc(args):
	conn=MySQLdb.connect(host='genome-mysql.soe.ucsc.edu',port=3306,user='genome',passwd='',db=args.database)
	cur=conn.cursor()
	sql_cmd='''select a.chrom,a.strand,a.exonStarts,a.exonEnds,b.geneSymbol,b.description from knownGene as a,kgXref as b where a.name=b.kgID and lower(geneSymbol)='{0}' order by a.txEnd-a.txStart desc limit 1'''
	count=cur.execute(sql_cmd.format(args.genesymbol))
	geneInfo=cur.fetchone()
	if cur:
		cur.close()	
	return geneInfo
	
		
def run_geneInfo(geneInfo, args):
	if not geneInfo:
		print "Please check your input, it's not a geneSymbol."
		exit()
	chrom,strand,starts,ends,genesymbol,description=geneInfo[0],geneInfo[1],geneInfo[2].strip(',').split(','),geneInfo[3].strip(',').split(','),geneInfo[4],geneInfo[5]
	if args.outfile:
		[args.outfile.write("{0}\t{1}\t{2}\t{3}\n".format(chrom,start,end,strand)) for start,end in zip(starts, ends)]
	else:
		for start,end in zip(starts, ends):
			print "{0}\t{1}\t{2}\t{3}".format(chrom,start,end,strand)
			
if __name__=='__main__':
	args=get_opt()
	geneInfo=conn_ucsc(args)
	run_geneInfo(geneInfo, args)
