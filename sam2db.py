#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author: WangYaheng
Date: Sep-30-2017 13:20
"""

import argparse, sys, os, subprocess
import pandas, sqlite3

def get_opt():
	'''Check and parsing the opts'''
	parser=argparse.ArgumentParser(
				prog='sam2db', 
				description='sam2db: A program to transfer sam to databases.', 
				usage='%(prog)s.py [options] -f sam'
			)
	parser.add_argument('-f','--sam_file',nargs='?',type=str,default=sys.stdin,help="The sam file path. [Required]",required=True)
	parser.add_argument('-d','--database',nargs='?',type=str,default='sqlite3.db',help='The name of the database to save sam data. default="sqlit3.db"',required=False)
	
	args=parser.parse_args()
	return args
  


def run(args):
	# sam header
	header=[
		'QNAME',
		'FLAG',
		'RNAME',
		'POS',
		'MAPQ',
		'CIGAR',
		'MENM',
		'MPOS',
		'ISIZE',
		'SEQ',
		'QUAL'
		]

	# reading sam with pandas
	sam=pandas.read_table(args.sam_file, 		
		header=None,					
		comment="@",					
		usecols=range(11),				
		error_bad_lines=False)

	sam.columns=header

	# 连接数据库
	conn=sqlite3.connect(args.database)
	print "Opened database successfully"

	# 将pandas数据转到sql中
	sam.to_sql(os.path.splitext(os.path.basename(args.sam_file))[0], 
		conn, 
		if_exists='replace', 
		index=True)
	print "ok"
	
	conn.commit()
	conn.close()


if __name__=='__main__':
	args=get_opt()
	run(args)
