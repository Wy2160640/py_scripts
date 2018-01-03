# -*- coding:utf-8 -*-
from __future__ import print_function
import primer3
import math

__author__ 	= 'WangYaheng'
__date__ 	= '2017-12-06'
__version__ = '1.0'

# read csv file
def parse_csv(fh):
	import csv
	file = open(fh, 'r')
	csvfile = csv.reader(file, dialect='excel')
	primer_dict = {}
	for row in csvfile:
		primer_dict.setdefault(row[0]+'_fp', row[1])
		primer_dict.setdefault(row[0]+'_rp', row[2])
	file.close()
	return primer_dict


def calcHairpin(seq, tm_Threshold=47, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.25, dna_conc=50.0, temp_c=37, max_loop=30):
	thermoresult = primer3.calcHairpin(seq, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop)
	if thermoresult.tm >= tm_Threshold:
		return thermoresult.tm, thermoresult.dg/1000
	else:
		return -1, -1


def calcDimer(seq1, seq2, dg_Threshold=8, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.25, dna_conc=50.0, temp_c=37, max_loop=30):
	thermoresult = primer3.calcHeterodimer(seq1, seq2, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c, max_loop)
	if math.fabs(thermoresult.dg/1000) >= dg_Threshold:
		return thermoresult.tm, thermoresult.dg/1000
	else:
		return -1, -1

def hairpin_analysis(seqs_dict):
	hairpin_list =[]
	keys = sorted(list(seqs_dict.keys()))
	import math
	for i in range(len(keys)):
		tm, dg = calcHairpin(seqs_dict[keys[i]])
		if tm!=-1:
			hairpin_list.append([keys[i], tm, dg])
	return hairpin_list


def dimer_analysis(seqs_dict):
	dimer_list = []
	keys = sorted(list(seqs_dict.keys()))
	import math
	for i in range(len(keys)):
		for j in range(i, len(keys)):
			tm, dg = calcDimer(seqs_dict[keys[i]], seqs_dict[keys[j]])
			if dg!=-1:
				dimer_list.append([keys[i], keys[j], tm, dg])
	return dimer_list


if __name__ == '__main__':
	import sys, csv, os
	primers_csv = sys.argv[1]
	primers_dict = parse_csv(primers_csv)

	dimer = dimer_analysis(primers_dict)
	check_file = open('_dimerchcek'.join(os.path.splitext(primers_csv)), 'w', newline='')
	csv_check_file = csv.writer(check_file, dialect='excel')
	for row in dimer:
		csv_check_file.writerow(row)