# -*- coding: utf-8 -*-
from __future__ import print_function
# snp画图和按阈值过滤

__author__ = 'Yaheng Wang'
__date__ = '2018-01-03'
__version_ = '1.2'

import csv
import os
import xlrd
from matplotlib import pyplot as plt

# 画单个rs号散点图
def plot(alle_percent, figname):
	plt.plot(alle_percent, 'o', markersize='2')
	plt.xlim(xmin=0)
	plt.ylim(ymin=0, ymax=1)
	plt.savefig(figname)
	plt.close()


# 读取测序分析后文件
def read_ngs(fh):
	of = open(fh, 'r')
	line_dict = csv.DictReader(of, dialect='excel-tab')
	fieldnames = line_dict.fieldnames+ ['alle_percent']
	# 按rsID排序
	line_list = sorted(line_dict, key=lambda x:x['rsID'])
	# 得到所有rs号
	rs_nums = list(set(map(lambda line:line['rsID'], line_list)))
	# 给每个信息添加alle_percent项
	for line in line_list:
		line['alle_percent'] = float(line[line['reference base']])/float(line['total mapped reads'])
	of.close()
	return fieldnames, rs_nums, line_list


# 读取过滤条件文件excel
def read_filterexcel(fh):
	table = xlrd.open_workbook(fh)
	data = table.sheet_by_index(0)
	filter_dict = {}
	for i in range(data.nrows):
		row = data.row_values(i)
		filter_dict[row[0]] = list(map(lambda x:float(x), row[1:]))
	return filter_dict


# 画过滤前snp散点分布图
def before_plt(path, rs_nums, line_list):
	for rs in rs_nums:
		alle_percent = list(map(lambda y:y['alle_percent'], filter(lambda x:x['rsID']==rs, line_list)))
		figname = os.path.join(path, rs+'.jpg')
		plot(alle_percent, figname)


# 画过滤后snp散点分布图
def after_plt(path, rs_nums, line_list, filter_dict):
	filter_list = list()
	for rs in rs_nums:
		threshold = filter_dict.get(rs, [0, 0, 0, 0])
		filter_line = list(filter(lambda x:x['rsID']==rs and not (threshold[0]<=x['alle_percent']<=threshold[1] or threshold[2]<=x['alle_percent']<=threshold[3]), line_list))
		alle_percent = list(map(lambda y:y['alle_percent'], filter_line))
		figname = os.path.join(path, rs+'.jpg')
		plot(alle_percent, figname)
		filter_list += filter_line
	return filter_list


# 输出过滤后的文件
def writer_csv(path, fieldnames, filter_list):
	with open(os.path.join(path, 'filter.csv'), 'w', newline='') as fh:
		line_dict = csv.DictWriter(fh, fieldnames=fieldnames)
		line_dict.writeheader()
		line_dict.writerows(filter_list)


def main(ngsf, filterf=None):
	if not filterf:
		bpath = os.path.join(os.path.dirname(ngsf), 'before_jpg')
		if not os.path.exists(bpath):
			os.mkdir(bpath)
		fieldnames, rs_nums, line_list = read_ngs(ngsf)
		before_plt(bpath, rs_nums, line_list)
	else:
		fieldnames, rs_nums, line_list = read_ngs(ngsf)
		filter_dict = read_filterexcel(filterf)
		apath = os.path.join(os.path.dirname(ngsf), 'after_jpg')
		if not os.path.exists(apath):
			os.mkdir(apath)
		filter_list = after_plt(apath, rs_nums, line_list, filter_dict)
		writer_csv(apath, fieldnames, filter_list)


if __name__ == '__main__':
	import sys
	if len(sys.argv) == 3:
		ngsf = sys.argv[1]
		filterf = sys.argv[2]
		main(ngsf, filterf)
	elif len(sys.argv) == 2:
		ngsf = sys.argv[1]
		main(ngsf)
	else:
		print("Error, please check your input args.")