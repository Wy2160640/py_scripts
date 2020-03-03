#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   identity.py
@Time    :   2020/03/03 10:15:13
@Author  :   jjbioinfo
@Version :   1.0
@Contact :   jjbioinfo@163.com
@License :   (C)Copyright 2019-2020, 生信小栈
@Desc    :   None
'''

# here put the import lib
import argparse
import csv
import os
import random
import re
import string
import sys
import subprocess
from collections import OrderedDict


class Alignment(object):
    def __init__(self, q_id, s_id, idy, align_len, mis, gap, q_s, q_e, s_s, s_e, ev, score):
        """
        q_id, query sequence ID
        s_id, target sequence ID identification on the alignment
        idy, sequence alignment percent identity
        align_len, the length of the matching region that matches
        mis, the mismatch number of the matching region
        gap, the number of gaps in the comparison area
        q_s, starting position of the alignment region on the query sequence
        q_e, termination point of the alignment region on the query sequence
        s_s, the start position of the alignment region on the subject sequence
        s_e, the end position of the alignment region on the subject sequence
        ev, the expected value of the comparison result is explained as how many random comparisons can occur once this score, the smaller the Evalue, indicating that this situation is less likely to occur, and then it means that this is more likely to be a true similar sequence
        score, bit score of the alignment result
        """
        self.q_id = q_id
        self.s_id = s_id
        self.idy = idy
        self.align_len = align_len
        self.mis = mis
        self.gap = gap
        self.q_s = q_s
        self.q_e = q_e
        self.s_s = s_s
        self.s_e = s_e
        self.score = score

        self.name, self.start, self.end = self._parse_query_id()

    def _parse_query_id(self):
        items = self.q_id.rsplit(":", 1)
        name = items[0]
        start, end = items[1].split("-")
        return name, int(start), int(end)

    def __str__(self):
        return "{self.q_id}:{self.q_s}-{self.q_e}({self.idy})".format(self=self)

    def __repr__(self):
        return "{self.q_id}:{self.q_s}-{self.q_e}({self.idy})".format(self=self)

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and self.q_id == other.q_id)


def wrap_blat(query, db):
    """
    query, query is each either a .fa, .nib or .2bit file
    db, db is each either a .fa, .nib or .2bit file
    """
    blat = os.path.join(os.path.dirname(__file__), 'bin', 'blat')
    if not os.path.exists(blat):
        blat = 'blat'
    cmd = "{blat} -noHead -t=dna -q=dna -out=blast8  -minIdentity=0 {db} {query} stdout".format(blat=blat, db=db, query=query)
    stdout, stderr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    content = stdout.decode()
    data = []
    if content:
        for line in content.split(os.linesep):
            if line:
                row = line.split("\t")
                data.append(Alignment(*row))
    
    data.sort(key=lambda x:(x.name, x.start, x.end, x.align_len), reverse=False)
    return data


def filter_blat(data, out=None):
    outfh = open(out, "w") if out else sys.stdout
    new_data = []
    for align in data:
        if align not in new_data:
            new_data.append(align)
            line = "{align.name}\t{align.start}\t{align.end}\t{align.idy}\n".format(align=align)
            outfh.write(line)
    outfh.close()


class Reader():
    def __init__(self, infasta, window, step, outfasta=None):
        """
        infasta, genome FASTA file
        outfasta, kmer FTASA file
        window, kmer size
        step, step size
        """

        self.infasta = infasta
        self.outfasta = open(outfasta, "w") if outfasta else sys.stdout
        self.window = window
        self.step = step

        self.data = self._load_data()

    def _load_data(self):
        data = OrderedDict()
        with open(self.infasta) as infasta:
            for i, line in enumerate(infasta):
                if line[0] == ">":
                    seq = line[1:].split()[0]
                    if seq not in data:
                        data[seq] = '' 
                    else:
                        assert i == 0, 'There may be a header without a sequence at line {}'.format(i)
                else:
                    data[seq] += line.strip()
        return data

    def _kmer_genome(self):
        try:
            for header, seq in self.data.items():
                for i in range(0, len(seq), self.step):
                    subseq = seq[i : i+self.window]
                    subheader = ">{0}:{1}-{2}\n".format(header, i, i+self.window)
                    self.outfasta.write(subheader + subseq + "\n")
            sys.stdout.write("Genome {self.infasta} have been kmer with window {self.window}bp, step {self.step}bp!\n".format(self=self))
        except Exception as e:
            sys.stderr.write(e)
        finally:
            self.outfasta.close()

    def run(self):
        self._kmer_genome()
        

def main():
    args = argparse.ArgumentParser(description='Important, enviroment must have installed <blat>, calculate query genome identity to the reference genome')
    args.add_argument('-d', '--db', type=str, help='database file, a .fa, .nib or .2bit file', required=True)
    args.add_argument('-q', '--query', type=str, help='query file, a .fa, .nib or .2bit file', required=True)
    args.add_argument('-w', '--window', type=int, help='sliding window algorithm, window size(bp)', default=1000)
    args.add_argument('-s', '--step', type=int, help='sliding window algorithm, step size(bp)', default=200)
    args.add_argument('-o', '--out', type=argparse.FileType(mode="w"), help='the output file', default=sys.stdout)

    params = args.parse_args()

    salt =  ''.join(random.sample(string.ascii_letters + string.digits, 24))
    outfasta = os.path.join("/tmp", salt)
    Reader(params.query, params.window, params.step, outfasta).run()

    os.remove(outfasta)
    alignments = wrap_blat(outfasta, params.db)
    #print(alignments)
    filter_blat(alignments, out)


if __name__ == "__main__":
    main()
