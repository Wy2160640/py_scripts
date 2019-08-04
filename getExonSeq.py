#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File	:	getExonSeq
@Time	:	2019/8/4
@Author	:	Yaheng	Wang 
@Version	:	1.0
@Contact	:	yhwang91@163.com
@License 	:	(C)Copyright 2018-2019, yhwang
@Link    : https://www.cnblogs.com/yahengwang/
'''


import argparse
import os
import sys
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import BLOB, CHAR, Column, Enum,String, Integer, create_engine, or_
from sqlalchemy.orm import sessionmaker

PY3 = sys.version_info > (3, )

parser = argparse.ArgumentParser(prog=os.path.basename(__file__))
parser.add_argument("-dialect", "--dialect", nargs='?', type=str, default="mysql")
parser.add_argument("-driver", "--driver", nargs='?', type=str, default="mysqldb")
parser.add_argument("-username", "--username", nargs='?', type=str, default="genome")
parser.add_argument("-password", "--password", nargs='?', type=str, default="")
parser.add_argument("-host", "--host", nargs='?', type=str, default="genome-mysql.soe.ucsc.edu")
parser.add_argument("-port", "--port", nargs='?', type=str, default="3306")
parser.add_argument("-database", "--database", nargs='?', type=str, default="hg38")
parser.add_argument("-gene", "--gene", nargs='+', type=str, required=True)
parser.add_argument("-flank", "--flank", nargs='?', type=int, default=0)
parser.add_argument("-format", "--format", nargs='?', choices=["FASTA", "BED"], default="FASTA")
params = parser.parse_args()

DB_URI = "{params.dialect}+{params.driver}://{params.username}:{params.password}@{params.host}:{params.port}/{params.database}?charset=utf8".format(params=params)
SEQ_URI = "http://genome.ucsc.edu/cgi-bin/das/{params.database}/dna?segment=".format(params=params)
engine = create_engine(DB_URI)
Base = declarative_base()

if PY3:
    from urllib.request import urlopen
    revcomp = lambda seq: seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
else:
    from urllib import urlopen
    from string import maketrans
    revcomp = lambda seq: seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


def _seq_from_xml(xml):
    start = xml.find(">", xml.find("<DNA")) + 1
    end = xml.rfind("</DNA>")
    return xml[start:end].replace(' ', '').replace('\n', '').strip()


def sequence(url, chrom, start, end):
    url = url + "%s:%i,%i" % (chrom, start, end)
    xml = urlopen(url).read().decode("utf-8")
    return _seq_from_xml(xml)


class Interval(object):
    def __init__(self, refname, symbol, chrom, start, end, strand):
        self.refname = refname
        self.symbol = symbol
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

        self._sequence = None
        self._revcomp = None

    @property
    def sequence(self):
        if self._sequence is None:
            self._sequence = sequence(SEQ_URI, self.chrom, self.start, self.end)
        return self._sequence

    @property
    def revcomp(self):
        if self._revcomp is None:
            self._revcomp = revcomp(self.sequence)
        return self._revcomp

    def __repr__(self):
        return "{x.symbol} {x.chrom}:{x.start}{x.strand}{x.end}".format(x=self)

    def __str__(self):
        return self.__repr__()


class RefGene(Base):
    __tablename__ = "refGene"

    bin = Column(Integer)
    name = Column(String(20), primary_key=True)
    chrom = Column(String(10))
    strand = Column(CHAR)
    txStart = Column(Integer)
    txEnd = Column(Integer)
    cdsStart = Column(Integer)
    cdsEnd = Column(Integer)
    exonCount = Column(Integer)
    exonStarts = Column(BLOB)
    exonEnds = Column(BLOB)
    score = Column(Integer)
    name2 = Column(String(20))
    cdsStartStat = Column(Enum("none", "unk", "incmpl", "cmpl"))
    cdsEndStat = Column(Enum("none", "unk", "incmpl", "cmpl"))
    exonFrames = Column(BLOB)

    def exonstarts(self, flank=0):
        return [int(i) - flank for i in self.exonStarts.rstrip(b",").split(b",")]

    def exonends(self, flank=0):
        return [int(i) + flank for i in self.exonEnds.rstrip(b",").split(b",")]

    def exonseqs(self, flank=0):
        for start, end in zip(self.exonstarts(flank), self.exonends(flank)):
            seq = Interval(self.name, self.name2, self.chrom, start, end, self.strand)
            yield seq


Session = sessionmaker(bind=engine)
session = Session()


def single_gene_query(gene):
    q = session.query(RefGene).filter(or_(RefGene.name == gene, RefGene.name2 == gene))\
        .order_by(RefGene.txEnd - RefGene.txStart).first()
    if q is None:
        sys.stderr.write("%s is not found in the database!\n" % gene)
        sys.exit()
    return q


def main(params):
    objs = [single_gene_query(gene) for gene in params.gene]
    for obj in objs:
        for exon in obj.exonseqs(params.flank):
            if params.format == "BED":
                line = "{x.chrom}\t{x.start}\t{x.end}\t{x.strand}\t{x.sequence}\n" \
                    if exon.strand == "+" else "{x.chrom}\t{x.start}\t{x.end}\t{x.strand}\t{x.revcomp}"
            else:
                line = ">{x.symbol} {x.chrom}:{x.start}{x.strand}{x.end}\n{x.sequence}\n"  \
                    if exon.strand == "+" else ">{x.symbol} {x.chrom}:{x.start}{x.strand}{x.end}\n{x.revcomp}\n"
            line = line.format(x=exon)
            sys.stdout.write(line)


if __name__ == "__main__":
    main(params)
