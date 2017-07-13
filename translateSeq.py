#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author:    Zhihao Xie  \(#^o^)/
# Date:      2017.05.31
# Version:   v1.0.0
# CopyRight: Copyright ©Zhihao Xie, All rights reserved.

import os, sys, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

if len(sys.argv) < 3:
    sys.stderr.write("Usage: python3 %s <gene_seq> <translate_outSeq> [code_table]\n" % sys.argv[0])
    sys.exit()

seqFile = os.path.abspath(sys.argv[1])
outSeq = os.path.abspath(sys.argv[2])
if len(sys.argv) == 4:
    code_table = sys.argv[3]
else:
    code_table = 11

handle = open(outSeq, "w")
for seq_record in SeqIO.parse(seqFile, "fasta"):
    seq_id = seq_record.id
    seq_desc = seq_record.description.replace(seq_record.id, "")
    seq_seq = seq_record.seq
    protein_seq = seq_seq.translate(table=code_table)
    pro_length = len(protein_seq)-1
    if re.search(r"\d+_?nt", seq_desc):
        seq_desc = re.sub(r"\d+_?nt", str(pro_length)+"_aa", seq_desc)
    protein_record = SeqRecord(Seq(str(protein_seq).rstrip("*"), IUPAC.protein), id=seq_id, description=seq_desc)
    SeqIO.write(protein_record, handle, "fasta")

handle.close()

