#!/usr/bin/env python3

import Bio
import sys
import glob, os

from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.codonalign import build
from Bio.codonalign import CodonSeq, CodonAlignment
from Bio.codonalign.codonseq import _get_codon_list, CodonSeq, cal_dn_ds
from Bio.Data import CodonTable
import codonalignmentJDY

def MK_orthologs(pop1,pop2,outgroup=None):
     results=[]
     for file in glob.glob("*.gap.fasta"):
         print(file)
         seqs=AlignIO.read(file, 'fasta')
         codon_aln = CodonAlignment()
         codon_aln = codon_aln.from_msa(seqs)
         print(codon_aln)
         ###https://github.com/biopython/biopython/blob/master/Doc/doc.rst
         #dn_matrix, ds_matrix = codon_aln.get_dn_ds_matrix()
         #print(dn_matrix)
         #### population info
         ptt=[]
         ptm=[]
         for record in codon_aln:
             if record.id.startswith("M"):
                 ptt.append(record)
             else:
                 ptm.append(record)
         ptt=codonalignmentJDY.CodonAlignment(ptt)
         ptm=codonalignmentJDY.CodonAlignment(ptm)
         ### set order of populations to be examined
         if pop1=="ptt":
             result=codonalignmentJDY.mktest([ptt,ptm])
         elif pop1=="ptm":
             result=codonalignmentJDY.mktest([ptm,ptt])
         result1="\t".join([str(s) for s in result])
         result2="\t".join([str(s) for s in [file,result1]])
         results.append(result2)
     print(results)
     try:
         filename = "MK_results." + pop1+pop2+outgroup + ".txt"
         with open(filename, 'w') as f:
             for item in results:
                 f.write("%s\n" % item)
     except:
         filename = "MK_results." + pop1+pop2 + ".txt"
         with open(filename, 'w') as f:
             for item in results:
                 f.write("%s\n" % item)

print("Population1 PTT and Population2 PTM")
MK_orthologs(pop1="ptt",pop2="ptm")

print("Population1 PTM and Population2 PTT")
MK_orthologs(pop1="ptm",pop2="ptt")

