import sys
import Bio
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from decimal import *
getcontext().prec = 3

def PA_cgMLST_Alleles(input_fasta, output_folder):
    """Makes multifasta files for cgMLST alleles"""
    Gene_List = list(SeqIO.parse(input_fasta, 'fasta'))
    for genes in Gene_List:
        File_Name = output_folder + genes.id[0:-2] + '_alleles.fasta'
        output_handle = open(File_Name, 'w')
        SeqIO.write(genes, output_handle, 'fasta')
        output_handle.close()

Fasta = sys.argv[1]
Folder = sys.argv[2]
if Folder[-1] != '/':
    Folder = Folder + '/'

PA_cgMLST_Alleles(Fasta, Folder)
