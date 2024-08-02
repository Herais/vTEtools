import os, sys
import subprocess
from io import StringIO
import Bio
import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from fuc import pybed


class DNA(object):
    """
    """
 
    def __init__(self):
        """
        Parameters
        ----------

        """
        super(DNA, self).__init__()
    
    @staticmethod
    def get_tmpfp_bed4(S, fp="default"):

        if fp == "default":
            dname = "tmp_vTEtools"
            os.makedirs(dname, exist_ok=True)
            fp = "./tmp_vTEtools/tmp.bed"
    
        df_coord = S.str.split('[:-]', expand=True)
        df_coord.columns = ['Chromosome', 'Start', 'End']
        df_coord['repName'] = df_coord.index
        
        bed = pybed.BedFrame.from_frame(meta=[], data=df_coord)
        bed.to_file(fp)
        return fp
    
    def get_tmpfp_fasta(S, fp="default"):
            
            if fp == "default":
                dname = "tmp_vTEtools"
                os.makedirs(dname, exist_ok=True)
                fp = "./tmp_vTEtools/tmp.fasta"
            

            return fp

    @staticmethod
    def get_sequence_from_coordinates(coords, fp_fasta):
        
        myDNA = DNA()
        fp = myDNA.get_tmpfp_bed4(coords)
        cmd = "bedtools getfasta -fi {} -bed {} -s -bedOut".format(fp_fasta, fp)
        res =  subprocess.run(cmd, shell=True, capture_output=True, text=True)
        dftmp = pd.read_csv(StringIO(res.stdout), sep='\t', header=None).set_index(3)

        return dftmp[4]

    @staticmethod
    def get_flop_coordinates(S, fp_genome:str, len_flop=50, fptmp="default"):

        myDNA = DNA()
        
        fp_bed = myDNA.get_tmpfp_bed4(S, fp=fptmp)
        cmd = "bedtools slop -b {} -i {} -g {}".format(len_flop, fp_bed, fp_genome)
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        df_flop_coord = pd.read_csv(StringIO(res.stdout), sep='\t', header=None).set_index(3).sort_index()
        df_flop_coord.index.names = ['names']

        df_flop_coord.columns = [['meta']*3, ['chr_flop', 'start_flop', 'end_flop']]
        return df_flop_coord
    
    @staticmethod
    def get_reverse_complement(seq):
        """
        """
        return str(Seq(seq).reverse_complement())
    

class RNA(object):
    """
    """
    def __init__(self):
        """
        Parameters
        ----------

        """
        super(RNA, self).__init__()