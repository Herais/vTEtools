import os, sys
import pandas as pd
import numpy as np
import re
import collections
import Bio
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from fuc import pybed


class SGRNA(object):
    """
    """

    def __init__(self):
        """
        Parameters
        ----------

        """
        super(SGRNA, self).__init__()


    @staticmethod
    def get_seq_upstream_of_PAM(seq:str, len_guide=20, pam="NGG")->list:
        """
        @param seq: DNA sequence
        @param len_guide: length of the guide sequence
        @param pam: PAM sequence
        @return: list of guide
        """
        len_seq = len(seq)
        len_pam = len(pam)
        pat_pam = re.sub("N", "[ATCG]", pam)
        #print(pat_pam)
        
        ls_seq = []
        for i in range(len_seq-len_pam):
            _pam = seq[i:i+len_pam]
            #print(_pam)
            if re.match(pat_pam, _pam):
                if (i - len_guide) > -1:
                    pat_seq = seq[i-len_guide:i]
                    #print(i, re.match(pat_pam, _pam), pat_seq)
                    ls_seq.append((pat_seq, _pam))

        return ls_seq
    
    @staticmethod
    def get_seq_without_PAM_restriction(seq, len_guide=17):
        L = len(seq)
        i=0
        ls_seq = []
        while i < (L-len_guide):
            ls_seq.append(seq[i:i+len_guide])
            i += 1
        return ls_seq

    @staticmethod
    def get_all_gRNA_with_PAM(
        S_seq,
        len_guide=17,
        fp_genome="/gladstone/alexanian/datasets-online/VX_reference_and_indexes/hg38/hg38.chrom.sizes",
        fp_fasta = "/gladstone/alexanian/datasets-online/VX_reference_and_indexes/hg38/GRCh38.primary_assembly.genome.fa"
    ) -> pd.DataFrame:

        """
        """
        mySGRNA = SGRNA()

        S_gRNA_watson = S_seq \
            .apply(lambda x: mySGRNA.get_seq_upstream_of_PAM(
                x,
                len_guide=len_guide,)) \
            .apply(lambda x: [y[0] for y in x])
        
        S_gRNA_cricket = S_seq \
            .apply(myDNA.get_reverse_complement) \
            .apply(lambda x: 
                mySGRNA.get_seq_upstream_of_PAM(x, len_guide=len_guide,)) \
            .apply(lambda x: [y[0] for y in x])

        df_sgRNA_all = pd.DataFrame([S_gRNA_watson, S_gRNA_cricket]).T
        df_sgRNA_all.columns = ['gRNA_watson', 'gRNA_cricket']
        #df_sgRNA_all['name_unique'] = df_loci['name_unique']
        
        df_sgRNA_all['gRNA_bothStrand'] = df_sgRNA_all[['gRNA_watson', 'gRNA_cricket']] \
            .apply(lambda x: x['gRNA_watson'] + x['gRNA_cricket'], axis=1)
        
        df_sgRNA_all['gRNA_bothStrand_unique'] = df_sgRNA_all['gRNA_bothStrand'].apply(set)
        
        df_sgRNA_all['counter_gRNA_bothStrand'] = df_sgRNA_all['gRNA_bothStrand'] \
            .apply(Counter)
        
        df_sgRNA_all['counter_gRNA_bothStrand_unique'] = df_sgRNA_all['gRNA_bothStrand_unique'] \
            .apply(Counter)
        
        #Counter_sgRNA_places = df_sgRNA_all['counter_gRNA_bothStrand'].sum()
        #Counter_sgRNA_teloci =  df_sgRNA_all['counter_gRNA_bothStrand_unique'].sum()

        return df_sgRNA_all