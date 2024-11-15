import os, subprocess
import re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
from anndata import AnnData


class TE(object):
    """
    """

    def __init__(self):
        """
        Parameters
        ----------

        """
        super(TE, self).__init__()

    @staticmethod
    def perform_grouping(df, print_summary=True):
        df = df.copy()
        # GROUP 0
        # ungrouped
        df[('meta', 'selection_group')] = 0

        # GROUP 1
        # in vivo, 'TAC_vs_Sham_up' > 0
        df.loc[df['ranking']['TAC_vs_Sham_up'] >
            0, ('meta', 'selection_group')] = 1

        # GROUP 2
        # in vivo, 'TAC_vs_Sham_down', and
        # in vivo, TACJQ1w vs. TACJQ1_up, and
        # in vitro, 'TGFb vs. Unstim H3K27ac_up'
        df_unassigned = df[df['meta']['selection_group'] == 0]
        idx_grp2a = df_unassigned[df_unassigned['ranking']
                                ['TACJQ1w_vs_TACJQ1_up'] == 1].index
        # print(len(idx_grp3a))
        df_grp2a = df_unassigned.loc[idx_grp2a]
        idx_grp2b = df_grp2a[df_grp2a['ranking']
                            ['TGFb vs. Unstim_h3k27ac'] == 1].index
        # print(len(idx_grp3b))
        df_grp2b = df_grp2a.loc[idx_grp2b]
        df.loc[idx_grp2b, ('meta', 'selection_group')] = 2

        # GROUP 3
        # in vivo, Brd4 FC >= 10, and
        # in vivo H3K27ac FC >= 0
        df_unassigned = df[df['meta']['selection_group'] == 0]
        idx_grp3a = df_unassigned[df_unassigned['TGFb vs. Unstim_Brd4']
                                ['logfoldchanges'] > 10].index
        df_grp3a = df_unassigned.loc[idx_grp3a,]  # shape (60, 70)
        idx_grp3b = df_grp3a[df_grp3a['TGFb vs. Unstim_h3k27ac']
                            ['logfoldchanges'] > 0].index
        df.loc[idx_grp3b, ('meta', 'selection_group')] = 3

        # GROUP 4
        # proximity, closest gene up is in screened genes, and
        # proximity, closest gene downstream is in screened genes
        df_unassigned = df[df['meta']['selection_group'] == 0]
        idx_grp4a = df_unassigned[df_unassigned['proximity']
                                ['closest_gene_us_is_screened']].index
        df_grp4a = df_unassigned.loc[idx_grp4a,]
        idx_grp4b = df_grp4a[df_grp4a['proximity']
                            ['closest_gene_ds_is_screened']].index
        idx_grp4 = list(idx_grp4a) + list(idx_grp4b)
        df.loc[idx_grp4, ('meta', 'selection_group')] = 4

        # GROUP 5
        # proximal to genes of relevant GO-term

        # GROUP 0
        df_unassigned = df[df['meta']['selection_group'] == 0]

        # print summary
        if print_summary:
            print('group 0: {}'.format(df_unassigned.shape[0]))
            print('group 1: {}'.format((df['meta']['selection_group'] == 1).sum()))
            print('group 2: {} (exclusive of previous group)'.format(
                df_grp2b.shape[0]))
            print('group 3: {} (exclusive of previous group)'.format(len(idx_grp3b)))
            print('group 4: {} (exclusive of previous group)'.format(len(idx_grp4)))

        return df

    @staticmethod
    def create_SeqRecord(row):
        seq = Seq(row['seg_for_alignment']) # 'sequence'
        id = str(row['name'])
        name = str(row['name'])
        quality = row['quality']
        record = SeqRecord(seq,id, name, quality)
        return record

    def create_record_seq(row, col='Sequence'):
        seq = Seq(row[col])
        id = "id_{}".format(row['index'])
        name = "id_{}".format(row['index'])
        description = ''
        record = SeqRecord(seq,id, name, description)
        return record


    @staticmethod
    def multialign_clustalw2(fp_fasta):
        """
        """
        cpath_clustal = "/wynton/home/alexanian/veexu/env_py3.10b/bin/clustalw2"
        cmd = ClustalwCommandline(cpath_clustal, infile=fp_fasta)
        print(cmd)
        stdout, stderr = cmd()
        align = AlignIO.read(fp_fasta.replace('.fasta', '.aln'), "clustal")
        alignment = AlignIO.read(open(fp_fasta.replace('.fasta', '.aln')), "clustal")
        print("Alignment length %i" % alignment.get_alignment_length())
        ls_seq = []

        for record in alignment:
            #print(record.seq + " " + record.id)
            ls_seq.append(str(record.seq))

        return ls_seq

    @staticmethod
    def get_seq_upstream_of_PAM(seq: str, len_guide=20, pam="*GG"):

        """
        """
        len_seq = len(seq)
        len_pam = len(pam)
        pat_pam = re.sub(r"\*", "[ATCG]", pam)
        ls_gRNA_candidates = []

        for i in range(len_seq-len_pam):
            _pam = seq[i:i+len_pam]
            if re.match(pat_pam, _pam):
                if (i - len_guide) > -1:
                    pat_seq = seq[i-len_guide:i]
                    # print(i, re.match(pat_pam, _pam), pat_seq)
                    ls_gRNA_candidates.append((pat_seq, _pam))

        return ls_gRNA_candidates
    
    @staticmethod
    def get_count_gRNA_loci(df, col='gRNA_20_bp'):
        """
        """
        df_gRNA = df[df['meta']['num_{}'.format(col)] > 0].copy()  # 295
        S = df_gRNA['meta'][col].apply(lambda x:
                                    Counter(set([g[0] for g in x])))
        return S.sum()

    @staticmethod
    def get_gcount_from_sequence(
        seq, 
        fp_genome='/gladstone/alexanian/datasets-online/VX_reference_and_indexes/hg38//Homo_sapiens.GRCh38.dna.primary_assembly_noScaffold.fa'
        ):
        
        cmd = 'grep -c "{}" {}'.format(seq, fp_genome)
        #print(cmd)
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        gcount = (res.stdout).rstrip()
        gcount = int(gcount)
    
        return gcount

    @staticmethod
    def generate_gRNA_by_length(seq, len_guide=20):
        """
        """

        assert (len_guide > len_seq,
                "The sequence of length {} is shorter than gRNA of length {}".format(len_seq, len_guide))

        ls_gRNA = []
        for i in range(len_seq-len_guide):
            ls_gRNA.append(seq[i:i+len_guide])

        return ls_gRNA

    @staticmethod
    def runBigWigAverageOverBed(
        ls_fp_bw,
        fp_bed,
        dpath_result=None, 
        cmd_bigWigAverageOverBed = "/gladstone/alexanian/datasets-online/VX_pyProjects/pjScAtacSeq/tools/bigWigAverageOverBed",
        bw_extension='bigWig',
        show_progress=False,):

        """
        """
        os.makedirs(dpath_result, exist_ok=True)

        for fp_bw in ls_fp_bw:
            bw = fp_bw.split('/')[-1]
            fout_tab = "{}/{}".format(dpath_result, re.sub(bw_extension, "tab", bw))
            fout_ra = "{}/{}".format(dpath_result, re.sub(bw_extension, "ra", bw))
            fout_bed = "{}/{}".format(dpath_result, re.sub(bw_extension, "bed", bw))
            job_name = re.sub(bw_extension, "", bw)
            
            cmd="{} {} {} {} -stats={} -bedOut={} -minMax | qsub -S /bin/bash -cwd -N {}".format(
                cmd_bigWigAverageOverBed,
                fp_bw,
                fp_bed,
                fout_tab, 
                fout_ra,
                fout_bed,
                job_name
            )
        
            res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if show_progress:
                print(res.stderr)
        print('done.')

    @staticmethod
    def split_dfTE_to_groups(df, repName=None, return_all=False):

        ret = {}

        if repName:
            df = df[df['repName'] == repName].copy()

        if return_all:
            ret['all'] = df.copy()

        df_distal = df[df['isin_windowed_gap']]
        ret['distal_to_pcgene'] = df_distal.copy()

        df_pcgene = df[df['isin_gene']]
        ret['pcgene'] = df_pcgene.copy()

        df_pcgene_exonic = df[df['isin_gene']][df['isin_exon']]
        ret['pcgene_exonic'] = df_pcgene_exonic.copy()
        df_pcgene_intronic = df[df['isin_gene']][~df['isin_exon']]
        ret['pcgene_intronic'] = df_pcgene_intronic.copy()

        df_transcripts_nonpc = df[df['isin_transcript']][~df['isin_gene']]
        ret['transcripts_nonpc'] = df_transcripts_nonpc.copy()

        df_transcript_nonpc_exonics = df_transcripts_nonpc[df_transcripts_nonpc['isin_exon']]
        ret['transcript_nonpc_exonic'] = df_transcript_nonpc_exonics.copy()

        df_transcript_nonpc_intronic = df_transcripts_nonpc[~df_transcripts_nonpc['isin_exon']]
        ret['transcript_nonpc_intronic'] = df_transcript_nonpc_intronic.copy()

        df_buffer = df[~df['isin_windowed_gap']][~df['isin_transcript']]
        ret['in_buffer'] = df_buffer.copy()

        return ret
