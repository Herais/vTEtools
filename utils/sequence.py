import Bio
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
    def get_seq_upstream_of_PAM(seq:str, len_guide=20, pam="*GG"):
        """
        """
        len_seq = len(seq)
        len_pam = len(pam)
        pat_pam = re.sub(r"\*", "[ATCGN]", pam)
        
        ls_seq = []
        for i in range(len_seq-len_pam):
            _pam = seq[i:i+len_pam]
            if re.match(pat_pam, _pam):
                if (i - len_guide) > -1:
                    pat_seq = seq[i-len_guide:i]
                    #print(i, re.match(pat_pam, _pam), pat_seq)
                    ls_seq.append((pat_seq, _pam))

        return ls_seq
    
    @staticmethod
    def get_tmpfp_bed4(S, fp="./tmp.bed"):
    
        df_coord = S.str.split('[:-]', expand=True)
        df_coord.columns = ['Chromosome', 'Start', 'End']
        df_coord['repName'] = df_coord.index
        
        bed = pybed.BedFrame.from_frame(meta=[], data=df_coord)
        bed.to_file(fp)
        return fp

    @staticmethod
    def get_flop_coordinates(S, fp_genome:str, len_flop=50, fptmp="./tmp.bed"):

        myDNA = DNA()
        
        fp_bed = myDNA.get_tmpfp_bed4(S, fp=fptmp)
        cmd = "bedtools slop -b {} -i {} -g {}".format(len_flop, fp_bed, fp_genome)
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        df_flop_coord = pd.read_csv(StringIO(res.stdout), sep='\t', header=None).set_index(3).sort_index()
        df_flop_coord.index.names = ['names']

        df_flop_coord.columns = [['meta']*3, ['chr_flop', 'start_flop', 'end_flop']]
        return df_flop_coord
    

class RNA(object):
    """
    """
    def __init__(self):
        """
        Parameters
        ----------

        """
        super(RNA, self).__init__()