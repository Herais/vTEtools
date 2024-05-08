import Bio
from Bio import SeqIO
from Bio.Seq import Seq


class DNA(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(DNA, self).__init__()

    @staticmethod
    def get_seq_upstream_of_PAM(seq:str, len_guide=20, pam="*GG"):
        len_seq = len(seq)
        len_pam = len(pam)
        pat_pam = re.sub(r"\*", "[ATCG]", pam)
        ls_gRNA_candidates = []

        for i in range(len_seq-len_pam):
            _pam = seq[i:i+len_pam]
            if re.match(pat_pam, _pam):
                if (i - len_guide) > -1:
                    pat_seq = seq[i-len_guide:i]
                    #print(i, re.match(pat_pam, _pam), pat_seq)
                    ls_gRNA_candidates.append((pat_seq, _pam))
        return ls_gRNA_candidates