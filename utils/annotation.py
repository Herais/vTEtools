import utils
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqR


class ANNOTATION(object):
    """
    """

    def __init__(self):
        """
        Parameters
        ----------

        """
        super(ANNOTATION, self).__init__()

    def get_entresummary(gene_name, organism="Mus Musculus", email=""):
        """
        """
        Entrez.email = email
        query = f'{gene_name}[Gene Name] AND {organism}[Organism]'
        handle = Entrez.esearch(db='gene', term=query)
        record = Entrez.read(handle)
        handle.close()

        NCBI_ids = record['IdList']
        for id in NCBI_ids:
            handle = Entrez.esummary(db='gene', id=id)
            record = Entrez.read(handle)
            # print(record)
        # print(record['Summary'])
        return record


    Entrez.email = "vee.xu.x@gmail.com"
    # stream = Entrez.einfo(db="gene")
    # record = Entrez.read(stream)
