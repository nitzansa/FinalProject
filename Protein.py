import Bio
from Bio.Seq import Seq


class Protein:
    index = 0
    locus_tag = ""

    def __init__(self, index, locus_tag):
        self.index = index
        self.locus_tag = locus_tag

    # def test(self):

# a = ProteinFilesManager()
# a.test()
