import Bio
from Bio.Seq import Seq


class Protein:
    global index, locus_tag, name_y

    def __init__(self, index, locus_tag, name_y):
        self.index = index
        self.locus_tag = locus_tag
        self.name_y = name_y

    # def test(self):

# a = ProteinFilesManager()
# a.test()
