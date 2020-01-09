import Bio
from Bio.Seq import Seq
from Protein import Protein


class Strain:

    index = 0
    name = ""
    proteins = {}

    def __init__(self, index, name, dict_proteins):
        self.index = index
        self.name = name
        self.proteins = dict_proteins

# a = Strain(0, "aa", None)

