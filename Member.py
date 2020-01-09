from typing import Any

import Bio
from Bio.Seq import Seq


class Member:

    global strain_index
    global protein_index
    global length
    global identity
    global represntative

    def __init__(self, strain_index, protein_index, length, identity, represntative):
        self.strain_index = strain_index
        self.protein_index = protein_index
        self.length = length
        self.identity = identity
        self.represntative = represntative


    @property
    def getStrainInd(self):
        return self.strain_index

    @property
    def getProteinInd(self):
        return self.protein_index

    @property
    def getLength(self):
        return self.length

    @property
    def getIdentity(self):
        return self.identity

    @property
    def getRepresntative(self):
        return self.represntative


# a = Strain(0, "aa", None)

