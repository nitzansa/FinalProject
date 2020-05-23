import Bio
from Bio.Seq import Seq
from Protein import Protein


class Strain:


    index = 0
    name = ""
    proteins = {}
    numOfGenes = 0
    numOfCoreGenes = 0
    numOfSingletons = 0
    clusterList = set()

    def __init__(self, index, name, dict_proteins):
        self.index = index
        self.name = name
        self.proteins = dict_proteins


    def setNumOfSingleton(self, singletonsCount):
        self.numOfSingletons = self.numOfSingletons + singletonsCount

    def addToClusterList(self, clusterID):
        self.clusterList.add(clusterID)
# a = Strain(0, "aa", None)

    def getClusterList(self):
        return self.clusterList

    def getNumOfSingletons(self):
        return self.numOfSingletons

    def getNumOfGenes(self):
        return self.numOfGenes

    def increaseNumOfGenes(self):
        self.numOfGenes = self.numOfGenes + 1

    # def getStrainName(self):
    #     return self.name