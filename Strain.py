import Bio
from Bio.Seq import Seq
from Protein import Protein


class Strain:
    # strain!!

    index = 0
    name = ""
    proteins = {}
    numOfCoreGenes = 0
    global clusterList, numOfOutliers_class4_length, numOfSingletons, numOfGenes, numOfOutliers_class0_length, \
        numOfOutliers_class0_members
    #clusterList = {}

    def __init__(self, index, name, dict_proteins):
        self.index = index
        self.name = name
        self.proteins = dict_proteins
        self.numOfGenes = 0
        self.numOfCoreGenes = 0
        self.numOfSingletons = 0
        self.numOfOutliers_class4_length = 0
        self.numOfOutliers_class0_length = 0
        self.numOfOutliers_class0_members = 0
        self.clusterList = set()


    def setNumOfSingleton(self, singletonsCount):
        self.numOfSingletons = self.numOfSingletons + singletonsCount

    def addToClusterList(self, clusterID):
        self.clusterList.add(clusterID)
# a = Strain(0, "aa", None)

    def deleteFromClusterList(self, clusterID):
        self.clusterList.remove(clusterID)

    def getClusterList(self):
        return self.clusterList

    def getNumOfSingletons(self):
        return self.numOfSingletons

    def getNumOfGenes(self):
        return self.numOfGenes

    # def increaseNumOfGenes(self):
    #     self.numOfGenes = self.numOfGenes + 1

    def increaseNumOfCoreGenes(self):
        self.numOfCoreGenes = self.numOfCoreGenes + 1

    def getProteins (self):
        return self.proteins

    def increaseNumOfOutlierGenes_class4_length(self):
        self.numOfOutliers_class4_length = self.numOfOutliers_class4_length + 1

    def increaseNumOfOutlierGenes_class0_length(self):
            self.numOfOutliers_class0_length = self.numOfOutliers_class0_length + 1

    def increaseNumOfOutlierGenes_class0_members(self):
            self.numOfOutliers_class0_members = self.numOfOutliers_class0_members + 1

