from Protein import Protein
from Member import Member
from Artifacts import Artifact
from Strain import Strain
from ProteinFilesManager import ProteinFilesManager




class clusterCorrection:
    global artifacts, listOfClusters, indexOfNextCluster


    def __init__(self, artifacts):
        self.artifacts = artifacts
        self.listOfClusters = self.artifacts.getClusterList()
        self.indexOfNextCluster = len(self.listOfClusters)
        print(self.indexOfNextCluster)