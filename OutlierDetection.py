from Protein import Protein
from Member import Member
from Artifacts import Artifact
from Strain import Strain



class OutlierDetection:
    global artifacts, flagPerCluster, listOfStrains, most_common_length_dict, listOfStrainPerCluster

    def __init__(self, artifacts):
        self.artifacts = artifacts
        self.flagPerCluster = {}
        self.most_common_length_dict = {}
        self.listOfStrains = {}
        self.listOfStrainPerCluster = {}
        self.flagPerCluster = self.artifacts.flagPerCluster
        self.most_common_length_dict = self.artifacts.most_common_length_dict
        self.listOfStrains = self.artifacts.listOfStrains
        self.listOfStrainPerCluster = self.artifacts.getStrainsPerCluster()



    def countOfNonOutlier(self):
        nonOutlier = 0
        distance = 30
        atLeastStrainsInCLuster = 50
        most_common_percent = 80
        for cluster in self.flagPerCluster.keys():
            if self.flagPerCluster[cluster] == 4 and self.most_common_length_dict[cluster]['%_1'] > most_common_percent and \
                    len(self.listOfStrainPerCluster[cluster]) >= atLeastStrainsInCLuster:
                nonOutlier = nonOutlier + 1

        print(nonOutlier)
