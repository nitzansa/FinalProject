from Protein import Protein
from Member import Member
from Artifacts import Artifact
from Strain import Strain



class OutlierDetection:
    global artifacts, flagPerCluster, listOfStrains, most_common_length_dict, listOfStrainPerCluster, \
        listOfClass4, list_of_clusters_to_check_outlier

    def __init__(self, artifacts):
        self.artifacts = artifacts
        self.flagPerCluster = {}
        self.most_common_length_dict = {}
        self.listOfStrains = {}
        self.listOfStrainPerCluster = {}
        self.listOfClass4 = []
        self.list_of_clusters_to_check_outlier = []

        self.listOfClass4 = self.artifacts.listOfClass4
        self.flagPerCluster = self.artifacts.flagPerCluster
        self.most_common_length_dict = self.artifacts.most_common_length_dict
        self.listOfStrains = self.artifacts.listOfStrains
        self.listOfStrainPerCluster = self.artifacts.getStrainsPerCluster()

    def detectOutlier(self):
        self.countOfNonOutlier()
        # list_of_clusters_to_check_outlier = self.countOfNonOutlier()
        for cluster in self.list_of_clusters_to_check_outlier:
            self.checkOutliersInClusters(cluster)



    def countOfNonOutlier(self):
        # list_of_clusters_to_check_outlier = []
        nonOutlier = 0
        atLeastStrainsInCLuster = 50
        most_common_percent = 80
        for cluster in self.listOfClass4:
            if self.most_common_length_dict[cluster]['%_1'] > most_common_percent and \
                    len(self.listOfStrainPerCluster[cluster]) >= atLeastStrainsInCLuster:
                nonOutlier = nonOutlier + 1
                self.list_of_clusters_to_check_outlier.append(cluster)
        # print(self.list_of_clusters_to_check_outlier)
        print(nonOutlier)
        return self.list_of_clusters_to_check_outlier

    def checkOutliersInClusters(self, cluster):
        threshold = 0.5
        common_length = self.most_common_length_dict[cluster]['length_1']
        additionalRange = threshold * common_length

        clusterMembers = self.artifacts.listOfClusters.getClusterMembers(cluster)
        for member in clusterMembers.values():
            if member.getLength != common_length:
                if member.getLength < common_length - additionalRange or member.length > common_length + additionalRange:
                    if member.getStrainInd in self.listOfStrains.keys():
                        self.listOfStrains.get(member.getStrainInd).increaseNumOfOutlierGenes()
            #         outlier!!!!!!
            #        add outlier to the strain of this member
