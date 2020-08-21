from Protein import Protein
from Member import Member
from Artifacts import Artifact
from Strain import Strain



class OutlierDetection:
    global artifacts, flagPerCluster, listOfStrains, most_common_length_dict, listOfStrainPerCluster, \
        list_of_clusters_outlier_class4, list_of_clusters_outlier_class0
    # Outlier detection check!!!!!!!!!
    def __init__(self, artifacts):
        self.artifacts = artifacts
        # self.flagPerCluster = {}
        self.most_common_length_dict = {}
        # self.listOfStrains = {}
        self.listOfStrainPerCluster = {}
        # self.listOfClass4 = []
        # self.list_of_clusters_to_check_outlier = []
        self.list_of_clusters_outlier_class4 = []
        self.list_of_clusters_outlier_class0 = []

        # self.listOfClass4 = self.artifacts.listOfClass4
        # self.flagPerCluster = self.artifacts.flagPerCluster
        # self.most_common_length_dict = self.artifacts.most_common_length_dict
        # self.listOfStrains = self.artifacts.listOfStrains
        self.listOfStrainPerCluster = self.artifacts.getStrainsPerCluster()

    def detectOutlier(self):
        self.createListOfOutliers()
        # list_of_clusters_to_check_outlier = self.countOfNonOutlier()
        for cluster in self.list_of_clusters_outlier_class4:
            self.checkOutliersInClusters_length(cluster, 4)

        for cluster in self.list_of_clusters_outlier_class0:
            self.checkOutliersInClusters_length(cluster, 0)





    def createListOfOutliers(self):
        # list_of_clusters_to_check_outlier = []
        nonOutlier = 0
        atLeastStrainsInCLuster = 50
        most_common_percent = 80
        for cluster in self.artifacts.listOfClass4:
            if self.artifacts.most_common_length_dict[cluster]['%_1'] > most_common_percent and \
                    len(self.listOfStrainPerCluster[cluster]) >= atLeastStrainsInCLuster:
                nonOutlier = nonOutlier + 1
                self.list_of_clusters_outlier_class4.append(cluster)
        # print(self.list_of_clusters_to_check_outlier)
        # print(nonOutlier)

        for cluster in self.artifacts.listOfClass0:
            if self.artifacts.most_common_length_dict[cluster]['%_1'] > most_common_percent and \
                    len(self.listOfStrainPerCluster[cluster]) >= atLeastStrainsInCLuster:
                nonOutlier = nonOutlier + 1
                self.list_of_clusters_outlier_class0.append(cluster)

    def checkOutliersInClusters_length(self, cluster, classID):
        threshold = 0.3
        common_length = self.artifacts.most_common_length_dict[cluster]['length_1']
        additionalRange = threshold * common_length

        clusterMembers = self.artifacts.listOfClusters.getClusterMembers(cluster)
        for member in clusterMembers.values():
            if member.getLength != common_length:
                if member.getLength < common_length - additionalRange or member.length > common_length + additionalRange:
                    if member.getStrainInd in self.artifacts.listOfStrains.keys():
                        if classID == 4:
                            self.artifacts.listOfStrains.get(member.getStrainInd).increaseNumOfOutlierGenes_class4_length()

                        if classID == 0:
                            self.artifacts.listOfStrains.get(member.getStrainInd).increaseNumOfOutlierGenes_class0_length()

            #         outlier!!!!!!
            #        add outlier to the strain of this member
