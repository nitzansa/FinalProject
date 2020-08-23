import matplotlib.pyplot as plt
import statistics
import numpy as np
from CDHIT_Parser import CDHIT_Parser
import collections
import pandas as pd
from collections import Counter


class Artifact:
    # artifacts!!!!!
    global listOfClusters, minMemberLength, maxMemberLength, mean, std, strainsPerCluster, genesPerCluster, \
        avgMembersPerCluster, listOfStrains, flagPerCluster, most_common_length_dict, listOfClass0, listOfClass2, \
        listOfClass3, listOfClass4, listOfClass5, singletons, sum_of_core_clusters

    def __init__(self, path, strains):
        self.listOfClusters = CDHIT_Parser(path, strains)
        self.listOfStrains = strains
        self.sum_of_core_clusters = 0
        self.minMemberLength = {}
        self.maxMemberLength = {}
        self.mean = {}
        self.std = {}
        self.genesPerCluster = {}
        self.strainsPerCluster = {}
        self.avgMembersPerCluster = {}
        self.flagPerCluster = {}
        self.listOfClass0 = []
        self.listOfClass2 = []
        self.listOfClass3 = []
        self.listOfClass4 = []
        self.listOfClass5 = []
        self.most_common_length_dict = self.calculatingLengthDistributionOfEachCluster()
        self.variableLength()
        self.getGenesPerCluster()
        self.singletons = self.getSingleClusters()
        self.getStrainsPerCluster()
        self.calcAverageMemberPerCluster()
        self.calcFlagPerCluster()


    """
    This is a part of the first steps about the statistic.
    To identify the amount of members in each cluster, this function return a dictionary with 
    key- cluster ID and value- the number of proteins members of this cluster  
    """
    def getGenesPerCluster(self):
        for cluster in self.listOfClusters.clusters.keys():
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            self.genesPerCluster[cluster] = len(dict_members)

        return self.genesPerCluster

    """
    This is a part of the first steps about the statistic.
    To identify the amount of different strains in each cluster, this function return a dictionary with 
    key- cluster ID and value- list with size of 2, the first cell is the strain ID and the second cell is the frequency 
    of this strain in this cluster 
    """
    def getStrainsPerCluster(self):
        for cluster in self.listOfClusters.clusters.keys():
            listOfDifferentStrains = []
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            for member in dict_members.values():
                listOfDifferentStrains.append(member.getStrainInd)
            x = np.array(listOfDifferentStrains)
            unique, counts = np.unique(x, return_counts=True)
            self.strainsPerCluster[cluster] = np.asarray((unique, counts)).T

        return self.strainsPerCluster

    def getMinMembersPerStrainPerCluster(self, cluster):

        numOfMemberPerStrain = []
        for x in self.strainsPerCluster[cluster]:
            numOfMemberPerStrain.append(x[1])

        return min(numOfMemberPerStrain)

    def getMaxMembersPerStrainPerCluster(self, cluster):

        numOfMemberPerStrain = []
        for x in self.strainsPerCluster[cluster]:
            numOfMemberPerStrain.append(x[1])

        return max(numOfMemberPerStrain)

    def getSTDMembersPerStrainPerCluster(self, cluster):

        numOfMemberPerStrain = []
        for x in self.strainsPerCluster[cluster]:
            numOfMemberPerStrain.append(x[1])
            # print(x)
        # print(numOfMemberPerStrain)
        # print(len(numOfMemberPerStrain))
        # a = [1,2]
        # print(statistics.stdev(a))
        # print(statistics.stdev(1,1))
        if len(numOfMemberPerStrain) < 2:
            return 0

        return statistics.stdev(numOfMemberPerStrain)

    """
    This is a part of the first steps about the statistic.
    To identify all the clusters with one member, this function return a list with all clusters ID's that contains a 
    single member 
    """
    def getSingleClusters(self):
        singletons = []
        for cluster in self.listOfClusters.clusters.keys():
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            if len(dict_members) == 1:
                singletons.append(cluster)

        return singletons

    """
    This function export to a csv file the amount of clusters who contains the same amount of different strains.
    """
    def clustersPerCountOfStrains(self):
        frequency = []
        for strainFreq in self.strainsPerCluster.values():
            frequency.append(len(strainFreq))
        x = np.array(frequency)
        #unique- count of different strains. counts- frequencies
        unique, counts = np.unique(x, return_counts=True)

    """
    statistic about the first artifact - Variable Length.
    Show graph plotting to the std about the length of proteins inside cluster.
    Show graph plotting to the mean about the length of proteins inside cluster.
    This graphs for further study.
    """
    def variableLength(self):
        for cluster in self.listOfClusters.clusters.keys():
            data = []
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            #if for this cluster exist only one member
            if len(dict_members) < 2:
                continue
            for member in dict_members.values():
                data.append(member.getLength)
            self.minMemberLength[cluster] = min(data)
            self.maxMemberLength[cluster] = max(data)
            self.mean[cluster] = statistics.mean(data)
            self.std[cluster] = statistics.stdev(data)
        self.mean = {k: v for k, v in sorted(self.mean.items(), key=lambda item: item[1])}
        self.std = {k: v for k, v in sorted(self.std.items(), key=lambda item: item[1])}

    """
    This function export to a csv file the average of cluster members for each cluster. 
    """
    def calcAverageMemberPerCluster(self):

        for cluster in self.genesPerCluster.keys():
            geneCount = self.genesPerCluster.get(cluster)
            strainCount = len(self.strainsPerCluster[cluster])
            self.avgMembersPerCluster[cluster] = geneCount / strainCount

    def isCoreCluster(self, cluster):

        counter = 0
        for i in self.strainsPerCluster[cluster]:
            if i[1] == 1:
              counter = counter + 1

        if counter / len(self.listOfStrains) >= 0.9:
            return True

        return False

    def avgIdentity(self, cluster):

        dict_member = self.listOfClusters.getClusterMembers(cluster)
        list_identity = []
        for member in dict_member.values():
            list_identity.append(member.getIdentity)

        return statistics.mean(list_identity)

    def stdIdentity(self, cluster):

        dict_member = self.listOfClusters.getClusterMembers(cluster)
        list_identity = []
        for member in dict_member.values():
            list_identity.append(member.getIdentity)

        return statistics.stdev(list_identity)

    #% of members with %score 70-80
    def PercentOfMembersWithC_Score(self, cluster):

        counter = 0
        dict_member = self.listOfClusters.getClusterMembers(cluster)
        for member in dict_member.values():
            if 70 <= member.getIdentity < 80:
                counter = counter + 1

        return (counter / len(dict_member)) * 100

    # % of members with %score 80-90
    def PercentOfMembersWithB_Score(self, cluster):

        counter = 0
        dict_member = self.listOfClusters.getClusterMembers(cluster)
        for member in dict_member.values():
            if 80 <= member.getIdentity < 90:
                counter = counter + 1

        return (counter / len(dict_member)) * 100

    # % of members with %score 90-100
    def PercentOfMembersWithA_Score(self, cluster):

        counter = 0
        dict_member = self.listOfClusters.getClusterMembers(cluster)
        for member in dict_member.values():
            if 90 <= member.getIdentity <= 100:
                counter = counter + 1

        return (counter / len(dict_member)) * 100


    def updateNumOfCoreGeneInStrains(self):
        # for cluster in self.listOfClusters.clusters.keys():
        #     if self.isCoreCluster(cluster):
        #         for strain in self.getStrainsPerCluster():
        for cluster in self.listOfClusters.clusters.keys():
            if self.isCoreCluster(cluster):
                strainListPerCluster = self.strainsPerCluster[cluster]
                for strain in strainListPerCluster:
                    self.listOfStrains.get(strain[0]).increaseNumOfCoreGenes()
                    # print(self.listOfStrains.get(strain[0]).numOfCoreGenes)

    def getNeighbours(self, member):
        neighbours_dict = {} # key- StrainInd/ProteinInd, val-locus_tag
        for i in range(5):
            i = i + 1
            if self.listOfStrains.get(member.getStrainInd).getProteins()['locus_tag'][member.getProteinInd + i] is not None:
                neighbours_dict[str(member.getStrainInd) + '/' + str(member.getProteinInd + i)] = \
                    self.listOfStrains.get(member.getStrainInd).getProteins()['locus_tag'][member.getProteinInd + i]
            if self.listOfStrains.get(member.getStrainInd).getProteins()['locus_tag'][member.getProteinInd - i] is not None:
                neighbours_dict[str(member.getStrainInd) + '/' + str(member.getProteinInd - i)] = \
                    self.listOfStrains.get(member.getStrainInd).getProteins()['locus_tag'][member.getProteinInd - i]

        return neighbours_dict


    def getNeighborsClusters(self, _member):
        neighbors_clusters_dict = {}
        neighbors_dict = self.getNeighbours(_member)

        for neighbour in neighbors_dict.keys():
            neighbors_clusters_dict[neighbour] = set()

        for cluster in self.listOfClusters.clusters.keys():
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            for member in dict_members.values():
                possible_key = str(member.getStrainInd) + '/' + str(member.getProteinInd)
                if possible_key in neighbors_dict:
                    neighbors_clusters_dict[possible_key].add(cluster)
        return neighbors_clusters_dict

    def getClusterList(self):
        return self.listOfClusters.clusters

    def calcFlagPerCluster(self):
        for cluster in self.listOfClusters.clusters.keys():
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            # if for this cluster exist only one member
            if len(dict_members) >= 2:
                flag = 0
                if len(self.strainsPerCluster[cluster]) == 1 and len(
                        self.listOfClusters.getClusterMembers(cluster)) > 1:
                    flag = 2
                if self.getMaxMembersPerStrainPerCluster(cluster) == 1:
                    flag = 3
                if flag == 3 and 30 <= self.most_common_length_dict[cluster]['%_1'] < 100:
                    flag = 4
                if flag == 3 and self.most_common_length_dict[cluster]['%_1'] < 30:
                    flag = 5
                # self.flagPerCluster[cluster] = flag
                # print(cluster)
                if flag == 0:
                    self.listOfClass0.append(cluster)
                if flag == 2:
                    self.listOfClass2.append(cluster)
                if flag == 3:
                    self.listOfClass3.append(cluster)
                if flag == 4:
                    self.listOfClass4.append(cluster)
                if flag == 5:
                    self.listOfClass5.append(cluster)

   # for cluster in self.artifacts.listOfClusters.clusters.keys():

    # dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster)
    # # if for this cluster exist only one member
    # if len(dict_members) < 2:
    #     True
    #     # self.reportToClustersWithOneMember(cluster)
    # # elif len(self.artifacts.strainsPerCluster[cluster]) == 1 and len(
    # #             self.artifacts.listOfClusters.getClusterMembers(cluster)) > 1:
    # # #     flag ==2
    # #     True
    # else:
    #     flag = 0
    #     # write clusters from class 2 in the report, to check if it can be deleted?
    #     if len(self.artifacts.strainsPerCluster[cluster]) == 1 and len(
    #             self.artifacts.listOfClusters.getClusterMembers(cluster)) > 1:
    #         flag = 2
    #     if self.artifacts.getMaxMembersPerStrainPerCluster(cluster) == 1:
    #         flag = 3
    #     if flag == 3 and 30 <= self.most_common_length_dict[cluster]['%_1'] < 100:
    #         flag = 4
    #     if flag == 3 and self.most_common_length_dict[cluster]['%_1'] < 30:
    #         flag = 5

    def calculatingLengthDistributionOfEachCluster(self): #top 3
        most_common_length_dict = {}
        for cluster in self.listOfClusters.clusters:
            length_freq = []
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            if len(dict_members) > 1:
                for member in dict_members.values():
                    length_freq.append(member.getLength)
                df = pd.DataFrame(length_freq, columns=['strain index'])
                counts = df['strain index'].value_counts().to_dict()
                c = Counter(counts)
                top3 = c.most_common(3)
                if len(top3) == 1:
                    most_common_length_dict[cluster] = {'length_1': top3[0][0],
                                                        '%_1': (top3[0][1] / len(dict_members)) * 100,
                                                        'length_2': 0, '%_2': 0,
                                                        'length_3': 0, '%_3': 0}
                elif len(top3) == 2:
                    most_common_length_dict[cluster] = {'length_1': top3[0][0],
                                                        '%_1': (top3[0][1] / len(dict_members)) * 100,
                                                        'length_2': top3[1][0],
                                                        '%_2': (top3[1][1] / len(dict_members)) * 100,
                                                        'length_3': 0, '%_3': 0}
                elif len(top3) == 3:
                    most_common_length_dict[cluster] = {'length_1': top3[0][0],
                                                        '%_1': (top3[0][1] / len(dict_members)) * 100,
                                                        'length_2': top3[1][0],
                                                        '%_2': (top3[1][1] / len(dict_members)) * 100,
                                                        'length_3': top3[2][0],
                                                        '%_3': (top3[2][1] / len(dict_members)) * 100}
        return most_common_length_dict

    def updateSingletonsStrainCount(self):
        singleton_strains = []
        clustersFromClass1 = self.singletons
        # countOfSingletonsClass2 = len(clustersFromClass2)
        # print(singletons)
        # print(clustersFromClass2)
        # print(merge)

        for cluster in clustersFromClass1:
            members = self.listOfClusters.getClusterMembers(cluster)
            for member in members.values():
                if member.getStrainInd in self.listOfStrains.keys():
                    self.listOfStrains.get(member.getStrainInd).setNumOfSingleton(1)
                break

    def increase_sum_of_core_clusters(self):
        self.sum_of_core_clusters = self.sum_of_core_clusters + 1