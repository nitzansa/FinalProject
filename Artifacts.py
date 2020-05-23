import matplotlib.pyplot as plt
import statistics
import numpy as np
from CDHIT_Parser import CDHIT_Parser


class Artifact:

    global listOfClusters, minMemberLength, maxMemberLength, mean, std, strainsPerCluster, genesPerCluster, avgMembersPerCluster, listOfStrains

    def __init__(self, path, strains):
        self.listOfClusters = CDHIT_Parser(path, strains)
        self.listOfStrains = strains
        self.minMemberLength = {}
        self.maxMemberLength = {}
        self.mean = {}
        self.std = {}
        self.genesPerCluster = {}
        self.strainsPerCluster = {}
        self.avgMembersPerCluster = {}
        self.variableLength()
        self.getGenesPerCluster()
        self.getSingleClusters()
        self.getStrainsPerCluster()
        self.calcAverageMemberPerCluster()

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

        import csv
        with open('reports\\forHist.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            for i in range(1, len(counts)):
                writer.writerow([unique[i], counts[i]])

        # # plotting the points
        # plt.plot(unique, counts)
        #
        # # naming the x axis
        # plt.xlabel('x - axis')
        # # naming the y axis
        # plt.ylabel('y - axis')
        #
        # # giving a title to my graph
        # plt.title('My first graph!')
        #
        # # function to show the plot
        # plt.show()

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
        #graph for STD
        # self.showGraphPlotting("Variable length of proteins inside cluster", "std", "cluster", self.std)

        #graph for mean
        # self.showGraphPlotting("Mean length of proteins inside cluster", "mean", "cluster", self.mean)

    """
    This function export to a csv file the average of cluster members for each cluster. 
    """
    def calcAverageMemberPerCluster(self):

        for cluster in self.genesPerCluster.keys():
            geneCount = self.genesPerCluster.get(cluster)
            strainCount = len(self.strainsPerCluster[cluster])
            self.avgMembersPerCluster[cluster] = geneCount / strainCount

        # import csv
        # with open('reports\\averageMembers.csv', 'w', newline='') as file:
        #     writer = csv.writer(file)
        #     for cluster in self.avgMembersPerCluster.keys():
        #         writer.writerow([cluster, self.avgMembersPerCluster.get(cluster)])

    """
    Show graph plotting of selected statistic.
    name- the name of the plot
    xlabel- the of x-axis
    ylabel- the of y-axis
    data- a dictionary with values for x-axis and y-axis
    """
    def showGraphPlotting(self, name, xlabel, ylabel, data):
        keys = []
        values = []

        for key, value in data.items():
            keys.append(float(key))
            values.append(value)

        # plotting the points
        plt.plot(values, keys, color='grey', linestyle='dashed', linewidth=1, marker='o', markerfacecolor='blue', markersize=5)
        plt.xlim(0, max(values))
        plt.ylim(0, max(keys))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(name)
        plt.show()

    def isCoreCluster(self, cluster):

        counter = 0
        for i in self.strainsPerCluster[cluster]:
            if i[1] == 1:
              counter = counter + 1

        if counter / len(self.strainsPerCluster[cluster]) >= 0.9:
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
                strainListPerCluster = np.array(self.getStrainsPerCluster()[cluster])
                for i in range(0, len(strainListPerCluster.__array__())):
                    self.listOfStrains.get(strainListPerCluster.__array__()[i][0]).increaseNumOfCoreGenes()





# CD_output = CDHIT_Parser("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
# # a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
# a = Artifact("resources/23cluster")
# a.variableLength()
# a.getGenesPerCluster()
# a.getStrainsPerCluster()
# a.getMinStrainsPerCluster(2)
# print(a.getSingleClusters())
# print(len(a.getSingleClusters()))
# a.getStrainsPerCluster()
# a.calcAverageMemberPerCluster()
# a.clustersPerCountOfStrains()
# print(a.strainsPerCluster)
# a.downloadReport()
# a.classifyCluster("resources/report")
# a.downloadClassifyReport()
# a.downloadStrainSingletonsReport()
# print(a.isCoreCluster(6))
# print(a.avgIdentity(6))
# print(a.stdIdentity(6))
# print(a.PercentOfMembersWithA_Score(6))
# print(a.PercentOfMembersWithB_Score(6))
# print(a.PercentOfMembersWithC_Score(6))