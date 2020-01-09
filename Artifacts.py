import matplotlib.pyplot as plt
import statistics
import numpy as np
from CDHIT_Parser import CDHIT_Parser


class Artifact:

    global listOfClusters
    global mean
    global std
    global strainsPerCluster
    global genesPerCluster

    def __init__(self, path):
        self.listOfClusters = CDHIT_Parser(path)
        self.mean = {}
        self.std = {}
        self.genesPerCluster = {}
        self.strainsPerCluster = {}

    def getGenesPerCluster(self):
        for cluster in self.listOfClusters.clusters.keys():
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            self.genesPerCluster[cluster] = len(dict_members)

        return self.genesPerCluster

    def getStrainsPerCluster(self):
        for cluster in self.listOfClusters.clusters.keys():
            listOfDifferentStrains = []
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            for member in dict_members.values():
                listOfDifferentStrains.append(member.getStrainInd)
            x = np.array(listOfDifferentStrains)
            unique, counts = np.unique(x, return_counts=True)
            # y = np.bincount(x)
            # ii = np.nonzero(y)[0]
            self.strainsPerCluster[cluster] = np.asarray((unique, counts)).T
            # self.strainsPerCluster[cluster] = len(listOfDifferentStrains)

        return self.strainsPerCluster

    def getSingleClusters(self):
        singletons = []
        for cluster in self.listOfClusters.clusters.keys():
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            if len(dict_members) == 1:
                singletons.append(cluster)

        return singletons

    # graph for STD
    def variableLength(self):
        for cluster in self.listOfClusters.clusters.keys():
            data = []
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            #if for this cluster exist only one member
            if len(dict_members) < 2:
                continue
            for member in dict_members.values():
                data.append(member.getLength)
            self.mean[cluster] = statistics.mean(data)
            self.std[cluster] = statistics.stdev(data)
        keys = [] #cluster
        values = [] #std for each cluster
        for key, value in self.std.items():
            keys.append(float(key))
            values.append(value)
        # plotting the points
        plt.plot(keys, values, color='grey', linestyle='dashed', linewidth=1,
                 marker='o', markerfacecolor='blue', markersize=5)
        plt.ylim(1, max(values))
        plt.xlim(1, len(keys))
        plt.xlabel('cluster')
        plt.ylabel('std')
        plt.title('Variable length of proteins inside cluster')
        plt.show()



# CD_output = CDHIT_Parser("cd_test")
a = Artifact("clusters_output")
# a.variableLength()
# print(a.getGenesPerCluster())
# print(a.getStrainsPerCluster())
# print(a.getSingleClusters())
a.getStrainsPerCluster()