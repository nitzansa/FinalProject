import matplotlib.pyplot as plt
import statistics

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
            listOfDifferentStrains = set()
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            for member in dict_members.values():
                listOfDifferentStrains.add(member.getStrainInd)
            self.strainsPerCluster[cluster] = len(listOfDifferentStrains)

        return self.strainsPerCluster

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
        keys = []
        values = []
        # for key, value in self.mean.items():
        #     keys.append(float(key))
        #     values.append(value)
        for key, value in self.std.items():
            keys.append(float(key))
            values.append(value)
        plt.plot(keys, values)
        # naming the x axis
        plt.xlabel('cluster')
        # naming the y axis
        plt.ylabel('mean')
        # giving a title to my graph
        plt.title('Variable length of proteins inside cluster')
        plt.show()





# CD_output = CDHIT_Parser("cd_test")
a = Artifact("clusters_output")
# a.variableLength()
print(a.getGenesPerCluster())
print(a.getStrainsPerCluster())

