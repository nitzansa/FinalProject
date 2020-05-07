import csv
import pandas as pd

import matplotlib.pyplot as plt
import statistics
import numpy as np
from CDHIT_Parser import CDHIT_Parser


class Artifact:

    global listOfClusters, minMemberLength, maxMemberLength, mean, std, strainsPerCluster, genesPerCluster, avgMembersPerCluster

    def __init__(self, path):
        self.listOfClusters = CDHIT_Parser(path)
        self.minMemberLength = {}
        self.maxMemberLength = {}
        self.mean = {}
        self.std = {}
        self.genesPerCluster = {}
        self.strainsPerCluster = {}
        self.avgMembersPerCluster = {}

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

    def downloadReport(self):

        with open('report one member.csv', mode='w') as report_one_member_csv: # TODO: change the file name
            report_one_member_writer = csv.writer(report_one_member_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            #the first row in file
            report_one_member_writer.writerow(['cluster num',
                                               'mean length',
                                               'std length',
                                               '# strains in each cluster',
                                               '# of members in each cluster',
                                               'mean of the per strain members',
                                               'min length',
                                               'max length',
                                               'min number of members per strains',
                                               'max number of members per strain',
                                               'flag'])
        report_one_member_csv.close()

        with open('report.csv', mode='w') as report_csv: # TODO: change the file name
            report_writer = csv.writer(report_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            #the first row in file
            report_writer.writerow(['cluster num',
                                    'mean length',
                                    'std length',
                                    '# strains in each cluster',
                                    '# of members in each cluster',
                                    'mean of the per strain members',
                                    'min length',
                                    'max length',
                                    'min number of members per strains',
                                    'max number of members per strain',
                                    'flag'])

            for cluster in self.listOfClusters.clusters.keys():
                dict_members = self.listOfClusters.getClusterMembers(cluster)
                # if for this cluster exist only one member
                if len(dict_members) < 2:
                    self.reportToClustersWithOneMember(cluster)
                else:
                    flag = 0
                    if len(self.strainsPerCluster[cluster]) == 1 and len(self.listOfClusters.getClusterMembers(cluster)) > 1:
                        flag = 2
                    if self.getMaxMembersPerStrainPerCluster(cluster) == 1:
                        flag = 3
                    report_writer.writerow([cluster,
                                            self.mean[cluster],
                                            self.std[cluster],
                                            len(self.strainsPerCluster[cluster]),
                                            len(self.listOfClusters.getClusterMembers(cluster)),
                                            self.avgMembersPerCluster[cluster],
                                            self.minMemberLength[cluster],
                                            self.maxMemberLength[cluster],
                                            self.getMinMembersPerStrainPerCluster(cluster),
                                            self.getMaxMembersPerStrainPerCluster(cluster),
                                            flag])
        report_csv.close()

    def reportToClustersWithOneMember(self, cluster):

        with open('report one member.csv', mode='a') as report_csv: # TODO: change the file name
            report_writer = csv.writer(report_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            dict_members = self.listOfClusters.getClusterMembers(cluster)
            member_length = 0
            for member in dict_members.values():
                member_length = member.getLength
            report_writer.writerow([cluster,
                                    member_length,
                                    '0',
                                    '1',
                                    '1',
                                    '1',
                                    member_length,
                                    member_length,
                                    '1',
                                    '1',
                                    '1'])
        report_csv.close()

    def classifyCluster(self):
        curr_path = "resources/report.csv"
        try:
            df = pd.read_csv(curr_path, usecols=['flag'])
            counts = df['flag'].value_counts().to_dict()
            counts[1] = len(self.getSingleClusters())
            return counts
        except IOError:
            print("could not open the file")
            return


    def downloadClassifyReport(self):
        with open('cluster reports/classify.csv', mode='w') as classify__csv:  # TODO: change the file name
            classify_writer = csv.writer(classify__csv, delimiter=',', quotechar='"',
                                                  quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            classify_writer.writerow(['class type', 'count'])
            classes_dict = self.classifyCluster()
            for key, val in classes_dict.items():
                classify_writer.writerow([key, val])
            classify__csv.close()

    def downloadStrainSingletonsReport(self):
        singleton_strains = []
        singletons = self.getSingleClusters()
        countOfSingletons = len(singletons)
        for cluster in singletons:
            members = self.listOfClusters.getClusterMembers(cluster)
            for member in members.values():
                singleton_strains.append(member.getStrainInd)

        df = pd.DataFrame(singleton_strains, columns=['strain index'])
        counts = df['strain index'].value_counts().to_dict()

        with open('resources/singletonStrain.csv', mode='w') as singletons_starin_csv:  # TODO: change the file name
            singleton_strains_writer = csv.writer(singletons_starin_csv, delimiter=',', quotechar='"',
                                                  quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            singleton_strains_writer.writerow(['strain ind', 'count of singletons', '% of singletons'])
            for key, val in counts.items():
                singleton_strains_writer.writerow([key, val, (val / countOfSingletons) * 100])

            singletons_starin_csv.close()

# CD_output = CDHIT_Parser("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
# a = Artifact("resources/23cluster")
a.variableLength()
a.getGenesPerCluster()
a.getStrainsPerCluster()
# a.getMinStrainsPerCluster(2)
# print(a.getSingleClusters())
# print(len(a.getSingleClusters()))
# a.getStrainsPerCluster()
a.calcAverageMemberPerCluster()
# a.clustersPerCountOfStrains()
# print(a.strainsPerCluster)
# a.downloadReport()
# a.classifyCluster("resources/report")
# a.downloadClassifyReport()
a.downloadStrainSingletonsReport()