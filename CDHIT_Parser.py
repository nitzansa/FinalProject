import csv

from Member import Member
from Strain import Strain


class CDHIT_Parser:

    global dict_clusters
    #  parser chaeck!!!
    def __init__(self, path, strains):

        self.dict_clusters = {}
        self.dict_strains = strains
        self.readFile(path)

    # a getter function
    @property
    def clusters(self):
        return self.dict_clusters

    def dict_clusters_setter(self, key, val):
       self.dict_clusters[key] = val

    def removeCluster(self, cluster_index):
       del self.dict_clusters[cluster_index]

    def readFile(self, path):
        try:
            file = open(path + ".txt", "r")
            with file:
                line = file.readline().split(" ")
                while line.__len__() > 1:
                    if line[0][0] == '>':  # new cluster
                        key = int(line[1])
                        value = {}  # dict of members
                        line = file.readline().split(" ")
                        while line.__len__() > 2:
                            strain_index = int(line[1].split("|")[0].replace(">", ""))
                            protein_index = int(line[1].split("|")[1].replace("...", ""))
                            length = int(line[0].split("\t")[1].replace("aa,", ""))
                            locus_tag = self.dict_strains.get(strain_index).getProteins()['locus_tag'][protein_index]
                            if line[2] == '*\n' or line[2] == '*':
                                represntative = True
                                identity = 100.0
                            else:
                                identity = float(line[3].split("/")[1].replace("%", "").replace("\n", ""))
                                represntative = False
                            new_member = Member(strain_index, protein_index, locus_tag, length, identity, represntative)
                            if strain_index in self.dict_strains.keys():
                                self.dict_strains.get(strain_index).addToClusterList(key)

                            value[int(line[0].split("\t")[0])] = new_member
                            line = file.readline().split(" ")
                        self.dict_clusters[key] = value
            file.close()
        except IOError:
            print("could not open the file")

    def getClusterMembers(self, cluster_index):
        return self.dict_clusters[cluster_index]

    # def downloadClusterFrenqLength(self, cluster_index):
    #     clusterFreqLength = {}
    #     countOfMembers = 0
    #     dict_members = self.getClusterMembers(int(cluster_index))
    #     for member in dict_members.values():
    #         if member.getLength in clusterFreqLength.keys():
    #             clusterFreqLength[member.getLength] = clusterFreqLength[member.getLength] + 1
    #         else:
    #             clusterFreqLength[member.getLength] = 1
    #         countOfMembers = countOfMembers + 1
    #
    #     with open('cluster reports/frequency length cluster_' + cluster_index + '.csv',
    #               mode='w') as cluster_freq__csv:  # TODO: change the file name
    #         cluster_freq_length_writer = csv.writer(cluster_freq__csv, delimiter=',', quotechar='"',
    #                                                 quoting=csv.QUOTE_MINIMAL)
    #
    #         # the first row in file
    #         cluster_freq_length_writer.writerow(['length', 'count of members', '% of members'])
    #         for length in clusterFreqLength.keys():
    #             cluster_freq_length_writer.writerow(
    #                 [length, clusterFreqLength[length], (clusterFreqLength[length] / countOfMembers) * 100])
    #
    #         cluster_freq__csv.close()


