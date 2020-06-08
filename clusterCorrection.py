from Protein import Protein
from Member import Member
from Artifacts import Artifact
from Strain import Strain
from ProteinFilesManager import ProteinFilesManager




class clusterCorrection:
    global artifacts, listOfClusters, indexOfNextCluster, listOfStrains


    def __init__(self, artifacts):
        self.artifacts = artifacts
        self.listOfClusters = self.artifacts.getClusterList()
        self.indexOfNextCluster = len(self.listOfClusters) + 1
        self.listOfStrains = self.artifacts.listOfStrains

    def getTheBiggestParalog(self, clusterInd):
        strainInCluster = {}
        clusterMembers = self.listOfClusters[clusterInd]

        for member in clusterMembers:
            strainInd = clusterMembers[member].strain_index
            if strainInd in strainInCluster.keys():
                strainInCluster[strainInd] = strainInCluster[strainInd] + 1
            else:
                strainInCluster[strainInd] = 1

        biggestParalog = max(strainInCluster, key=strainInCluster.get)
        return biggestParalog

    def getMembersFromStrain(self,strain_index, cluster_index):
        membersFromStrain = {}
        clusterMembers = self.listOfClusters[cluster_index]
        for member in clusterMembers:
            if clusterMembers[member].strain_index == strain_index:
                membersFromStrain[member] = clusterMembers[member]
        return membersFromStrain

    def getMembersToReClustering(self,strain_index, cluster_index):
        # insert the strain of the paralog, and the cluster index to split. -> and return the rest of the members in the cluster
        # to split to news.

        membersFromOtherStrains = {}
        clusterMembers = self.listOfClusters[cluster_index]
        for member in clusterMembers:
            if clusterMembers[member].strain_index != strain_index:
                membersFromOtherStrains[member] = clusterMembers[member]
        return membersFromOtherStrains



    def createNewCluster(self,paralog_dict, oldClusterID):


        for paralog_member in paralog_dict:
            value = {}  # dict of members
            paralog_member.represntative = True
            paralog_member.identity = 100.0

            value[0] = paralog_member
            self.listOfClusters[self.indexOfNextCluster] = value
            if paralog_member.strain_index in self.listOfStrains.keys():
                self.listOfStrains.get(paralog_member.strain_index).deleteFromClusterList(oldClusterID)
                self.listOfStrains.get(paralog_member.strain_index).addToClusterList(self.indexOfNextCluster)
            value[0] = paralog_member
            self.listOfClusters[self.indexOfNextCluster] = value
            self.indexOfNextCluster = self.indexOfNextCluster + 1



        #     new_member = Member(strain_index, protein_index, locus_tag, length, identity, represntative)
        #     if strain_index in self.dict_strains.keys():
        #         self.dict_strains.get(strain_index).addToClusterList(key)
        #
        #     value[int(line[0].split("\t")[0])] = new_member
        #     line = file.readline().split(" ")
        # self.dict_clusters[key] = value
