from Protein import Protein
from Member import Member
from Artifacts import Artifact
from Strain import Strain
from ProteinFilesManager import ProteinFilesManager


class clusterCorrection:
    global artifacts, cluster_index, cluster_members, indexOfNextCluster, listOfStrains, list_of_new_clusters

    def __init__(self, artifacts, cluster_index):
        self.artifacts = artifacts
        self.cluster_index = cluster_index
        self.cluster_members = self.artifacts.listOfClusters.getClusterMembers(cluster_index)
        self.indexOfNextCluster = len(self.artifacts.getClusterList()) + 1
        self.listOfStrains = self.artifacts.listOfStrains
        self.list_of_new_clusters = []

    def getTheBiggestParalog(self):
        strain_in_cluster = self.artifacts.strainsPerCluster[self.cluster_index]
        strain_and_freq_dict = {}

        for x in strain_in_cluster:
            strain_and_freq_dict[x[0]] = x[1]

        # strain index of paralog
        biggest_paralog = max(strain_and_freq_dict, key=strain_and_freq_dict.get)

        return biggest_paralog

    # get the members from the paralog to open new clusters and be their representative
    def getMembersFromStrain(self, strain_index):
        # dict of all members that contain to the paralog
        members_from_strain = {}
        for member in self.cluster_members.keys():
            if self.cluster_members[member].strain_index == strain_index:
                members_from_strain[member] = self.cluster_members[member]
        return members_from_strain

    # insert the strain of the paralog, and the cluster index to split. -> and return the rest of the members in the cluster
    # to split to news.
    def getMembersToReClustering(self, strain_index):
        members_from_other_strains = {}
        for member in self.cluster_members.keys():
            if self.cluster_members[member].strain_index != strain_index:
                members_from_other_strains[member] = self.cluster_members[member]
        return members_from_other_strains

    def createNewCluster(self, paralog_dict, paralog_strain):
        self.listOfStrains.get(paralog_strain).deleteFromClusterList(self.cluster_index)  # remove cluster from the strain list
        self.listOfStrains.get(paralog_strain).addToClusterList(self.indexOfNextCluster)  # add cluster to the strain list

        # members_from_strain (of paralog)
        for member in paralog_dict.keys():
            value = {}  # dict of members
            self.cluster_members[member].represntative = True
            self.cluster_members[member].identity = 100.0

            value[0] = self.cluster_members[member]
            self.artifacts.listOfClusters.dict_clusters_setter(self.indexOfNextCluster, value)# add the new cluster
            self.list_of_new_clusters.append(self.indexOfNextCluster)
            self.indexOfNextCluster = self.indexOfNextCluster + 1
            # del self.cluster_members[member]  # remove member from the cluster list

        self.artifacts.listOfClusters.removeCluster(self.cluster_index)
        return self.list_of_new_clusters

    # members_dict- members_from_other_strains
    def selectNewClusters(self, members_dict):
        maaping_members_to_new_clusters = {}
        for m in members_dict:
            member_to_associate = members_dict[m]
            member_neighbors_in_old_cluster = []
            num_of_common_neighbors = {}

            for x in self.artifacts.getNeighborsClusters(member_to_associate).values():
                member_neighbors_in_old_cluster.extend(x) # clusters of neighbors of member_to_associate

            for cluster_index in self.list_of_new_clusters:
                members_neighbors_in_new_cluster = []
                members_in_new_cluster = self.artifacts.listOfClusters.getClusterMembers(cluster_index)
                # the neighbors of the members in new clusters
                for member in members_in_new_cluster.values():
                    for list in self.artifacts.getNeighborsClusters(member).values():
                        members_neighbors_in_new_cluster.extend(list) # clusters of neighbors of members in new cluster

                members_neighbors_in_new_cluster_set = set(members_neighbors_in_new_cluster)
                member_neighbors_in_old_cluster_set = set(member_neighbors_in_old_cluster)
                intersection = set.intersection(members_neighbors_in_new_cluster_set, member_neighbors_in_old_cluster_set)
                num_of_common_neighbors[cluster_index] = len(intersection)

            max_common_neighbors = max(num_of_common_neighbors, key=num_of_common_neighbors.get)
            if num_of_common_neighbors[max_common_neighbors] == 0:
                self.newClusterForMemberWithoutCommonNeighbors(m)
                continue
            maaping_members_to_new_clusters[member_to_associate] = max_common_neighbors

        return maaping_members_to_new_clusters

    def addingToNewCluster(self, members_dict):
        maaping_members_to_new_clusters = self.selectNewClusters(members_dict)
        for member in maaping_members_to_new_clusters:
            if member is not None:
                size_members_befor = len(self.artifacts.listOfClusters.getClusterMembers(maaping_members_to_new_clusters[member]))
                self.artifacts.listOfClusters.getClusterMembers(maaping_members_to_new_clusters[member])[size_members_befor] = member

    def newClusterForMemberWithoutCommonNeighbors(self, member):
        self.listOfStrains.get(self.cluster_members[member].getStrainInd).addToClusterList(self.indexOfNextCluster)  # add cluster to the strain list
        value = {}  # dict of members
        self.cluster_members[member].represntative = True
        self.cluster_members[member].identity = 100.0
        value[0] = self.cluster_members[member]
        self.artifacts.listOfClusters.dict_clusters_setter(self.indexOfNextCluster, value)  # add the new cluster
        self.list_of_new_clusters.append(self.indexOfNextCluster)
        self.indexOfNextCluster = self.indexOfNextCluster + 1
