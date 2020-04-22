from Member import Member


class CDHIT_Parser:

    global dict_clusters

    def __init__(self, path):
        self.dict_clusters = {}
        self.readFile(path)

    # a getter function
    @property
    def clusters(self):
        return self.dict_clusters

    def readFile(self, path):
        try:
            file = open(path + ".txt", "r")
        except IOError:
            print("could not open the file")
        with file:
            line = file.readline().split(" ")
            while line.__len__() > 1:
                if line[0][0] == '>':
                    key = int(line[1])
                    value = {} #dict of members
                    line = file.readline().split(" ")
                    while line.__len__() > 2:
                        strain_index = int(line[1].split("|")[0].replace(">", ""))
                        protein_index = int(line[1].split("|")[1].replace("...", ""))
                        length = int(line[0].split("\t")[1].replace("aa,", ""))
                        if line[2] == '*\n' or line[2] == '*':
                            represntative = True
                            identity = 100.0
                        else:
                            identity = float(line[3].split("/")[1].replace("%", "").replace("\n", ""))
                            represntative = False
                        new_member = Member(strain_index, protein_index, length, identity, represntative)
                        value[int(line[0].split("\t")[0])] = new_member
                        line = file.readline().split(" ")
                    self.dict_clusters[key] = value
        file.close()

    def getClusterMembers(self, cluster_index):
        return self.dict_clusters[cluster_index]


a = CDHIT_Parser("clusters_output")
a.clusters
membersOfOneCluster = a.getClusterMembers(14)
#print("hj")

