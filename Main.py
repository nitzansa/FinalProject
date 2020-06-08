from Artifacts import Artifact
from CDHIT_Parser import CDHIT_Parser
from Member import Member
from ProteinFilesManager import ProteinFilesManager
from StrainReport import StrainsReports
from Reports import Reports
from clusterCorrection import clusterCorrection


class Main:

    def main(self):
        protein_file_manager = ProteinFilesManager()
        # a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")

        # strains = protein_file_manager.read_strains_file("seq_index_new", "Dataset")
        strains = protein_file_manager.read_strains_file("seq_index_new", "/home/local/BGU-USERS/sabagnit/Data_project/Dataset")
        # cdhit_Parser = CDHIT_Parser("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster", strains)
        # clusters = cdhit_Parser.clusters

        # artifacts = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # artifacts = Artifact("cluster6", strains)

        # clusterC = clusterCorrection(artifacts)
        # clusterC.getTheBiggestParalog(6)
        # c = clusterC.getMembersFromStrain(145,6)
        # print(c)
        # del c[0]
        # print(c)
        # for a in c.items():
        #     print(a)

        # artifacts = Artifact("resources/cluster6", strains)
        # member = Member(145, 2937, 'Q003_RS36535', 415, 100, False)
        # testt = artifacts.getNeighboursClusters(member)
        # print()
        # artifacts.calcAverageMemberPerCluster()

        r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()

        # strain_report = StrainsReports(strains, artifacts)
        # strain_report.downloadStrainSingletonsReport()
        # strain_report.downloadStrainReport()

        # r = Reports("cluster6", strains)
        r.downloadLengthDistributionForCluster(22724)
        r.downloadLengthDistributionForCluster(15746)
        r.downloadLengthDistributionForCluster(6573)

        # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()
        # i = 0
        # flag = 0
        # for cluster in r.artifacts.listOfClusters.clusters.keys():
        #     i = i + 1
        #     if i <= 1000:
        #         dict_members = r.artifacts.listOfClusters.getClusterMembers(cluster)
        #         # if for this cluster exist only one member
        #         if len(dict_members) < 2:
        #             True
        #         else:
        #             if r.artifacts.getMaxMembersPerStrainPerCluster(cluster) == 1:
        #                 flag = 3
        #             if flag == 3 and len(r.artifacts.strainsPerCluster[cluster]) >= 50 and 30 <= r.most_common_length_dict[cluster]['%_1'] < 100:
        #                 r.downloadClusterInfo(str(cluster))
        # r.downloadClusterInfo('6')



        #a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
        # artifacts = Artifact("resources/23cluster")
        # artifacts = Artifact("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster")



main = Main()
main.main()