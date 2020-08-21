from Artifacts import Artifact
from CDHIT_Parser import CDHIT_Parser
from Member import Member
from ProteinFilesManager import ProteinFilesManager
from StrainReport import StrainsReports
from Reports import Reports
from clusterCorrection import clusterCorrection
from OutlierDetection import OutlierDetection


class Main:
    # main check!!!!!!!!!
    def main(self):
        # protein_file_manager = ProteinFilesManager()
        #         # a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")
        #
        #         # strains = protein_file_manager.read_strains_file("seq_index_new", "Dataset")
        #         # strains = protein_file_manager.read_strains_file("resources/seq_index_new", "/home/local/BGU-USERS/sabagnit/Data_project/Dataset")
        #         # cdhit_Parser = CDHIT_Parser("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster", strains)
        #         # clusters = cdhit_Parser.clusters
        #
        #         # artifacts = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        #         # artifacts = Artifact("cluster6", strains)
        #
        #         # clusterC = clusterCorrection(artifacts)
        #         # clusterC.getTheBiggestParalog(6)
        #         # c = clusterC.getMembersFromStrain(145,6)
        #         # print(c)
        #         # del c[0]
        #         # print(c)
        #         # for a in c.items():
        #         #     print(a)
        #
        #         # artifacts = Artifact("resources/cluster6", strains)
        #         # member = Member(145, 2937, 'Q003_RS36535', 415, 100, False)
        #         # testt = artifacts.getNeighborsClusters(member)
        #         # print(testt)
        #         # artifacts.calcAverageMemberPerCluster()
        #
        #         # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        #         # r.downloadReport()
        #
        #         # strain_report = StrainsReports(strains, artifacts)
        #         # strain_report.downloadStrainSingletonsReport()
        #         # strain_report.downloadStrainReport()
        #
        #         # r = Reports("cluster6", strains)
        #         # r.downloadLengthDistributionForCluster(22724)
        #         # r.downloadLengthDistributionForCluster(15746)
        #         # r.downloadLengthDistributionForCluster(6573)
        #
        #         # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        #         # r.downloadReport()
        #         # i = 0
        #         # flag = 0
        #         # for cluster in r.artifacts.listOfClusters.clusters.keys():
        #         #     i = i + 1
        #         #     if i <= 1000:
        #         #         dict_members = r.artifacts.listOfClusters.getClusterMembers(cluster)
        #         #         # if for this cluster exist only one member
        #         #         if len(dict_members) < 2:
        #         #             True
        #         #         else:
        #         #             if r.artifacts.getMaxMembersPerStrainPerCluster(cluster) == 1:
        #         #                 flag = 3
        #         #             if flag == 3 and len(r.artifacts.strainsPerCluster[cluster]) >= 50 and 30 <= r.most_common_length_dict[cluster]['%_1'] < 100:
        #         #                 r.downloadClusterInfo(str(cluster))
        #         # r.downloadClusterInfo('6')



        #a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
        # artifacts = Artifact("resources/23cluster")
        # artifacts = Artifact("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster")
        # first = [1,2,3]
        # second = [2,3,4,5,6]
        # join = set(first + second)
        # print(first)
        # print(second)
        # print(join)

        #######################################################
        protein_file_manager = ProteinFilesManager()
        strains = protein_file_manager.read_strains_file("resources/seq_index_new", "/home/local/BGU-USERS/sabagnit/Data_project/Dataset")
        # strains = protein_file_manager.read_strains_file("resources/seq_index_new", "Dataset")

        artifacts = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # artifacts = Artifact("resources/cluster6", strains)
        # artifacts = Artifact("resources/clusterNew", strains)
        r = Reports(artifacts)
        # r.downloadClusterInfo('6573')
        # clusterC = clusterCorrection(artifacts, 6573)

        r.downloadReport()
        # clusterC = clusterCorrection(artifacts, 0)
        # paralog_index = clusterC.getTheBiggestParalog()
        # print('get biggest paralog')
        # print(paralog_index)
        # other_strains = clusterC.getMembersToReClustering(paralog_index)
        # list = clusterC.createNewCluster(clusterC.getMembersFromStrain(paralog_index), paralog_index)
        # clusterC.addingToNewCluster(other_strains)
        # print('before download all new clusters info')
        # for new_cluster in list:
        #     # print('clusterInfo')
        #     # print(new_cluster)
        #     r.downloadClusterInfo(str(new_cluster))

        outlierDetection = OutlierDetection(artifacts)
        outlierDetection.detectOutlier()

        strain_report = StrainsReports(strains, artifacts)
        strain_report.downloadStrainSingletonsReport()
        strain_report.downloadStrainReport()

        # print(len(artifacts.listOfClass0))
        # print(len(artifacts.singletons))
        # print(len(artifacts.listOfClass2))
        # print(len(artifacts.listOfClass3))
        # print(len(artifacts.listOfClass4))
        # print(len(artifacts.listOfClass5))
        # print(artifacts.listOfClass4)
        # print(artifacts.listOfClass0)
        # print(artifacts.listOfClass2)
        # print(artifacts.listOfClass5)
        # print(artifacts.singletons)


main = Main()
main.main()