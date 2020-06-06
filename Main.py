from Artifacts import Artifact
from CDHIT_Parser import CDHIT_Parser
from ProteinFilesManager import ProteinFilesManager
from StrainReport import StrainsReports
from Reports import Reports


class Main:

    def main(self):
        protein_file_manager = ProteinFilesManager()
        # a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")

        # strains = protein_file_manager.read_strains_file("seq_index_new", "Dataset")
        strains = protein_file_manager.read_strains_file("seq_index_new", "/home/local/BGU-USERS/sabagnit/Data_project/Dataset")
        # cdhit_Parser = CDHIT_Parser("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster", strains)
        # clusters = cdhit_Parser.clusters

        # artifacts = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # artifacts = Artifact("resources/23cluster", strains)

        # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()

        # strain_report = StrainsReports(strains, artifacts)
        # strain_report.downloadStrainSingletonsReport()
        # strain_report.downloadStrainReport()

        # r = Reports("resources/23cluster", strains)

        r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()
        # i = 0
        # flag = 0
        # for cluster in r.artifacts.listOfClusters.clusters.keys():
        #     flag = 0
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
        r.downloadClusterInfo('566')



        #a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
        # artifacts = Artifact("resources/23cluster")
        # artifacts = Artifact("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster")



main = Main()
main.main()