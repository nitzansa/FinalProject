from Artifacts import Artifact
from CDHIT_Parser import CDHIT_Parser
from ProteinFilesManager import ProteinFilesManager
from StrainReport import StrainsReports
from Reports import Reports


class Main:

    def main(self):
        protein_file_manager = ProteinFilesManager()
        # a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")

        strains = protein_file_manager.read_strains_file("seq_index_new", "/home/local/BGU-USERS/sabagnit/Data_project/Dataset")
        # cdhit_Parser = CDHIT_Parser("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster", strains)
        # clusters = cdhit_Parser.clusters

        artifacts = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # artifacts = Artifact("23cluster", strains)

        # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()

        strain_report = StrainsReports(strains, artifacts)
        strain_report.downloadStrainSingletonsReport()
        strain_report.downloadStrainReport()

        # r = Reports("23cluster", strains)

        # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()



        #a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
        # artifacts = Artifact("resources/23cluster")
        # artifacts = Artifact("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster")



main = Main()
main.main()