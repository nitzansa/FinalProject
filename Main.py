from Artifacts import Artifact
from CDHIT_Parser import CDHIT_Parser
from ProteinFilesManager import ProteinFilesManager
from StrainReport import StrainsReports
from Reports import Reports
from zipfile import ZipFile
import subprocess, sys


class Main:

    def main(self):
        return
        # try:
        #     with ZipFile("/home/local/BGU-USERS/sabagnit/Data_project.zip", 'r') as zip:
        #         # printing all the contents of the zip file
        #         zip.printdir()
        #
        #         # extracting all the files
        #         print('Extracting all the files now...')
        #         zip.extractall("/home/local/BGU-USERS/sabagnit/Dataset")
        #         print('Done!')
        # except:
        #     zipFile = "/home/local/BGU-USERS/sabagnit/Data_project.zip"
        #     destinationDirectory = '/home/local/BGU-USERS/sabagnit/Dataset'
        #     print("An exception occurred extracting with Python ZipFile library.")
        #     print("Attempting to extract using 7zip")
        #     subprocess.Popen(["7z", "e", f"{zipFile}", f"-o{destinationDirectory}", "-y"])
        # protein_file_manager = ProteinFilesManager()
        # a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")

        # strains = protein_file_manager.read_strains_file("seq_index_new", "/home/local/BGU-USERS/sabagnit/Data_project")
        # cdhit_Parser = CDHIT_Parser("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster", strains)
        # clusters = cdhit_Parser.clusters

        # artifacts = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # artifacts = Artifact("23cluster", strains)

        # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()

        # strain_report = StrainsReports(strains, artifacts)
        # strain_report.downloadStrainSingletonsReport()
        # strain_report.downloadStrainReport()

        # r = Reports("23cluster", strains)

        # r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence", strains)
        # r.downloadReport()



        #a = Artifact("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
        # artifacts = Artifact("resources/23cluster")
        # artifacts = Artifact("C:\\Users\\Paz\\Desktop\\test for project\\FinalProject\\23cluster")

        # for v in strains.values():
        #     print(v.index)
        #     for c in v.getClusterList():
        #         print(c)


        # print('clusters:')
        # print(clusters)
        # print()
        # print()
        # print('strains: ')
        # print(strains)


main = Main()
main.main()