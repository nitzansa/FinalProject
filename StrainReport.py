import csv
import pandas as pd
from collections import Counter
from Artifacts import Artifact
from Strain import Strain
from ProteinFilesManager import ProteinFilesManager


class StrainsReports:
    global artifacts, dict_strains


    def __init__(self, strains, artifacts):
        self.dict_strains = strains
        self.artifacts = artifacts

    def downloadStrainSingletonsReport(self):
        singleton_strains = []
        singletons = self.artifacts.getSingleClusters()
        countOfSingletons = len(singletons)
        for cluster in singletons:
            members = self.artifacts.listOfClusters.getClusterMembers(cluster)
            for member in members.values():
                singleton_strains.append(member.getStrainInd)
                if member.getStrainInd in self.dict_strains.keys():
                    self.dict_strains.get(member.getStrainInd).setNumOfSingleton(1)


        df = pd.DataFrame(singleton_strains, columns=['strain index'])
        counts = df['strain index'].value_counts().to_dict()

        with open('strain reports/singleton strains.csv', mode='w') as singletons_starin_csv:  # TODO: change the file name
            singleton_strains_writer = csv.writer(singletons_starin_csv, delimiter=',', quotechar='"',
                                                  quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            singleton_strains_writer.writerow(['strain ind', 'count of singletons', '% of singletons'])
            for key, val in counts.items():
                singleton_strains_writer.writerow([key, val, (val / countOfSingletons) * 100])
            singletons_starin_csv.close()




    # def read_strains_file(self, strains_file_path):
    #
    #     try:
    #         file = open(strains_file_path + ".txt", "r")
    #     except IOError:
    #         print("could not open the file")
    #     with file:
    #         line = file.readline()
    #         strain_line = line.split(", ")
    #         for s in strain_line:
    #             split_strain = s.split(": ")
    #             strain_index = int(split_strain[0].replace("{", "").replace("}", "").replace("\"", ""))
    #             strain_name = split_strain[1].replace("{", "").replace("}", "").replace("\"", "")

    def updateNumOfCoreGeneOfStrain(self):
        print()

    def downloadStrainReport(self):
        self.artifacts.updateNumOfCoreGeneInStrains()
        with open('strain reports/strain_Report.csv', mode='w') as strain_report_csv:  # TODO: change the file name
            strain_report_writer = csv.writer(strain_report_csv, delimiter=',', quotechar='"',
                                              quoting=csv.QUOTE_MINIMAL)
            # the first row in file
            strain_report_writer.writerow(['strain index', 'strain name', '# genes', '# core genes', '# singletons',
                                           '# outlier genes', 'Recommend to exclude'])
            # the last 2 columns will filled later


            for key, val in self.dict_strains.items():
                strain = self.dict_strains.get(key)
                strain_report_writer.writerow([key, strain.name, strain.getNumOfGenes(), strain.numOfCoreGenes, strain.getNumOfSingletons(),
                                               'later2', 'later3'])
            strain_report_csv.close()




# rr = StrainsReports()
# rr.read_strains_file("C:\\Users\\Paz\\Desktop\\פרוייקט\\FinalProject\\seq_index_new")