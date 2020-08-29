import csv
import pandas as pd
from collections import Counter
from Artifacts import Artifact
from Strain import Strain
from ProteinFilesManager import ProteinFilesManager


class StrainsReports:
    global artifacts, dict_strains

    # starin report check!!!!!!!!!
    def __init__(self, strains, artifacts):
        self.dict_strains = strains
        self.artifacts = artifacts

    def downloadStrainSingletonsReport(self):
        singleton_strains = []
        singletons = self.artifacts.getSingleClusters()
        clustersFromClass2 = self.artifacts.listOfClass2
        merge = set(singletons + clustersFromClass2)
        countOfSingletons = len(merge)
        # print(singletons)
        # print(clustersFromClass2)
        # print(merge)

        for cluster in merge:
            members = self.artifacts.listOfClusters.getClusterMembers(cluster)
            for member in members.values():
                singleton_strains.append(member.getStrainInd)
                # if member.getStrainInd in self.dict_strains.keys():
                    # self.dict_strains.get(member.getStrainInd).setNumOfSingleton(1)
                break

        # for cluster in singletons:
        #     members = self.artifacts.listOfClusters.getClusterMembers(cluster)
        #     for member in members.values():
        #         singleton_strains.append(member.getStrainInd)
        #         if member.getStrainInd in self.dict_strains.keys():
        #             self.dict_strains.get(member.getStrainInd).setNumOfSingleton(1)
        #
        # for cluster in clustersFromClass2:
        #     members = self.artifacts.listOfClusters.getClusterMembers(cluster)
        #     for member in members.values():
        #         singleton_strains.append(member.getStrainInd)
        #         if member.getStrainInd in self.dict_strains.keys():
        #             self.dict_strains.get(member.getStrainInd).setNumOfSingleton(1)

        # print(singleton_strains)
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

    # def updateSingletonsStrainCount(self):
    #     singleton_strains = []
    #     clustersFromClass1 = self.artifacts.getSingleClusters()
    #     # countOfSingletonsClass2 = len(clustersFromClass2)
    #     # print(singletons)
    #     # print(clustersFromClass2)
    #     # print(merge)
    #
    #     for cluster in clustersFromClass1:
    #         members = self.artifacts.listOfClusters.getClusterMembers(cluster)
    #         for member in members.values():
    #             if member.getStrainInd in self.dict_strains.keys():
    #                 self.dict_strains.get(member.getStrainInd).setNumOfSingleton(1)
    #             break


    def count_of_singeltons_class2_per_strain(self):
        singleton_class2_strains = []
        clustersFromClass2 = self.artifacts.listOfClass2

        for cluster in clustersFromClass2:
            members = self.artifacts.listOfClusters.getClusterMembers(cluster)
            for member in members.values():
                singleton_class2_strains.append(member.getStrainInd)
                # if member.getStrainInd in self.dict_strains.keys():
                #     self.dict_strains.get(member.getStrainInd).setNumOfSingleton(1)
                break

        df = pd.DataFrame(singleton_class2_strains, columns=['strain index'])
        counts = df['strain index'].value_counts().to_dict()
        return counts

    def downloadStrainReport(self):
        self.artifacts.updateNumOfCoreGeneInStrains()
        self.artifacts.updateSingletonsStrainCount()
        countOfclass2 = self.count_of_singeltons_class2_per_strain()
        with open('strain reports/strain_Report.csv', mode='w') as strain_report_csv:  # TODO: change the file name
            strain_report_writer = csv.writer(strain_report_csv, delimiter=',', quotechar='"',
                                              quoting=csv.QUOTE_MINIMAL)
            # the first row in file
            strain_report_writer.writerow(['strain index', 'strain name', '# genes', '# core genes',
                                           '% core genes', '# singletons class 1', '# singletons class 2',
                                           '# outlier genes (length, class 4)', '# outlier genes (length, class 0)',
                                           '# outlier genes (members, class 0)', 'Recommend to exclude'])
            # the last 2 columns will filled later

            # print(countOfclass2)
            # print(countOfclass2.items())
            # print(countOfclass2)

            for key, val in self.dict_strains.items():
                countOfSingletonsClass2 = 0
                if (key in countOfclass2.keys()):
                    # print(countOfclass2[key])
                    countOfSingletonsClass2 = countOfclass2[key]

                strain = self.dict_strains.get(key)
                if self.artifacts.sum_of_core_clusters == 0:
                    percent_core_clusters = 0
                elif self.artifacts.sum_of_core_clusters != 0:
                    percent_core_clusters = (strain.numOfCoreGenes / self.artifacts.sum_of_core_clusters) * 100
                strain_report_writer.writerow([key, strain.name, strain.getNumOfGenes(), strain.numOfCoreGenes,
                                               percent_core_clusters,
                                               strain.getNumOfSingletons(), countOfSingletonsClass2,
                                               strain.numOfOutliers_class4_length, strain.numOfOutliers_class0_length,
                                               strain.numOfOutliers_class0_members, 'recommend to exclude_later'])
            strain_report_csv.close()




# rr = StrainsReports()
# rr.read_strains_file("C:\\Users\\Paz\\Desktop\\פרוייקט\\FinalProject\\seq_index_new")