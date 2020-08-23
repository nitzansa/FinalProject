import csv
import statistics

import pandas as pd
from collections import Counter
from Artifacts import Artifact

class Reports:

    global artifacts, most_common_length_dict
    # report check!!!!!!
    def __init__(self, artifacts):

        self.artifacts = artifacts
        # self.most_common_length_dict = self.calculatingLengthDistributionOfEachCluster()

    def downloadReport(self):
        self.most_common_length_dict = self.artifacts.most_common_length_dict
        # with open('report one member.csv', mode='w') as report_one_member_csv:  # TODO: change the file name
        #     report_one_member_writer = csv.writer(report_one_member_csv, delimiter=',', quotechar='"',
        #                                           quoting=csv.QUOTE_MINIMAL)
        #
        #     # the first row in file
        #     report_one_member_writer.writerow(['cluster num',
        #                                        'mean length',
        #                                        'std length',
        #                                        '# strains in each cluster',
        #                                        '# of members in each cluster',
        #                                        'mean of the per strain members',
        #                                        'min length',
        #                                        'max length',
        #                                        'min number of members per strains',
        #                                        'max number of members per strain',
        #                                        'flag'])
        # report_one_member_csv.close()

        with open('cluster reports/clusters report.csv', mode='w') as report_csv:  # TODO: change the file name
            report_writer = csv.writer(report_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            report_writer.writerow(['cluster num',
                                    'mean length',
                                    'std length',
                                    '# strains in each cluster',
                                    '# of members in each cluster',
                                    'mean of the per strain members',
                                    'std of the per strain members',
                                    'min length',
                                    'max length',
                                    'min number of members per strains',
                                    'max number of members per strain',
                                    'core cluster',
                                    'average identity score',
                                    'std identity score',
                                    '% of members with score 70-80%',
                                    '% of members with score 80-90%',
                                    '% of members with score 90-100%',
                                    'most common length_1',
                                    '% of members_1',
                                    'most common length_2',
                                    '% of members_2',
                                    'most common length_3',
                                    '% of members_3',
                                    'flag'])

            for cluster in self.artifacts.listOfClusters.clusters.keys():
                dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster)
                # if for this cluster exist only one member
                if len(dict_members) < 2:
                    True
                    # self.reportToClustersWithOneMember(cluster)
                # elif len(self.artifacts.strainsPerCluster[cluster]) == 1 and len(
                #             self.artifacts.listOfClusters.getClusterMembers(cluster)) > 1:
                # #     flag ==2
                #     True
                else:
                    flag = 0
                    # write clusters from class 2 in the report, to check if it can be deleted?
                    if len(self.artifacts.strainsPerCluster[cluster]) == 1 and len(
                            self.artifacts.listOfClusters.getClusterMembers(cluster)) > 1:
                        flag = 2
                    if self.artifacts.getMaxMembersPerStrainPerCluster(cluster) == 1:
                        flag = 3
                    if flag == 3 and 30 <= self.most_common_length_dict[cluster]['%_1'] < 100:
                        flag = 4
                    if flag == 3 and self.most_common_length_dict[cluster]['%_1'] < 30:
                        flag = 5
                    is_core = False
                    if self.artifacts.isCoreCluster(cluster):
                        is_core = True
                        self.artifacts.increase_sum_of_core_clusters()
                    report_writer.writerow([cluster,
                                            self.artifacts.mean[cluster],
                                            self.artifacts.std[cluster],
                                            len(self.artifacts.strainsPerCluster[cluster]),
                                            len(self.artifacts.listOfClusters.getClusterMembers(cluster)),
                                            self.artifacts.avgMembersPerCluster[cluster],
                                            self.artifacts.getSTDMembersPerStrainPerCluster(cluster),
                                            self.artifacts.minMemberLength[cluster],
                                            self.artifacts.maxMemberLength[cluster],
                                            self.artifacts.getMinMembersPerStrainPerCluster(cluster),
                                            self.artifacts.getMaxMembersPerStrainPerCluster(cluster),
                                            str(is_core),
                                            str(self.artifacts.avgIdentity(cluster)),
                                            str(self.artifacts.stdIdentity(cluster)),
                                            str(self.artifacts.PercentOfMembersWithC_Score(cluster)),
                                            str(self.artifacts.PercentOfMembersWithB_Score(cluster)),
                                            str(self.artifacts.PercentOfMembersWithA_Score(cluster)),
                                            self.most_common_length_dict[cluster]['length_1'],
                                            self.most_common_length_dict[cluster]['%_1'],
                                            self.most_common_length_dict[cluster]['length_2'],
                                            self.most_common_length_dict[cluster]['%_2'],
                                            self.most_common_length_dict[cluster]['length_3'],
                                            self.most_common_length_dict[cluster]['%_3'],
                                            flag])
        report_csv.close()


    def reportToClustersWithOneMember(self, cluster):
        with open('report one member.csv', mode='a') as report_csv:  # TODO: change the file name
            report_writer = csv.writer(report_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster)
            member_length = 0
            for member in dict_members.values():
                member_length = member.getLength
            report_writer.writerow([cluster,
                                    member_length,
                                    '0',
                                    '1',
                                    '1',
                                    '1',
                                    member_length,
                                    member_length,
                                    '1',
                                    '1',
                                    '1'])
        report_csv.close()


    def classifyCluster(self):
        curr_path = "resources/report.csv"
        try:
            df = pd.read_csv(curr_path, usecols=['flag'])
            counts = df['flag'].value_counts().to_dict()
            counts[1] = len(self.artifacts.getSingleClusters())
            return counts
        except IOError:
            print("could not open the file")
            return


    def downloadClassifyReport(self):
        with open('cluster reports/classify.csv', mode='w') as classify__csv:  # TODO: change the file name
            classify_writer = csv.writer(classify__csv, delimiter=',', quotechar='"',
                                         quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            classify_writer.writerow(['class type', 'count'])
            classes_dict = self.classifyCluster()
            for key, val in classes_dict.items():
                classify_writer.writerow([key, val])
            classify__csv.close()


    def downloadStrainSingletonsReport(self):
        singleton_strains = []
        singletons = self.artifacts.getSingleClusters()
        countOfSingletons = len(singletons)
        for cluster in singletons:
            members = self.artifacts.listOfClusters.getClusterMembers(cluster)
            for member in members.values():
                singleton_strains.append(member.getStrainInd)

        df = pd.DataFrame(singleton_strains, columns=['strain index'])
        counts = df['strain index'].value_counts().to_dict()

        with open('resources/singleton strains', mode='w') as singletons_starin_csv:  # TODO: change the file name
            singleton_strains_writer = csv.writer(singletons_starin_csv, delimiter=',', quotechar='"',
                                                  quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            singleton_strains_writer.writerow(['strain ind', 'count of singletons', '% of singletons'])
            for key, val in counts.items():
                singleton_strains_writer.writerow([key, val, (val / countOfSingletons) * 100])
            singletons_starin_csv.close()

    # def calculatingLengthDistributionOfEachCluster(self): #top 3
    #     most_common_length_dict = {}
    #     for cluster in self.artifacts.listOfClusters.clusters:
    #         length_freq = []
    #         dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster)
    #         if len(dict_members) > 1:
    #             for member in dict_members.values():
    #                 length_freq.append(member.getLength)
    #             df = pd.DataFrame(length_freq, columns=['strain index'])
    #             counts = df['strain index'].value_counts().to_dict()
    #             c = Counter(counts)
    #             top3 = c.most_common(3)
    #             if len(top3) == 1:
    #                 most_common_length_dict[cluster] = {'length_1': top3[0][0],
    #                                                     '%_1': (top3[0][1] / len(dict_members)) * 100,
    #                                                     'length_2': 0, '%_2': 0,
    #                                                     'length_3': 0, '%_3': 0}
    #             elif len(top3) == 2:
    #                 most_common_length_dict[cluster] = {'length_1': top3[0][0],
    #                                                     '%_1': (top3[0][1] / len(dict_members)) * 100,
    #                                                     'length_2': top3[1][0],
    #                                                     '%_2': (top3[1][1] / len(dict_members)) * 100,
    #                                                     'length_3': 0, '%_3': 0}
    #             elif len(top3) == 3:
    #                 most_common_length_dict[cluster] = {'length_1': top3[0][0],
    #                                                     '%_1': (top3[0][1] / len(dict_members)) * 100,
    #                                                     'length_2': top3[1][0],
    #                                                     '%_2': (top3[1][1] / len(dict_members)) * 100,
    #                                                     'length_3': top3[2][0],
    #                                                     '%_3': (top3[2][1] / len(dict_members)) * 100}
    #     return most_common_length_dict

    def downloadClusterInfo(self, cluster_index):
        with open('cluster reports/cluster_' + cluster_index + '.csv', mode='w') as cluster_info__csv:  # TODO: change the file name
            cluster_info_writer = csv.writer(cluster_info__csv, delimiter=',', quotechar='"',
                                                  quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            cluster_info_writer.writerow(['member number',
                                          'identity',
                                          'length',
                                          'protein index',
                                          'representative',
                                          'strain index',
                                          'locus_tag',
                                          'name_y'])
            dict_members = self.artifacts.listOfClusters.getClusterMembers(int(cluster_index))
            num_of_member = 0
            for member in dict_members.values():
                cluster_info_writer.writerow([num_of_member,
                                              member.getIdentity,
                                              member.getLength,
                                              member.getProteinInd,
                                              member.getRepresntative,
                                              member.getStrainInd,
                                              self.artifacts.listOfStrains.get(member.getStrainInd).getProteins()['locus_tag'][member.getProteinInd],
                                              self.artifacts.listOfStrains.get(member.getStrainInd).getProteins()['name_y'][member.getProteinInd]])
                num_of_member = num_of_member + 1
            cluster_info__csv.close()

    def downloadLengthDistributionForCluster(self,cluster_index):
        length_freq = []
        dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster_index)

        if len(dict_members) > 1:
            with open('resources/length_distribution_cluster_' + str(cluster_index) + '.csv',
                      mode='w') as cluster_info__csv:  # TODO: change the file name
                cluster_length_info_writer = csv.writer(cluster_info__csv, delimiter=',', quotechar='"',
                                                        quoting=csv.QUOTE_MINIMAL)
                # the first row in file
                cluster_length_info_writer.writerow(['length', '# members', '% members'])
                for member in dict_members.values():
                    length_freq.append(member.getLength)
                df = pd.DataFrame(length_freq, columns=['strain index'])
                counts = df['strain index'].value_counts().to_dict()
                c = Counter(counts)
                for i in c.keys():
                    cluster_length_info_writer.writerow([i, c[i], c[i] / len(dict_members) * 100])

                cluster_info__csv.close()



# r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
# r.downloadReport()

# r = Reports('resources/23cluster')

# r.downloadClassifyReport()
# r.downloadStrainSingletonsReport()
# r.calculatingLengthDistributionOfEachCluster()