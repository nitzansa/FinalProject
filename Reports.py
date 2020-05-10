import csv
import pandas as pd
from collections import Counter
from Artifacts import Artifact

class Reports:

    global artifacts

    def __init__(self, path):
        self.artifacts = Artifact(path)

    def downloadReport(self):
        with open('report one member.csv', mode='w') as report_one_member_csv:  # TODO: change the file name
            report_one_member_writer = csv.writer(report_one_member_csv, delimiter=',', quotechar='"',
                                                  quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            report_one_member_writer.writerow(['cluster num',
                                               'mean length',
                                               'std length',
                                               '# strains in each cluster',
                                               '# of members in each cluster',
                                               'mean of the per strain members',
                                               'min length',
                                               'max length',
                                               'min number of members per strains',
                                               'max number of members per strain',
                                               'flag'])
        report_one_member_csv.close()

        with open('report.csv', mode='w') as report_csv:  # TODO: change the file name
            report_writer = csv.writer(report_csv, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            report_writer.writerow(['cluster num',
                                    'mean length',
                                    'std length',
                                    '# strains in each cluster',
                                    '# of members in each cluster',
                                    'mean of the per strain members',
                                    'min length',
                                    'max length',
                                    'min number of members per strains',
                                    'max number of members per strain',
                                    'flag'])

            for cluster in self.artifacts.listOfClusters.clusters.keys():
                dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster)
                # if for this cluster exist only one member
                if len(dict_members) < 2:
                    self.reportToClustersWithOneMember(cluster)
                else:
                    flag = 0
                    if len(self.artifacts.strainsPerCluster[cluster]) == 1 and len(
                            self.artifacts.listOfClusters.getClusterMembers(cluster)) > 1:
                        flag = 2
                    if self.artifacts.getMaxMembersPerStrainPerCluster(cluster) == 1:
                        flag = 3
                    report_writer.writerow([cluster,
                                            self.artifacts.mean[cluster],
                                            self.artifacts.std[cluster],
                                            len(self.artifacts.strainsPerCluster[cluster]),
                                            len(self.artifacts.listOfClusters.getClusterMembers(cluster)),
                                            self.artifacts.avgMembersPerCluster[cluster],
                                            self.artifacts.minMemberLength[cluster],
                                            self.artifacts.maxMemberLength[cluster],
                                            self.artifacts.getMinMembersPerStrainPerCluster(cluster),
                                            self.artifacts.getMaxMembersPerStrainPerCluster(cluster),
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

    def calculatingLengthDistributionOfEachCluster(self): #top 3

        with open('cluster reports/frequency length.csv',
                  mode='w') as cluster_freq__csv:  # TODO: change the file name
            cluster_freq_length_writer = csv.writer(cluster_freq__csv, delimiter=',', quotechar='"',
                                                    quoting=csv.QUOTE_MINIMAL)

            # the first row in file
            cluster_freq_length_writer.writerow(['cluster', 'length', 'count of members', '% of members'])
            for cluster in self.artifacts.listOfClusters.clusters:
                length_freq = []
                dict_members = self.artifacts.listOfClusters.getClusterMembers(cluster)
                for member in dict_members.values():
                    length_freq.append(member.getLength)
                df = pd.DataFrame(length_freq, columns=['strain index'])
                counts = df['strain index'].value_counts().to_dict()
                c = Counter(counts)
                top3 = c.most_common(3)
                for i in top3:
                    cluster_freq_length_writer.writerow([cluster, i[0], i[1], (i[1] / len(dict_members)) * 100])

        cluster_freq__csv.close()


r = Reports("/home/local/BGU-USERS/sabagnit/CD_HIT_output_sqeuence")
# r.downloadReport()
# r.downloadClassifyReport()
# r.downloadStrainSingletonsReport()
r.calculatingLengthDistributionOfEachCluster()