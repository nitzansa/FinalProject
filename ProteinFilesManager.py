import zipfile

import Bio
from Bio.Seq import Seq
from Strain import Strain
import pandas as pd


class ProteinFilesManager:

    global strains

    def __init__(self):
        self.strains = {}

    def read_strains_file(self, strains_path, protein_path):
        try:
            file = open(strains_path + ".txt", "r")
        except IOError:
            print("could not open the file")
        with file:
            line = file.readline()
            strain_line = line.split(", ")
            for s in strain_line:
                split_strain = s.split(": ")
                strain_index = int(split_strain[0].replace("{", "").replace("}", "").replace("\"", ""))
                strain_name = split_strain[1].replace("{", "").replace("}", "").replace("\"", "")
                proteins_df = self.read_proteins_file(protein_path, strain_name)
                if proteins_df is not None:
                    proteins_dict = proteins_df.to_dict()
                    strain = Strain(strain_index, strain_name, proteins_dict)
                    strain.numOfGenes = len(proteins_df.index)
                    self.strains[strain.index] = strain

        # print(self.strains['4580'].proteins[0])


        return self.strains

    def read_proteins_file(self, path, strain_name):
        curr_path = path + "/" + strain_name + "/data.csv"
        try:
            df = pd.read_csv(curr_path, usecols=['locus_tag', 'name_y'])
        except IOError:
            # print("could not open the file")
            return

        return df

    def getStrain(self, strain_index):
        return self.strains[strain_index]






# a = ProteinFilesManager()
# a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")
# a.read_strains_file("seq_index_example", "Dataset")
