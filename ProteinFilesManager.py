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
                split_strein = s.split(": ")
                strain_name = split_strein[0].replace("{", "").replace("}", "").replace("\"", "")
                strain_ind = int(split_strein[1].replace("{", "").replace("}", "").replace("\"", ""))
                proteins = self.read_proteins_file(protein_path, strain_name)
                strain = Strain(strain_ind, strain_name, proteins)
                self.strains[strain.index] = strain

        # print(self.strains['4580'].proteins[0])
        return self.strains

    def read_proteins_file(self, path, strain_name):
        curr_path = path + "/" + strain_name + "/data.csv"
        try:
            df = pd.read_csv(curr_path, usecols=['locus_tag'])
        except IOError:
            print("could not open the file")
            return

        return df








a = ProteinFilesManager()
# a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")
a.read_strains_file("seq_index_example", "Dataset")
