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
                proteins = self.read_proteins_file(protein_path, strain_name)
                strain = Strain(strain_index, strain_name, proteins)
                if proteins is not None:
                    strain.numOfGenes = proteins.size
                self.strains[strain.index] = strain

        # print(self.strains['4580'].proteins[0])


        return self.strains

    def read_proteins_file(self, path, strain_name):
        with zipfile.ZipFile(path + ".zip") as z:
            name_of_folder = z.namelist()[0]
            with z.open(name_of_folder + strain_name + "/data.csv") as f:
                try:
                    df = pd.read_csv(f, usecols=['locus_tag'])
                    return df
                except IOError:
                    print("could not open the file")

        return "could not open the file"

    def getStrain(self, strain_index):
        return self.strains[strain_index]






a = ProteinFilesManager()
a.read_proteins_file("Dataset", "GCF_901472595.1_36340_C01")
# a.read_strains_file("seq_index_example", "Dataset")
