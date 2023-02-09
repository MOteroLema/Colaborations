import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase



def parser_atomtypes(line, desired_atoms, types_dic):
        fields = line.split()
        if fields[0] in desired_atoms:
            types_dic[fields[0]]=(fields[6], fields[7])

def fun_pass(line, desired_atoms, types_dic):
    pass
        
def vdw_reader(path, desired_atoms, types_dic = None):
    

    with open(path) as file:
        lines = file.readlines()
    if types_dic==None:
        types_dic = {}

    for line in lines:
        if "[" in line and "]" in line:
            section = line.split('[')[1].split(']')[0].strip()
            if section == "atomtypes":
                parser=parser_atomtypes
            else:
                parser=fun_pass
        elif line.split()[0]=="#include":
            path_include = line.split()[1][1:-2]
            print(path_include)
            types_dic = vdw_reader(path_include, desired_atoms, types_dic)
        elif line[0]==";" or line[0]=="#":
            pass
        else:
            parser(line, desired_atoms, types_dic)
    return types_dic
    
class ExcessPotential(AnalysisBase):

    def __init__(self, tpr_file, traj_file, path_vdw, n_insertions, inserted_radius, sigma, epsilon, mixing_rule = "geometric"):

        self.u = mda.Universe(tpr_file, traj_file)
        self.n_insertions = n_insertions
        self.inserted_radius = inserted_radius
        self.sigma = sigma
        self.epsilon = epsilon
        self.atomtypes = np.unique(self.u.atoms.types)
        self.atom_parameters = vdw_reader(path_vdw, self.atomtypes)
        
        @classmethod
        def LJ(r, eps, sigma):

            return 4 * eps * ( (sigma/r)**12 - (sigma/r)**6 )
        


        