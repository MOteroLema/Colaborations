import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
import pathlib

    
class ExcessPotential(AnalysisBase):

    def __init__(self, tpr_file, traj_file, top_file, n_insertions, inserted_radius, sigma, epsilon, 
                 T = 293.15, p = 1, mixing_rule_sigma = "geometric", mixing_rule_epsilon = "geometric", vdw_cutoff = 1.1):

        self.u = mda.Universe(tpr_file, traj_file)
        self.n_insertions = n_insertions
        self.inserted_radius = inserted_radius
        self.sigma = sigma
        self.epsilon = epsilon
        self.atomtypes = np.unique(self.u.atoms.types)
        self.vdw_cutoff = vdw_cutoff
        self.atom_parameters = ExcessPotential.vdw_reader(top_file, self.atomtypes)
        self.KT = T * 1.380649 * 1e-23 ## SI units
        self.BETA = 1 / self.KT ## SI units
        self.NA = 6.0221408e+23
        
        functions = {"geometric": ExcessPotential.geometric_mix, "arithmetic": ExcessPotential.arithmetic_mix}

        self.mix_sigma = functions[mixing_rule_sigma]
        self.mix_epsilon = functions[mixing_rule_epsilon]


    @staticmethod
    def LJ(r, eps, sigma):
        return 4 * eps * ( (sigma/r)**12 - (sigma/r)**6 )

    @staticmethod
    def geometric_mix(a, b):
        return np.sqrt(a*b)

    @staticmethod
    def arithmetic_mix(a, b):
        return (a+b)/2

    @staticmethod
    def _parser_atomtypes(line, desired_atoms, types_dic):
        fields = line.split()
        if fields[0] in desired_atoms:
            types_dic[fields[0]]=(fields[6], fields[7])

    @staticmethod
    def vdw_reader(path, desired_atoms, types_dic = None):
        if isinstance(path, str):
            path = pathlib.Path(path)
        folder = path.parent
        with open(path) as file:
            lines = file.readlines()
        if types_dic==None:
            types_dic = {}
        for line in lines:
            if "[" in line and "]" in line:
                section = line.split('[')[1].split(']')[0].strip()
                if section == "atomtypes":
                    parser = ExcessPotential._parser_atomtypes
                else:
                    parser = lambda *args, **kwargs: None

            elif len(line.split())==0 or line[0]==";":
                pass
            elif line.split()[0]=="#include":
                included_file = line.split('"')[1]
                path_include = folder / included_file
                types_dic = ExcessPotential.vdw_reader(path_include, desired_atoms, types_dic)
            elif line[0]=="#":
                pass
            else:
                parser(line, desired_atoms, types_dic)
        return types_dic

        
    def _prepare(self):
        
        self.potential_contributions = np.zeros(self.n_frames)
        self.volumes = np.zeros(self.n_frames)

    def _single_frame(self):

        sides = self.u.dimensions[:3]
        volume = np.prod(sides)
        self.volumes[self._frame_index] = volume

        
        for _ in range(self.n_insertions):

            insertion_coordinates = np.random.random(3) * sides

            indices, distances = mda.lib.distances.capped_distance(insertion_coordinates,
                                                                   self.u.atoms.positions,
                                                                   max_cutoff = self.vdw_cutoff,
                                                                   box = self._ts.dimensions)           
            energy = 0
            for atom_indexes, distance in zip(indices, distances):

                atom_index = atom_indexes[1]
                atomtype = self.u.atoms[atom_index].type
                sig, eps = self.atom_parameters[atomtype]

                mixed_sig = self.mix_sigma(sig, self.sigma) * 10 ## In angstroms
                mixed_eps = self.mix_epsilon(eps, self.epsilon) * 1000 / self.NA  ## SI units (J)

                energy += self.LJ(distance, mixed_eps, mixed_sig)
            
            exponential = np.exp(- energy * self.BETA)
            self.potential_contributions[self._frame_index] += volume * exponential
                                                                
    def _conclude(self):

        self.potential_contributions /= self.n_insertions
        
        


        