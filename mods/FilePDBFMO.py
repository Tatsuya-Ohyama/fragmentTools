#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from mods.molecule_object import Molecule
from mods.FragmentData import FragmentData



# =============== constant =============== #
RE_ATOMNAME_DECO = re.compile(r"^[\d'\s\*]")
RE_TERMH1 = re.compile(r"H[53]T")
RE_TERMH2 = re.compile(r"HO[53]")
RE_PYRIMIDINE = re.compile(r'^\s*[RD][UTC][53]?\s*$')
RE_PURINE = re.compile(r'^\s*[RD][GA][53]?\s*$')
RE_5TERM = re.compile(r'^[RD][AGCTU]5$')
RE_3TERM = re.compile(r'^[RD][AGCTU]3$')

RESIDUE_TYPES = {
	"AminoAcid": ["ACE", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIP", "HID", "HIE", "ILE", "LEU", "LYS", "MET", "NME", "PHE", "PRO", "SER", "SYM", "THR", "TRP", "TYR", "VAL"],
	"NucleicAcid": ["DA5", "DT5", "DG5", "DC5", "DA3", "DT3", "DG3", "DC3", "DA", "DT", "DG", "DC", "RA5", "RU5", "RG5", "RC5", "RA3", "RU3", "RG3", "RC3", "RA", "RU", "RG", "RC", "DNA_base"],
	"Water": ["SOL", "WAT", "HOH"],
	"Ion": ["Na", "Mg", "K", "Ca", "Cl", "Zn", "CAL", "ZIN"]
}


# =============== functions =============== #
def get_residue_type(residue_name):
	"""
	Function to return residue type

	Args:
		residue_name (str): residue name

	Returns:
		str: `AminoAcid`, `NucleicAcid`, `Water`, `Ion` or `Other`
	"""
	for key, value in residue_types.items():
		if residue_name in value:
			return key
	return "Other"



# =============== class =============== #
class FilePDBFMO:
	def __init__(self, pdb_file):
		self._obj_molecule = None
		self._fragments = []
		self._connections = []
		self._n_atom = None
		self._charge = None
		self._ligands = {}

		self.read(pdb_file)


	@property
	def n_atom(self):
		return self._n_atom

	@property
	def n_fragment(self):
		if len(self._fragments) == 0:
			self.create_fragments()
		return len(self._fragments)

	@property
	def charge(self):
		if len(self._fragments) == 0:
			self.create_fragments()
		return self._charge

	@property
	def fragments(self):
		if len(self._fragments) == 0:
			self.create_fragments()
		return self._fragments

	@property
	def connections(self):
		if len(self._fragments) == 0:
			self.create_fragments()
		return self._connections


	def read(self, input_file):
		"""
		Method to read .pdb file and create initial fragment objects

		Args:
			input_file (str): .pdb file

		Returns:
			self
		"""
		self._obj_molecule = Molecule(input_file)
		return self


	def create_fragments(self, option_nuc, option_peptide, charge_map={}):
		"""
		Method to create fragment objects

		Args:
			option_nuc (str): `+base`, `/base` or `/sugar`
			option_peptide (str): `+amino`, `/amino`, `+peptide` or `/peptide`

		Returns:
			self
		"""
		self._fragments = []
		self._connections = []

		atoms_remain = []
		fragment_idx = 0
		for res_idx, obj_residue in enumerate(self._obj_molecule.residues, 1):
			residue_type = get_residue_type(obj_residue.name)
			if residue_type == "AminoAcid":
				# Amino acid
				# charge
				if obj_residue.name in ["LYS", "ARG", "HIP", "SYM", "ACE"]:
					obj_fragment.set_charge(1)
				elif obj_residue.name in ["ASP", "GLU", "NME"]:
					obj_fragment.set_charge(-1)
				elif obj_residue.name in ["ACE", "ALA", "ASN", "CYS", "GLN", "GLY", "HIS", "HID", "HIE", "ILE", "LEU", "MET", "NME", "PHE", "PRO", "SER", "SYM", "THR", "TRP", "TYR", "VAL"]:
					obj_fragment.set_charge(0)

				# fragmentation
				new_atoms_remain = [obj_atom.number for obj_atom in obj_residue.atoms if obj_atom.name in ["C", "O"]]
				atoms_join = [obj_atom.number for obj_atom in obj_residue.atoms if obj_atom.name not in ["C", "O"]]
				obj_fragment.set_atoms(atoms_remain + atoms_join)
				atoms_remain = new_atoms_remain

				# connection
				atom_idx_c = [obj_atom.number for obj_atom in obj_residue.atoms if obj_atom.name == "C"][0][1]
				atom_idx_ca = [obj_atom.number for obj_atom in obj_residue.atoms if obj_atom.name == "CA"][0][1]
				self._connections.append([atom_idx_ca, atom_idx_c])

			elif residue_type == "NucleicAcid":
				new_atoms_remain = [obj_atom.number for obj_atom in obj_residue.atoms if obj_atom.name in "P" or obj_atom.name in ["C5'", "H5'", "H5''", "H5'1", "H5'2"]]
				if option_nuc == "/base":
					if len(self._fragments) >= 2:
						# move phosphate atoms to the previous backbone fragment
						self._fragments[-2].set_atoms(list(sorted(self._fragments[-2].atoms + new_atoms_remain)))

					atoms_bone = [obj_atom.number for obj_atom in obj_residue.atoms if "'" in obj_atom.name or "P" in obj_atom.name or RE_TERMH1.match(obj_atom.name) or RE_TERMH2.match(obj_atom.name)]
					fragment_idx += 1
					obj_fragment_bone = FragmentData()
					obj_fragment_bone.set_fragment_index(fragment_idx)
					if obj_residue.name.endswith("5"):
						obj_fragment_bone.set_charge(1)
					else:
						obj_fragment_bone.set_charge(0)
					self._connections.append([])

					obj_fragment_bone.set_atoms(atoms_bone)
					self._fragments.append(obj_fragment_bone)

					atoms_base = [obj_atom.number for obj_atom in obj_residue.atoms if not ("'" in obj_atom.name or "P" in obj_atom.name or RE_TERMH1.match(obj_atom.name) or RE_TERMH2.match(obj_atom.name))]
					fragment_idx += 1
					obj_fragment_base = FragmentData(fragment_idx)
					obj_fragment_base.set_fragment_index(fragment_idx)
					obj_fragment_base.set_charge(-1)
					obj_fragment_base.set_atoms(atoms_base)
					self._fragments.append(obj_fragment_base)

					self._connections.append([])



				elif option_nuc == "+base":
					if len(self._fragments) >= 1:
						# move phosphate atoms to the previous backbone fragment
						self._fragments[-1].set_atoms(list(sorted(self._fragments[-1].atoms + new_atoms_remain)))

					atoms = [obj_atom.number for obj_atom in obj_residue.atoms]
					fragment_idx += 1
					obj_fragment = FragmentData()
					obj_fragment.set_fragment_index(fragment_idx)
					obj_fragment.set_atoms(atoms)
					if obj_residue.name.endswith("5"):
						obj_fragment.set_charge(0)
					else:
						obj_fragment.set_charge(-1)
					self._fragments.append(obj_fragment)
					self._connections

				else:
					sys.stderr.write("ERROR: undefined option_nuc.\n")
					sys.exit(1)

			elif residue_type == "Water":
				obj_fragment.set_charge(0)

			elif residue_type == "Ion":
				if self.residue_name in ["Na", "K"]:
					obj_fragment.set_charge(1)
				elif self.residue_name in ["Zn", "Ca", "CAL", "ZIN"]:
					obj_fragment.set_charge(2)

			else:
				obj_fragment.set_charge("ERR")


			# overwrite charge information
			if obj_residue.name in charge_map.keys():
				obj_fragment.set_charge(charge_map[obj_residue.name])

		self._n_atom = sum([len(obj_fragment.atoms) for obj_fragment in self._fragments])
		self._charge = sum([obj_fragment.charge for obj_fragment in self._fragments if obj_fragment.charge != "ERR"])

		return self._fragments
