#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import parmed

from mods.FragmentData import FragmentData



# =============== classes =============== #
class MoleculeInformation:
	""" Molecular structure class """
	def __init__(self, input_file):
		# member
		self._obj_mol = None
		self._filename = None

		# initiation
		self._load_file(input_file)


	@property
	def obj_mol(self):
		return self._obj_mol

	@property
	def filename(self):
		return self._filename


	def _load_file(self, input_file):
		"""
		Method to read file

		Args:
			input_file (str): molecular file
		"""
		self._filename = input_file
		self._obj_mol = parmed.load_file(input_file)
		return self


	def get_info(self, data_type=None, idx=None):
		"""
		Method to return atomic information

		Args:
			data_type (str, optional): `residue_name`, `residue_idx`, `atom_name` or `atom_idx` (Default: None)
			idx (int, optional): target index (Default: None)

		Returns:
			list
		"""
		if data_type == "residue_name":
			if idx is not None:
				try:
					return self._obj_mol.atoms[idx].residue.name
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [obj_atom.residue.name for obj_atom in self._obj_mol.atoms]

		elif data_type == "residue_idx":
			if idx is not None:
				try:
					return self._obj_mol.atoms[idx].residue.number
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [obj_atom.residue.number for obj_atom in self._obj_mol.atoms]

		elif data_type == "atom_name":
			if idx is not None:
				try:
					return self._obj_mol.atoms[idx].name
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [obj_atom.name for obj_atom in self._obj_mol.atoms]

		elif data_type == "atom_idx":
			if idx is not None:
				try:
					return self._obj_mol.atoms[idx].number
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [obj_atom.number for obj_atom in self._obj_mol.atoms]

		else:
			if idx is not None:
				try:
					return [[res_name, res_idx, atom_name, atom_idx] for res_name, res_idx, atom_name, atom_idx in zip(self.get_info("residue_name"), self.get_info("residue_idx"), self.get_info("atom_name"), self.get_info("atom_idx"))][idx]
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [[res_name, res_idx, atom_name, atom_idx] for res_name, res_idx, atom_name, atom_idx in zip(self.get_info("residue_name"), self.get_info("residue_idx"), self.get_info("atom_name"), self.get_info("atom_idx"))]


	def convert_number(self, mask):
		"""
		Method to convert Ambermask to atom index

		Args:
			mask (str): Ambermask

		Returns:
			list
		"""
		mask = parmed.amber.AmberMask(self._obj_mol, mask)
		return [self.get_info("atom_idx", idx) for idx in mask.Selected()]


	def get_coord(self, idx=None):
		"""
		Method to return coordinates

		Args:
			idx (int, optional): atom index (Default: None)

		Returns:
			list
		"""
		if idx is not None:
			try:
				return [[obj_atom.xx, obj_atom.xy, obj_atom.xz] for obj_atom in self._obj_mol.atoms[idx]]
			except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
		else:
			return [[obj_atom.xx, obj_atom.xy, obj_atom.xz] for obj_atom in self._obj_mol.atoms]


	def output_fragmentdata(self, output_type, obj_base_molecule, flag_multi):
		"""
		Method to output fragment information

		Args:
			output_type (str): `text` or `object`
			obj_base_molecule (MoleculeInformation object): MoleculeInformation object
			flag_multi (bool): apply the same fragmentation to multiple molecules

		Returns:
			list
		"""
		fragment_records = [[]]

		if flag_multi:
			# apply the same fragmentation to multiple molecules
			mol_info = [list(x) for x in list(zip(self.get_info("residue_name"), self.get_info("atom_name")))]
			len_atom = len(self.get_info("atom_name"))
			check_list = [False] * len_atom

			for info in obj_base_molecule.get_info():
				# loop for whole structure
				if [info[0], info[2]] in mol_info:
					# when contain residue name and atom name
					match_idx = mol_info.index([info[0], info[2]])

					if check_list[match_idx]:
						# when dealing with the same atoms as before, create new fragment list
						fragment_records.append([])
						check_list = [False] * len_atom

					fragment_records[-1].append(info[3])
					check_list[match_idx] = True

		else:
			# when apply as single fragment, match by coordinates
			coords_base = obj_base_molecule.get_coord()

			for array_idx, coords in enumerate(self.get_coord()):
				# loop for coordinates and index
				if coords in coords_base:
					# when match with coordinates, get index
					match_idx = coords_base.index(coords)
					if self.get_info("atom_name", array_idx) == obj_base_molecule.get_info("atom_name", match_idx):
						# also check for atomic name, and add fragment information
						fragment_records[-1].append(obj_base_molecule.get_info("atom_idx", match_idx))
					else:
						sys.stderr.write("ERROR: wrong coordinates in MoleculeInformation.output_fragmentdata(). Perhaps mismatched atom names in `{0}`.\n".format(self._filename))
						sys.exit(1)

		if output_type == "text":
			return ["{0}|{1}|{2}|{3}".format(0, "ERR", "ERR", " ".join([str(v) for v in record])) for record in fragment_records]
		elif output_type == "object":
			return [FragmentData().set_fragment_index(0).set_charge(None).set_bda(None).set_atoms(record) for record in fragment_records]
		else:
			sys.stderr.write("ERROR: undefined output_type in MoleculeInformation class.\n")
			sys.exit(1)

		return fragment_records
