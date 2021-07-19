#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from mods.FragmentData import FragmentData



# =============== constant =============== #
INDENT = "  "



# =============== function =============== #
def split_n(line, length):
	"""
	Function to split by n-chars

	Args:
		line (str): target string
		length (int): split length

	Returns:
		list
	"""
	line = line.rstrip("\r\n")
	datas = []
	pos = 0
	while pos + length <= len(line):
		datas.append(line[pos : pos + length])
		pos += length

	if pos != len(line):
		datas.append(line[pos:len(line)])

	datas = list(map(lambda data:data.strip(), datas))

	return datas



# =============== class =============== #
class FileAJF:
	""" .ajf file class """
	def __init__(self, ajf_file):
		# member
		self._parameters = {}

		# init
		self.read(ajf_file)


	@property
	def parameters(self):
		return self._parameters


	def read(self, input_file):
		"""
		function to read .ajf file

		Args:
			ajf_file (str): .ajf file (input)

		Returns:
			self
		"""
		self._parameters = {"LIST_ORDER": []}
		is_list_type = False
		with open(input_file, "r") as obj_input:
			group_name = None
			for line_val in obj_input:
				line_orig = line_val
				line_val = line_val.strip()
				if len(line_val.strip()) == 0:
					continue

				if line_val.startswith("&"):
					# start of name list
					group_name = line_val
					self._parameters["LIST_ORDER"].append(group_name)

					if group_name in ["&XYZ", "&FRAGMENT", "&FRAGPAIR"]:
						is_list_type = True
						self._parameters[group_name] = []
					else:
						is_list_type = False
						self._parameters[group_name] = {}

				elif line_val.startswith("/"):
					# end of name list
					group_name = None
					is_list_type = False
					continue

				elif is_list_type:
					# in name list of list type (order is important)
					self._parameters[group_name].append(line_orig)

				else:
					# in namelist of non-list type (order is not important)
					tmp_values = [v.strip() for v in line_val.split("=", maxsplit = 1)]
					self._parameters[group_name][tmp_values[0]] = tmp_values[1]

		return self


	def set_parameters(self, parameters):
		"""
		Method to set parameters

		Args:
			parameters (dict): parameters

		Returns:
			self
		"""
		self._parameters = parameters
		return self


	def create_fragment_objects(self, pdb_file=None):
		"""
		Method to create fragment objects

		Args:
			pdb_file (str, optional): .pdb file (Default: None)

		Returns:
			list: fragment object list
			list: connection list
		"""
		# get number of atoms
		if pdb_file is None:
			if "&CNTRL" not in self.parameters["LIST_ORDER"]:
				sys.stderr.write("ERROR: .ajf file has no `&CNTRL` namelist.\n")
				sys.exit(1)

			if "ReadGeom" not in self.parameters["&CNTRL"].keys():
				sys.stderr.write("ERROR: .ajf file has no `ReadGeom` parameter.\n")
				sys.exit(1)

			pdb_file = self.parameters["&CNTRL"]["ReadGeom"]
			if pdb_file.startswith("'") and pdb_file.endswith("'") \
				or pdb_file.startswith('"') and pdb_file.endswith('"'):
				pdb_file = pdb_file[1:-1]

		n_atom = 0
		with open(pdb_file, "r") as obj_input:
			for line_val in obj_input:
				if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
					n_atom += 1

		# deal with fragment information
		if "&FRAGMENT" not in self.parameters["LIST_ORDER"]:
			sys.stderr.write("ERROR: .ajf file has not `FRAGMENT` namelist.\n")
			sys.exit(1)

		list_fragments = []
		list_connections = []
		fragment_idx_start = 0
		n_atom_curr = 0
		atom_idx = 0
		n_fragment = 0
		n_fragment_curr = 0
		fragment_idx = 0
		n_step = 0
		for line_val in self.parameters["&FRAGMENT"]:
			datas = [int(v) for v in split_n(line_val, 8)]

			if n_step == 0:
				# number of atoms in fragment
				if "0" in datas:
					# when fragment has 0 atoms
					sys.stderr.write("ERROR: zero atoms in fragment found.\n")
					sys.exit(1)

				list_fragments += [
					FragmentData().set_fragment_index(i).set_atoms([None for _ in range(v)])
					for i, v in enumerate(datas, fragment_idx_start)
				]
				fragment_idx_start += len(datas)
				n_atom_curr += sum(datas)
				if n_atom_curr == n_atom:
					# match number of atoms in system
					n_step = 1
					n_fragment = len(list_fragments)
				elif n_atom < n_atom_curr:
					# over number of atoms
					sys.stderr.write("ERROR: The number of atoms are different between .ajf ({0}) and .pdb {1} files.\n".format(
						n_atom_curr,
						n_atom
					))
					sys.exit(1)

			elif n_step == 1:
				# charge information
				for obj_fragment, charge in zip(list_fragments[n_fragment_curr : n_fragment_curr + len(datas)], datas):
					obj_fragment.set_charge(charge)
				n_fragment_curr += len(datas)

				if n_fragment_curr == n_fragment:
					# match number of fragments in system
					n_step = 2
					n_fragment_curr = 0

				elif n_fragment < n_fragment_curr:
					# over number of fragments
					sys.stderr.write("ERROR: The number of fragments are different ({0}, {1}) at {2} section.\n".format(
						n_fragment_curr,
						n_fragment,
						"charge"
					))
					sys.exit(1)

			elif n_step == 2:
				# BDA information
				for obj_fragment, bda in zip(list_fragments[n_fragment_curr : n_fragment_curr + len(datas)], datas):
					obj_fragment.set_bda(bda)
				n_fragment_curr += len(datas)

				if n_fragment_curr == n_fragment:
					# match number of fragments in system
					n_step = 3
					n_fragment_curr = 0

				elif n_fragment < n_fragment_curr:
					# over number of fragments
					sys.stderr.write("ERROR: The number of fragments are different ({0}, {1}) at {2} section.\n".format(
						n_fragment_curr,
						n_fragment,
						"BDA"
					))
					sys.exit(1)

			elif n_step == 3:
				# member of atoms in fragment
				obj_fragment = list_fragments[fragment_idx]
				list_atoms = obj_fragment.atoms
				list_atoms[atom_idx:atom_idx+len(datas)] = datas
				obj_fragment.set_atoms(list_atoms)
				if obj_fragment.atoms.count(None) == 0:
					atom_idx = 0
					fragment_idx += 1
				else:
					atom_idx += len(datas)

				if fragment_idx == n_fragment:
					# match number of fragments in system
					n_step = 4
					n_fragment_curr = 0

				elif fragment_idx > n_fragment:
					# over number of fragments
					sys.stderr.write("ERROR: The number of fragments are different ({0}, {1}) at {2} section.\n".format(
						n_fragment_curr,
						n_fragment,
						"atoms of fragment"
					))
					sys.exit(1)

			elif n_step == 4:
				# connection information
				list_connections.append(datas)

		return list_fragments, list_connections


	def parameter_list(self, omit_fragment=False):
		"""
		Method to output parameter list (not dict)

		Args:
			omit_fragment (bool, optional): omit &FRAGMENT information (Default: False)

		Returns:
			list
		"""
		parameters = []
		for group_name in self._parameters["LIST_ORDER"]:
			parameters.append("{0}\n".format(group_name))
			if group_name == "&FRAGMENT" and omit_fragment:
					parameters.append("{...}\n")
			elif isinstance(self._parameters[group_name], dict):
				for parameter_name, parameter_value in self._parameters[group_name].items():
					parameters.append("{0}{1}={2}\n".format(
						INDENT,
						parameter_name,
						parameter_value
					))
			else:
				for line_val in self._parameters[group_name]:
					parameters.append(line_val)
			parameters.append("/\n\n")

		return parameters


	def write(self, output_file):
		"""
		Method to write .ajf file

		Args:
			output_file (str): .ajf file

		Returns:
			self
		"""
		with open(output_file, "w") as obj_output:
			for group_name in self._parameters["LIST_ORDER"]:
				obj_output.write("{0}\n".format(group_name))
				if isinstance(self._parameters[group_name], dict):
					for parameter_name, parameter_value in self._parameters[group_name].items():
						obj_output.write("{0}{1}={2}\n".format(
							INDENT,
							parameter_name,
							parameter_value
						))
				else:
					for line_val in self._parameters[group_name]:
						obj_output.write(line_val)
				obj_output.write("/\n\n")
