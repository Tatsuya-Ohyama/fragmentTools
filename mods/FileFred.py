#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import re
import copy

from mods.FragmentData import FragmentData



# =============== variables =============== #
RE_FRAGMENT = re.compile(r"^[\s\t]*\d+[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|(?:[\s\t]*\d+)+")
RE_CONNECTION = re.compile(r"^(?:[\s\t]*\d+){2}[\s\t]*$")
RE_INT = re.compile(r"^-?\d+$")
INDENT = "  "



# =============== classes =============== #
class FileFred:
	""" fred File class """
	def __init__(self):
		self._n_atom = None
		self._charge = None
		self._fragments = []
		self._parameters = {}


	@property
	def n_atom(self):
		if self._n_atom is None:
			return sum([len(obj_fragment.atoms) for obj_fragment in self._fragments])
		return self._n_atom

	@property
	def charge(self):
		if self._charge is None:
			return sum([obj_fragment.charge for obj_fragment in self._fragments if obj_fragment.charge is not None])
		return self._charge

	@property
	def fragments(self):
		return self._fragments

	@property
	def n_fragment(self):
		return len(self._fragments)

	@property
	def parameters(self):
		return self._parameters

	@property
	def complete_parameters(self):
		parameters = copy.deepcopy(self._parameters)
		if "&CNTRL" in parameters["LIST_ORDER"]:
			parameters["&CNTRL"]["Natom"] = sum([len(v.atoms) for v in self._fragments])
			parameters["&CNTRL"]["Charge"] = sum([v.charge for v in self._fragments])

		if "&FMOCNTRL" in parameters["LIST_ORDER"]:
			parameters["&FMOCNTRL"]["NF"] = len(self._fragments)
			parameters["&FMOCNTRL"]["AutoFrag"] = "'OFF'"

		parameters["&FRAGMENT"] = ""
		parameters["&FRAGMENT"] += "\n".join([
			"".join(
				["{0:>8}".format(len(obj_fragment.atoms)) for obj_fragment in self._fragments[i:i+10]]
			) for i in range(0, len(self._fragments), 10)
		]) + "\n"
		parameters["&FRAGMENT"] += "\n".join([
			"".join(
				["{0:>8}".format(obj_fragment.charge) for obj_fragment in self._fragments[i:i+10]]
			) for i in range(0, len(self._fragments), 10)
		]) + "\n"
		parameters["&FRAGMENT"] += "\n".join([
			"".join(
				["{0:>8}".format(obj_fragment.bda) for obj_fragment in self._fragments[i:i+10]]
			) for i in range(0, len(self._fragments), 10)
		]) + "\n"
		parameters["&FRAGMENT"] += "\n".join([
			"\n".join([
				"".join([
					"{0:>8}".format(atom_idx) for atom_idx in obj_fragment.atoms[atom_i:atom_i+10]
				]) for atom_i in range(0, len(obj_fragment.atoms), 10)
			]) for obj_fragment in self._fragments
		]) + "\n"
		parameters["&FRAGMENT"] += "\n".join([
			"{0[0]:>8}{0[1]:>8}".format(connection)
			for obj_fragment in [obj_fragment for obj_fragment in self._fragments if len(obj_fragment.get_connections()) != 0]
			for connection in obj_fragment.get_connections()
		]) + "\n"
		return parameters


	def set_n_atom(self, n_atom):
		"""
		Method to set number of atoms in system

		Args:
			n_atom (int): number of atoms in system

		Returns:
			self
		"""
		self._n_atom = n_atom
		return self


	def set_charge(self, charge):
		"""
		Method to set system charge

		Args:
			charge (int): system charge

		Returns:
			self
		"""
		self._charge = charge
		return self


	def set_fragments(self, list_fragments):
		"""
		Method to set fragment object list

		Args:
			list_fragments (list): fragment object list

		Returns:
			self
		"""
		self._fragments = list_fragments
		return self


	def set_parameters(self, parameters):
		"""
		Method to set other information

		Args:
			parameters (list): other information

		Returns:
			self
		"""
		self._parameters = parameters
		return self


	def read(self, input_file):
		"""
		Method to read .fred file

		Args:
			input_file (str): .fred file

		Returns:
			self
		"""
		with open(input_file, "r") as obj_input:
			flag_read = 0
			is_list_type = False
			group_name = None
			self._parameters = {"LIST_ORDER": []}
			fragment_index = 1
			for line_idx, line_val in enumerate(obj_input, 1):
				line_val = line_val.strip()

				if line_idx == 1:
					flag_read = 1

				elif "connections" in line_val.lower():
					flag_read = 2

				elif "< namelist >" in line_val.lower():
					flag_read = 3

				elif flag_read == 1 and RE_FRAGMENT.search(line_val):
					elems = [v.strip() for v in line_val.strip().split("|")]

					obj_fragment = FragmentData()
					obj_fragment.set_fragment_index(fragment_index)
					fragment_index += 1

					# Charge
					if RE_INT.search(elems[1]):
						elems[1] = int(elems[1])
						self._charge = self.charge + elems[1]
					else:
						elems[1] = None
					obj_fragment.set_charge(elems[1])

					# BDA
					if RE_INT.search(elems[2]):
						elems[2] = int(elems[2])
					else:
						elems[2] = None
					obj_fragment.set_bda(elems[2])

					# atoms
					obj_fragment.set_atoms([int(v) for v in elems[3].split()])

					self._fragments.append(obj_fragment)
					self._n_atom = self.n_atom + len(obj_fragment.atoms)

				elif flag_read == 2 and RE_CONNECTION.search(line_val):
					# connection
					atom_i, atom_j = [int(x) for x in line_val.strip().split(maxsplit=1)]
					obj_fragment_j = [obj_fragment for obj_fragment in self._fragments if atom_j in obj_fragment.atoms][0]
					obj_fragment_j.append_connection([atom_i, atom_j])

				elif flag_read == 3:
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


	def add_fragment(self, obj_fragment):
		"""
		Method to add fragment and move exsited fragment to new fragment

		Args:
			obj_fragment (FragmentData object): FragmentData object

		Returns:
			self
		"""
		atoms = set(obj_fragment.atoms)
		new_list_fragment = []
		for fragment_idx, obj_fragment_registered in enumerate(self._fragments):
			atoms_registered = set(obj_fragment_registered.atoms)
			if len(atoms & atoms_registered) != 0:
				new_atoms = sorted(list(atoms_registered - atoms))

				if len(new_atoms) == 0:
					continue

				obj_fragment_registered.set_atoms(new_atoms)
			new_list_fragment.append(obj_fragment_registered)
		new_list_fragment.append(obj_fragment)

		self._fragments = new_list_fragment
		return self


	def check_fragments(self, verbose=False):
		"""
		Method to check fragments (charge, BDA and atoms)

		Args:
			verbose (bool): True=output to console

		Returns:
			bool
		"""
		has_problem_whole = False

		for i in range(len(self._fragments)):
			obj_fragment = self._fragments[i]
			has_problem = False

			conflict_fragments = []
			for j in range(len(self._fragments)):
				if i == j:
					continue
				obj_fragment_other = self._fragments[j]
				common_atom = set(obj_fragment.atoms) & set(obj_fragment_other.atoms)
				if len(common_atom) != 0:
					conflict_fragments.append(obj_fragment_other.fragment_index)
			if len(conflict_fragments) != 0:
				has_problem_whole = True
				if verbose:
					sys.stderr.write("WARNING: Fragment {0} has problem:\n".format(obj_fragment.fragment_index))
					sys.stderr.write("        * Conflicted atoms with Fragment {0}\n".format(", ".join([str(v) for v in conflict_fragments])))
				has_problem = True

			if obj_fragment.charge is None:
				has_problem_whole = True
				if verbose:
					if has_problem:
						sys.stderr.write("WARNING: Fragment {0} has problem:\n".format(obj_fragment.fragment_index))
					sys.stderr.write("        * No charge information.\n")
				has_problem = True

			if obj_fragment.bda is None:
				has_problem_whole = True
				if verbose:
					if has_problem:
						sys.stderr.write("WARNING: Fragment {0} has problem:\n".format(obj_fragment.fragment_index))
					sys.stderr.write("        * No BDA information.\n")

		return has_problem_whole


	def write(self, output_file, indent=INDENT):
		"""
		Method to write out file

		Args:
			output_file (str): output file
			indent (int, optional): indent length

		Returns:
			self
		"""
		self.check_fragments(verbose=True)
		with open(output_file, "w") as obj_output:
			obj_output.write("  FNo.  | Charge | BDA | Atoms of fragment\n")
			for obj_fragment in self._fragments:
				obj_output.write("{0:>7} |{1:>5}   |{2:>3}  | {3}\n".format(
					obj_fragment.fragment_index,
					obj_fragment.charge if obj_fragment.charge is not None else "ERR",
					obj_fragment.bda if obj_fragment.bda is not None else "ERR",
					" ".join(["{0:>8}".format(x) for x in obj_fragment.atoms])
				))
			obj_output.write("\n")

			obj_output.write("<< connections (ex. \"Next_fragment_atom   Prev_fragment_atom\") >>\n")
			for obj_fragment in self._fragments:
				list_connections = obj_fragment.get_connections()
				if len(list_connections) == 0:
					continue
				for connection in list_connections:
					obj_output.write("{0[0]:>9} {0[1]:>9}\n".format(connection))
			obj_output.write("\n")

			obj_output.write("===============< namelist >===============\n")

			if "&CNTRL" in self._parameters["LIST_ORDER"]:
				self._parameters["&CNTRL"]["Natom"] = self.n_atom
				self._parameters["&CNTRL"]["Charge"] = self.charge

			if "&FMOCNTRL" in self._parameters["LIST_ORDER"]:
				self._parameters["&FMOCNTRL"]["NF"] = len(self._fragments)

			for group_name in self._parameters["LIST_ORDER"]:
				obj_output.write("{0}\n".format(group_name))
				if isinstance(self._parameters[group_name], dict):
					for parameter_name, parameter_value in self._parameters[group_name].items():
						obj_output.write("{0}{1}={2}\n".format(
							indent,
							parameter_name,
							parameter_value
						))
				else:
					for line_val in self._parameters[group_name]:
						obj_output.write(line_val + "\n")
				obj_output.write("/\n\n")

		return self
