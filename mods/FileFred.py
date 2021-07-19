#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
from mods.FragmentData import FragmentData



# =============== variables =============== #
RE_FRAGMENT = re.compile(r"^[\s\t]*\d+[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|(?:[\s\t]*\d+)+")
RE_CONNECTION = re.compile(r"^(?:[\s\t]*\d+){2}[\s\t]*$")



# =============== classes =============== #
class FileFred:
	""" fred File class """
	def __init__(self):
		self._n_atom = 0
		self._charge = 0
		self._fragments = []
		self._connection = []
		self._others = []


	@property
	def n_atom(self):
		return self._n_atom

	@property
	def charge(self):
		return self._charge

	@property
	def fragments(self):
		return self._fragments

	@property
	def n_fragment(self):
		return len(self._fragments)


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


	def set_connections(self, list_connection):
		"""
		Method to set connection list

		Args:
			list_connection (list): connection list

		Returns:
			self
		"""
		self._connection = list_connection
		return self


	def set_other_info(self, other_info):
		"""
		Method to set other information

		Args:
			other_info (list): other information

		Returns:
			self
		"""
		self._others = other_info
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
			for line_idx, line_val in enumerate(obj_input, 1):
				line_val = line_val.strip()

				if line_idx == 2:
					flag_read = 1

				elif "connections" in line_val.lower():
					flag_read = 2

				elif "< namelist >" in line_val.lower():
					flag_read = 3

				if flag_read == 1 and RE_FRAGMENT.search(line_val):
					obj_fragment = FragmentData(line_val)
					if fragment.charge != "ERR":
						self._charge += obj_fragment.charge
					self._fragments.append(obj_fragment)
					self._n_atom += len(obj_fragment.atoms)

				elif flag_read == 2 and RE_CONNECTION.search(line_val):
					self._connection.append([int(x) for x in line_val.strip().split()])

				elif flag_read == 3:
					self._others.append(line_val)
		return self


	def add_fragment(self, obj_fragment):
		"""
		Method to add fragment

		Args:
			obj_fragment (FragmentData object): FragmentData object

		Returns:
			self
		"""
		for fragment_idx, obj_fragment_registered in enumerate(self._fragments):
			remain_atoms = sorted(list(set(obj_fragment_registered.atoms) - set(obj_fragment.atoms)))
			obj_fragment_registered.set_atoms(remain_atoms)
			if len(obj_fragment_registered.atoms) == 0:
				self._fragments.pop(fragment_idx)

		self._n_atom += len(obj_fragment.atoms)
		self._fragments.append(obj_fragment)
		return self


	def add_connection(self, connection):
		"""
		Method to add fragment connection

		Args:
			connection (int or list): connection information

		Returns:
			self
		"""
		if isinstance(connection, int):
			for i in range(connection):
				self._connection.append(["*", "*"])
		elif isinstance(connection, list):
			self._connection.append(connection)
		return self


	def write(self, output_file):
		"""
		Method to write out file

		Args:
			output_file (str): output file

		Returns:
			self
		"""
		re_NF = re.compile(r"NF=[\s\t]*-?\d+")
		tmp_fragments = sorted([[obj_fragment.min_index, obj_fragment] for obj_fragment in self._fragments], key=lambda x : x[0])
		self._fragments = [obj_fragment[1].set_fragment_index(idx) for idx, obj_fragment in enumerate(tmp_fragments, 1)]
		tmp_connection_int = sorted([x for x in self._connection if isinstance(x[0], int)], key=lambda x : x[0])
		tmp_connection_str = sorted([x for x in self._connection if isinstance(x[0], str)], key=lambda x : x[0])
		self._connection = tmp_connection_int + tmp_connection_str

		with open(output_file, "w") as obj_output:
			obj_output.write("  FNo.  | Charge | BDA | Atoms of fragment\n")
			for idx, fragment in enumerate(self._fragments, 1):
				fragment.set_fragment_index(idx)
				obj_output.write("{0:>7} |{1:>6}  |{2:>3}  |{3}\n".format(fragment.fragment_index, fragment.charge, fragment.bda, " ".join(["{0:>8}".format(x) for x in fragment.atoms])))
			obj_output.write("\n")

			obj_output.write("<< connections (ex. \"Next_fragment_atom   Prev_fragment_atom\") >>\n")
			for connection in self._connection:
				obj_output.write("{0}\n".format(" ".join(["{0:>9}".format(x) for x in connection])))
			obj_output.write("\n")

			obj_output.write("===============< namelist >===============\n")
			for line_val in self._others:
				redb_NF = re_NF.search(line_val)
				if redb_NF:
					line_val = line_val[: redb_NF.start()] + "NF={0}".format(len(self._fragments)) + line_val[redb_NF.end() :]
				obj_output.write(line_val)

		return self
