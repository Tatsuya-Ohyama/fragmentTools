#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
from mods.FragmentData import FragmentData



# =============== variables =============== #
RE_FRAGMENT = re.compile(r"^[\s\t]*\d+[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|(?:[\s\t]*\d+)+")
RE_CONNECTION = re.compile(r"^(?:[\s\t]*\d+){2}[\s\t]*$")



# =============== classes =============== #
class FredData:
	""" fred File class """
	def __init__(self, input_file):
		self._charge = 0
		self._atom = 0
		self._fragments = []
		self._connection = []
		self._others = []

		self._load_file(input_file)


	@property
	def fragment_number(self):
		return len(self._fragments)

	@property
	def last_fragment_index(self):
		return self._fragments[-1].fragment_index


	def _load_file(self, input_file):
		"""
		Method to read .fred file

		Args:
			input_file (str): .fred file
		"""
		with open(input_file, "r") as obj_input:
			cnt_line = 0
			flag_read = 0
			for line in obj_input:
				cnt_line += 1
				line = line.rstrip("\r\n")

				if cnt_line == 2:
					flag_read = 1

				elif "connections" in line:
					flag_read = 2

				elif "namelist" in line:
					flag_read = 3

				if flag_read == 1 and RE_FRAGMENT.search(line):
					fragment = FragmentData(line)
					if fragment.charge != "ERR":
						self._charge += fragment.charge
					self._fragments.append(fragment)
					self._atom += len(fragment.atoms)

				elif flag_read == 2 and RE_CONNECTION.search(line):
					self._connection.append([int(x) for x in line.strip().split()])

				elif flag_read == 3:
					self._others.append(line)


	def add_fragment(self, obj_fragment):
		"""
		Method to add fragment

		Args:
			obj_fragment (FragmentData object): FragmentData object

		Returns:
			self
		"""
		for idx, fragment in enumerate(self._fragments):
			fragment.delete_duplicate_atom(obj_fragment.atoms)
			if len(fragment.atoms) == 0:
				self._fragments.pop(idx)

		self._atom += len(obj_fragment.atoms)
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


	def write_file(self, output_file):
		"""
		Method to write out file

		Args:
			output_file (str): output file
		"""
		re_NF = re.compile(r"NF=[\s\t]*-?\d+")
		tmp_fragments = sorted([[obj_fragment.min_idx, obj_fragment] for obj_fragment in self._fragments], key=lambda x : x[0])
		self._fragments = [obj_fragment[1].update_fragment_index(idx + 1) for idx, obj_fragment in enumerate(tmp_fragments)]
		tmp_connection_int = sorted([x for x in self._connection if isinstance(x[0], int)], key=lambda x : x[0])
		tmp_connection_str = sorted([x for x in self._connection if isinstance(x[0], str)], key=lambda x : x[0])
		self._connection = tmp_connection_int + tmp_connection_str

		with open(output_file, "w") as obj_output:
			obj_output.write("  FNo.  | Charge | BDA | Atoms of fragment\n")
			for idx, fragment in enumerate(self._fragments):
				fragment.update_fragment_index(idx + 1)
				obj_output.write("{0:>7} |{1:>6}  |{2:>3}  |{3}\n".format(fragment.fragment_index, fragment.charge, fragment.bda, " ".join(["{0:>8}".format(x) for x in fragment.atoms])))
			obj_output.write("\n")

			obj_output.write("<< connections (ex. \"Next_fragment_atom   Prev_fragment_atom\") >>\n")
			for connection in self._connection:
				obj_output.write("{0}\n".format(" ".join(["{0:>9}".format(x) for x in connection])))
			obj_output.write("\n")

			for line in self._others:
				redb_NF = re_NF.search(line)
				if redb_NF:
					line = line[: redb_NF.start()] + "NF={0}".format(len(self._fragments)) + line[redb_NF.end() :]
				obj_output.write(line + "\n")
