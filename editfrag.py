#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
editfrag.py
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from py_module_basic import basic


# =============== class =============== #
class MoleculeInformation:
	""" 分子構造のクラス """
	def __init__(self, input_file):
		self.__atom_names = []
		self.__atom_idxs = []
		self.__coords = []
		self.__load_file(input_file)

	def __load_file(self, input_file):
		""" ファイルを読み込むメソッド """
		re_atom = re.compile(r"^(?:(?:ATOM)|(?:HETATM))")
		with open(input_file, "r") as obj_input:
			atom_idx = 0
			for line in obj_input:
				if re_atom.search(line):
					line = line.rstrip("\r\n")
					atom_idx += 1
					self.__atom_names.append(line[12:16])
					self.__atom_idxs.append(int(line[6:11]))
					self.__coords.append(line[30:54])

	def get_atom_index(self, idx = None):
		""" 原子順序番号を出力するメソッド """
		if idx is None:
			return self.__atom_idxs
		else:
			try:
				return self.__atom_idxs[idx]
			except IndexError:
				sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_atom_index\n")
				sys.exit(1)

	def get_atom_name(self, idx = None):
		""" 原子名を出力するメソッド """
		if idx is None:
			return self.__atom_names
		else:
			try:
				return self.__atom_names[idx]
			except IndexError:
				sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_atom_name\n")
				sys.exit(1)

	def get_coord(self, idx = None):
		""" 座標を出力するメソッド """
		if idx is None:
			return self.__coords
		else:
			try:
				return self.__coords[idx]
			except IndexError:
				sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_coord\n")
				sys.exit(1)

	def update_index(self, obj_base_molecule):
		""" 原子順序番号を更新するメソッド """
		coords_base = obj_base_molecule.get_coord()

		for array_idx, coords in enumerate(self.__coords):
			if coords in coords_base:
				match_idx = coords_base.index(coords)
				if self.__atom_names[array_idx] == obj_base_molecule.get_atom_name(match_idx):
					self.__atom_idxs[array_idx] = obj_base_molecule.get_atom_index(match_idx)
				else:
					sys.stderr.write("ERROR: wrong coordinates in MoleculeInformation.update_index()\n")
					sys.exit(1)

	def output_fragmentdata(self):
		return "{0}|{1}|{2}|{3}".format(0, "ERR", "ERR", " ".join([str(x) for x in self.__atom_idxs]))



class FragmentData:
	""" フラグメントデータのクラス """
	def __init__(self, line):
		line = line.strip("\r\n")
		datas = [x.strip() for x in line.split("|")]

		re_int = re.compile(r"^-?\d+$")
		self.__idx = int(datas[0])
		self.__charge = 0
		if re_int.search(datas[1]):
			self.__charge = int(datas[1])
		else:
			self.__charge = "ERR"

		self.__bda = 0
		if re_int.search(datas[2]):
			self.__bda = int(datas[2])
		else:
			self.__bda = "ERR"

		self.__atoms = sorted([int(x) for x in datas[3].split()])

	def update_fragment_index(self, index):
		""" フラグメント番号を更新するメソッド """
		self.__idx = index

	def duplicate_atom(self, atoms):
		""" フラグメント構成原子の重複を削除するメソッド """
		self.__atoms = sorted(list(set(self.__atoms) - set(atoms)))

	def get_fragment_index(self):
		""" フラグメント番号を返すメソッド """
		return self.__idx

	def get_charge(self):
		""" 電荷を返すメソッド """
		return self.__charge

	def get_bda(self):
		""" BDA を返すメソッド """
		return self.__bda

	def get_atoms(self):
		""" 構成原子を返すメソッド """
		return self.__atoms



class FredData:
	""" fred ファイルのクラス """
	def __init__(self, input_file):
		self.__charge = 0
		self.__atom = 0
		self.__fragments = []
		self.__connection = []
		self.__others = []

		self.__load_file(input_file)

	def __load_file(self, input_file):
		""" fred データを読み込むメソッド """
		with open(input_file, "r") as obj_input:
			cnt_line = 0
			re_fragment = re.compile(r"^[\s\t]*\d+[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|(?:[\s\t]*\d+)+")
			re_connection = re.compile(r"^(?:[\s\t]*\d+){2}[\s\t]*$")
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

				if flag_read == 1 and re_fragment.search(line):
					fragment = FragmentData(line)
					if fragment.get_charge() != "ERR":
						self.__charge += fragment.get_charge()
					self.__fragments.append(fragment)
					self.__atom += len(fragment.get_atoms())

				elif flag_read == 2 and re_connection.search(line):
					self.__connection.append([int(x) for x in line.strip().split()])

				elif flag_read == 3:
					self.__others.append(line)

	def add_fragment(self, obj_fragment):
		""" フラグメントを追加するメソッド """
		for idx, fragment in enumerate(self.__fragments):
			fragment.duplicate_atom(obj_fragment.get_atoms())
			if len(fragment.get_atoms()) == 0:
				self.__fragments.pop(idx)

		obj_fragment.update_fragment_index(self.__fragments[-1].get_fragment_index() + 1)

		self.__atom += len(obj_fragment.get_atoms())
		self.__fragments.append(obj_fragment)
		self.__connection.append(["*", "*"])

	def get_fragment_number(self):
		""" フラグメント数を返すメソッド """
		return len(self.__fragments)

	def get_last_fragment_index(self):
		""" 最後のフラグメント番号を返すメソッド """
		return self.__fragments[-1].get_fragment_index()

	def write_file(self, output_file):
		""" ファイルに出力するメソッド """
		re_NF = re.compile(r"NF=[\s\t]*-?\d+")

		with open(output_file, "w") as obj_output:
			obj_output.write("  FNo.  | Charge | BDA | Atoms of fragment\n")
			for idx, fragment in enumerate(self.__fragments):
				fragment.update_fragment_index(idx + 1)
				obj_output.write("{0:>7} |{1:>6}  |{2:>3}  |{3}\n".format(fragment.get_fragment_index(), fragment.get_charge(), fragment.get_bda(), " ".join(["{0:>8}".format(x) for x in fragment.get_atoms()])))
			obj_output.write("\n")

			obj_output.write("<< connections (ex. \"Next_fragment_atom   Prev_fragment_atom\") >>\n")
			for connection in self.__connection:
				obj_output.write("{0}\n".format(" ".join(["{0:>9}".format(x) for x in connection])))
			obj_output.write("\n")

			for line in self.__others:
				redb_NF = re_NF.search(line)
				if redb_NF:
					line = line[: redb_NF.start()] + "NF={0}".format(len(self.__fragments)) + line[redb_NF.end() :]
				obj_output.write(line + "\n")


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "editfrag.py", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-b", dest = "base_structure", metavar = "WHOLE.pdb", required = True, help = "whole structure for ABINIT-MP")
	parser.add_argument("-n", dest = "add_structures", metavar = "ADD.pdb", required = True, nargs = "+", help = "fragment structure")
	parser.add_argument("-f", dest = "fred", metavar = "BASE.fred", required = True, help = "original fred file")
	parser.add_argument("-o", dest = "output_file", metavar = "NEW.fred", required = True, help = "output fred")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	basic.check_exist(args.base_structure, 2)
	base_structure = MoleculeInformation(args.base_structure)

	basic.check_exist(args.fred, 2)
	obj_fred = FredData(args.fred)

	for new_structure_file in args.add_structures:
		basic.check_exist(new_structure_file, 2)
		new_structure = MoleculeInformation(new_structure_file)
		new_structure.update_index(base_structure)
		obj_fred.add_fragment(FragmentData(new_structure.output_fragmentdata()))
		sys.stderr.write("Added new fragment ({0}) to end of fragments\n".format(new_structure_file))

	obj_fred.write_file(args.output_file)
