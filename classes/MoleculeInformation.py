#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

# =============== classes =============== #
class MoleculeInformation:
	""" 分子構造のクラス """
	def __init__(self, input_file):
		self._atom_names = []
		self._atom_idxs = []
		self._coords = []
		self._load_file(input_file)

	def _load_file(self, input_file):
		""" ファイルを読み込むメソッド """
		re_atom = re.compile(r"^(?:(?:ATOM)|(?:HETATM))")
		with open(input_file, "r") as obj_input:
			atom_idx = 0
			for line in obj_input:
				if re_atom.search(line):
					line = line.rstrip("\r\n")
					atom_idx += 1
					self._atom_names.append(line[12:16])
					self._atom_idxs.append(int(line[6:11]))
					self._coords.append(line[30:54])

	def get_atom_index(self, idx = None):
		""" 原子順序番号を出力するメソッド """
		if idx is None:
			return self._atom_idxs
		else:
			try:
				return self._atom_idxs[idx]
			except IndexError:
				sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_atom_index\n")
				sys.exit(1)

	def get_atom_name(self, idx = None):
		""" 原子名を出力するメソッド """
		if idx is None:
			return self._atom_names
		else:
			try:
				return self._atom_names[idx]
			except IndexError:
				sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_atom_name\n")
				sys.exit(1)

	def get_coord(self, idx = None):
		""" 座標を出力するメソッド """
		if idx is None:
			return self._coords
		else:
			try:
				return self._coords[idx]
			except IndexError:
				sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_coord\n")
				sys.exit(1)

	def update_index(self, obj_base_molecule):
		""" 原子順序番号を更新するメソッド """
		coords_base = obj_base_molecule.get_coord()

		for array_idx, coords in enumerate(self._coords):
			if coords in coords_base:
				match_idx = coords_base.index(coords)
				if self._atom_names[array_idx] == obj_base_molecule.get_atom_name(match_idx):
					self._atom_idxs[array_idx] = obj_base_molecule.get_atom_index(match_idx)
				else:
					sys.stderr.write("ERROR: wrong coordinates in MoleculeInformation.update_index()\n")
					sys.exit(1)

	def output_fragmentdata(self):
		return "{0}|{1}|{2}|{3}".format(0, "ERR", "ERR", " ".join([str(x) for x in self._atom_idxs]))


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
