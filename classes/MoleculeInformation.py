#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

# =============== classes =============== #
class MoleculeInformation:
	""" 分子構造のクラス """
	def __init__(self, input_file):
		self._residue_names = []
		self._residue_idxs = []
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
					self._residue_names.append(line[17:20])
					self._residue_idxs.append(int(line[22:26]))
					self._atom_names.append(line[12:16])
					self._atom_idxs.append(int(line[6:11]))
					self._coords.append(line[30:54])


	def get_info(self, data_type = None, idx = None):
		""" 原子情報を返すメソッド """
		if data_type == "residue_name":
			# 残基名を返す
			if idx is not None:
				try:
					return self._residue_names[idx]
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return self._residue_names

		elif data_type == "residue_idx":
			# 残基インデックスを返す
			if idx is not None:
				try:
					return self._residue_idxs[idx]
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return self._residue_idxs

		elif data_type == "atom_name":
			if idx is not None:
				try:
					return self._atom_names[idx]
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return self._atom_names

		elif data_type == "atom_idx":
			if idx is not None:
				try:
					return self._atom_idxs[idx]
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return self._atom_idxs
		else:
			if idx is not None:
				try:
					return [[res_name, res_idx, atom_name, atom_idx] for res_name, res_idx, atom_name, atom_idx in zip(self._residue_names, self._residue_idxs, self._atom_names, self._atom_idxs)][idx]
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [[res_name, res_idx, atom_name, atom_idx] for res_name, res_idx, atom_name, atom_idx in zip(self._residue_names, self._residue_idxs, self._atom_names, self._atom_idxs)]


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

	def output_fragmentdata(self, obj_base_molecule, flag_multi):
		""" 該当するフラグメントを出力するメソッド """
		fragment_records = [[]]

		if flag_multi:
			# 同じ種類の残基の複数のフラグメントに適用する場合
			mol_info = [list(x) for x in list(zip(self.get_info("residue_name"), self.get_info("atom_name")))]
			len_atom = len(self._atom_names)
			check_list = [False] * len_atom

			for info in obj_base_molecule.get_info():
				# 全体構造の情報で回す
				if [info[0], info[2]] in mol_info:
					# 残基名と原子名のリストの情報が含まれる場合
					match_idx = mol_info.index([info[0], info[2]])

					if check_list[match_idx]:
						# 前回と同じ原子を処理している場合、新たなフラグメントリストを作成する
						fragment_records.append([])
						check_list = [False] * len_atom

					fragment_records[-1].append(info[3])
					check_list[match_idx] = True

		else:
			# 単一のフラグメントに適用する場合、座標でマッチさせる
			coords_base = obj_base_molecule.get_coord()

			for array_idx, coords in enumerate(self._coords):
				# 自身の座標とそのインデックスで回す
				if coords in coords_base:
					# 座標でマッチした場合、マッチしたインデックスを取得する
					match_idx = coords_base.index(coords)
					if self._atom_names[array_idx] == obj_base_molecule.get_info("atom_name", match_idx):
						# 原子名でもマッチを確認し、フラグメント情報に追加する
						fragment_records[-1].append(obj_base_molecule.get_info("atom_idx", match_idx))
					else:
						sys.stderr.write("ERROR: wrong coordinates in MoleculeInformation.output_fragmentdata().\n")
						sys.exit(1)

		fragment_records = ["{0}|{1}|{2}|{3}".format(0, "ERR", "ERR", " ".join([str(x) for x in record])) for record in fragment_records]

		return fragment_records


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
