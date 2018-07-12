#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import parmed

# =============== classes =============== #
class MoleculeInformation:
	""" 分子構造のクラス """
	def __init__(self, input_file):
		# member
		self._obj_mol = None

		# initiation
		self._load_file(input_file)


	def _load_file(self, input_file):
		""" ファイルを読み込むメソッド """
		self._obj_mol = parmed.load_file(input_file)


	def get_info(self, data_type = None, idx = None):
		""" 原子情報を返すメソッド """
		if data_type == "residue_name":
			# 残基名を返す
			if idx is not None:
				try:
					return self._obj_mol.atoms[idx].residue.name
				except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
			else:
				return [obj_atom.residue.name for obj_atom in self._obj_mol.atoms]

		elif data_type == "residue_idx":
			# 残基インデックスを返す
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
		""" マスクを原子順序番号に変換するメソッド """
		mask = parmed.amber.AmberMask(self._obj_mol, mask)
		return [self.get_info("atom_idx", idx) for idx in mask.Selected()]


	def get_coord(self, idx = None):
		""" 座標を出力するメソッド """
		if idx is not None:
			try:
				return [[obj_atom.xx, obj_atom.xy, obj_atom.xz] for obj_atom in self._obj_mol.atoms[idx]]
			except IndexError:
					sys.stderr.write("ERROR: Invalid index in MoleculeInformation.get_info.\n")
					sys.exit(1)
		else:
			return [[obj_atom.xx, obj_atom.xy, obj_atom.xz] for obj_atom in self._obj_mol.atoms]


	def output_fragmentdata(self, obj_base_molecule, flag_multi):
		""" 該当するフラグメントを出力するメソッド """
		fragment_records = [[]]

		if flag_multi:
			# 同じ種類の残基の複数のフラグメントに適用する場合
			mol_info = [list(x) for x in list(zip(self.get_info("residue_name"), self.get_info("atom_name")))]
			len_atom = len(self.get_info("atom_name"))
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

			for array_idx, coords in enumerate(self.get_coords()):
				# 自身の座標とそのインデックスで回す
				if coords in coords_base:
					# 座標でマッチした場合、マッチしたインデックスを取得する
					match_idx = coords_base.index(coords)
					if self.get_info("atom_name", array_idx) == obj_base_molecule.get_info("atom_name", match_idx):
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
