#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
autofrag
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import os
import re
import argparse
from mods.func_prompt_io import *



# =============== class =============== #
class FragmentData:
	def __init__(self, flag_sep=False):
		self._atom_idxs = {}
		self._residue_name = ""
		self._charge = 0
		self._type = ""
		self._flag_sep = flag_sep
		self._shift_atoms = {}
		self._sep_atoms = {}
		self._connectivity = [-1, -1]
		self._sep_connectivity = [-1, -1]
		self._bda = 0
		self._sep_bda = 0


	def append_data(self, atom_order, atom_type, residue_name):
		""" データを追加する関数 """
		self._atom_idxs[atom_order] = atom_type.replace("+", "").replace("-", "")
		if self._atom_idxs[atom_order] in residue_types["Ion"] and ("+" in residue_name or "-" in residue_name):
			residue_name = residue_name.replace("+", "").replace("-", "")
		self._residue_name = residue_name


	def terminate(self):
		""" 電荷などを計算する関数 """
		if self._type != "TER":
			self._type = check_residue(self._residue_name)

			if self._type == "AminoAcid":
				if self._residue_name in ["LYS", "ARG", "HIP", "SYM", "ACE"]:
					self._charge = 1
				elif self._residue_name in ["ASP", "GLU", "NME"]:
					self._charge = -1
				elif self._residue_name in ["ACE", "ALA", "ASN", "CYS", "GLN", "GLY", "HIS", "HID", "HIE", "ILE", "LEU", "MET", "NME", "PHE", "PRO", "SER", "SYM", "THR", "TRP", "TYR", "VAL"]:
					self._charge = 0

			elif self._type == "NucleicAcid":
				if re_5term.search(self._residue_name):
					if self._flag_sep:
						# 塩基と分けている場合
						self._charge = 1
					else:
						self._charge = 0
				elif re_3term.search(self._residue_name):
					if self._flag_sep:
						# 塩基と分けている場合
						self._charge = 0
					else:
						self._charge = -1

			elif self._residue_name in ["Na", "K"]:
				self._charge = 1
			elif self._residue_name in ["Zn", "Ca", "CAL", "ZIN"]:
				self._charge = 2
			elif self._type == "Water":
				self._charge = 0
			else:
				self._charge = "ERR"
				self._bda = "ERR"

			# BDA のため、フラグメント間で原子を移動させる準備 (他フラグメントに与える原子)
			if self._type == "AminoAcid":
				# アミノ酸の切り方
				delete_key = []
				for key, value in self._atom_idxs.items():
					if value in ["C", "O"]:
						self._shift_atoms[key] = value
						delete_key.append(key)
						# del(self._atom_idxs[key])
						if value == "C":
							self._connectivity[1] = key
					if value == "CA":
						self._connectivity[0] = key
				self._atom_idxs = {k: v for k, v in self._atom_idxs.items() if k not in delete_key}

			elif self._type == "NucleicAcid":
				# 核酸の切り方
				ref_atom_idx = self._atom_idxs.copy()
				if self._flag_sep:
					# 塩基原子を sep_atoms に移動させる
					for key, value in ref_atom_idx.items():
						if not ("'" in value or "P" in value or re_termH1.match(value) or re_termH2.match(value)):
							# 接続情報追加
							if re_pyrimidine.match(self._residue_name) and value == "N1":
								self._sep_connectivity[1] = key
							elif re_purine.match(self._residue_name) and value == "N9":
								self._sep_connectivity[1] = key
							# 原子の移動
							self._sep_atoms[key] = value
							del(self._atom_idxs[key])
						elif value == "C1'":
							# 接続情報追加
							self._sep_connectivity[0] = key

				ref_atom_idx = self._atom_idxs.copy()
				for key, value in ref_atom_idx.items():
					# フラグメントで移動させる
					if "5" not in self._residue_name and value in ["P", "O1P", "O2P", "OP1", "OP2", "C5'", "O5'", "C5*", "O5*", "H5'", "H5''", "H5'1", "H5'2"]:
						# 接続情報追加
						if value == "C5'":
							self._connectivity[0] = key
						# 原子の移動
						self._shift_atoms[key] = value
						del(self._atom_idxs[key])
					elif "5" not in self._residue_name and value == "C4'":
						# 接続情報追加
						self._connectivity[1] = key

			if -1 not in self._connectivity:
				# 接続相手がある場合
				self._bda = 1

			if self._flag_sep and -1 not in self._sep_connectivity:
				self._sep_bda = 1




# =============== variables =============== #
residue_types = {
	"AminoAcid": ["ACE", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HIP", "HID", "HIE", "ILE", "LEU", "LYS", "MET", "NME", "PHE", "PRO", "SER", "SYM", "THR", "TRP", "TYR", "VAL"],
	"NucleicAcid": ["DA5", "DT5", "DG5", "DC5", "DA3", "DT3", "DG3", "DC3", "DA", "DT", "DG", "DC", "RA5", "RU5", "RG5", "RC5", "RA3", "RU3", "RG3", "RC3", "RA", "RU", "RG", "RC", "DNA_base"],
	"Water": ["SOL", "WAT", "HOH"],
	"Ion": ["Na", "Mg", "K", "Ca", "Cl", "Zn", "CAL", "ZIN"]
}
re_sign = re.compile(r"[-\+\*\'\"]")
re_head = re.compile(r"^([\d'\s])")
re_ter = re.compile(r"^TER")
re_termH1 = re.compile(r"H[53]T")
re_termH2 = re.compile(r"HO[53]")
re_pyrimidine = re.compile(r'^\s*[RD][UTC][53]?\s*$')
re_purine = re.compile(r'^\s*[RD][GA][53]?\s*$')
re_5term = re.compile(r'^[RD][AGCTU]5$')
re_3term = re.compile(r'^[RD][AGCTU]3$')
program_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
template_files = [
	os.path.join(program_dir, "template", "autofrag_3.templ"),
	os.path.join(program_dir, "template", "autofrag_5.templ"),
	os.path.join(program_dir, "template", "autofrag_m.templ")
]



# =============== functions =============== #
def check_residue(residue_name):
	""" 残基の種類を調べる関数 """
	for key, value in residue_types.items():
		if residue_name in value:
			return key
	return "Other"



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="autofrag", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest="INPUT", metavar="INPUT.pdb", required=True, help="Input file")
	parser.add_argument("-o", dest="OUTPUT", metavar="OUTPUT.fred", required=True, help="Output file")
	parser.add_argument("-S", "--separate", dest="FLAG_SEP", action="store_true", default=False, help="Nucleotide was devided into sugar+phosphate group and base")
	parser.add_argument("-V", "--version", dest="VERSION", metavar="version", required=True, help="ajf version (3: ABINIT-MP 3, 5: ABINIT-MP 5 or later, m: mizuho ABINIT-MP) or template file path")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.INPUT, 2)

	template_file = ""
	if args.VERSION == "3":
		template_file = template_files[0]
	elif args.VERSION == "5":
		template_file = template_files[1]
	elif args.VERSION == "m":
		template_file = template_files[2]
	else:
		template_file = args.VERSION
	check_exist(template_file, 2)

	option_nuc = "+base"
	if args.FLAG_SEP:
		option_nuc = "/base"

	title = args.INPUT.replace(".pdb", "")
	total_atom = 0
	fragments = []
	fragment = None
	residue_info_old = ""
	idx = -1

	# ファイル読み込み
	with open(args.INPUT, "r") as obj_input:
		connectivity = [[], []]
		for line in obj_input:
			line = line.rstrip("\n")

			if len(line) != 3:
				record = line[0:6]
			else:
				record = line

			if "ATOM" in record or "HETATM" in record:
				# 原子の場合

				# 全原子数カウント
				total_atom += 1

				# 原子順序番号取得
				atom_order = int(line[6:11].strip())

				# アトムタイプ取得
				atom_type = line[12:16].strip().replace("*", "'")
				redb = re_head.search(atom_type)
				if redb:
					atom_type = re_head.sub("", atom_type) + redb.group()

				# 残基情報取得
				residue_info = line[17:26]
				residue_name = re_sign.sub("", line[17:20].strip())

				if residue_info != residue_info_old:
					# 前回の残基と異なる場合、最終処理とフラグメントの新規作成
					if fragment is not None:
						fragments[idx].terminate()
					idx += 1
					fragments.append(FragmentData(args.FLAG_SEP))
					fragments[idx].append_data(atom_order, atom_type, residue_name)
					residue_info_old = residue_info
				else:
					# 前回の残基と同じ場合、登録
					fragments[idx].append_data(atom_order, atom_type, residue_name)

			elif "TER" in record:
				# TER の場合
				if fragments[-1].type != "TER":
					# TER の重複でない場合
					fragments.append(FragmentData(args.FLAG_SEP))
					idx += 1
					fragments[idx].type = "TER"

	# 電荷計算やフラグメント間の原子移動の準備
	for idx in range(len(fragments)):
		fragments[idx].terminate()

	# フラグメント間の原子移動
	new_fragments = []
	atoms_give = {}
	new_idx = -1
	total_charge = 0
	ligands = {}
	for idx in range(len(fragments)):
		# フラグメント原子の受け渡し
		if fragments[idx].type == "AminoAcid":
			if idx < len(fragments) - 1:
				if fragments[idx + 1].type != "TER":
					# 次のフラグメントが TER でない場合
					fragments[idx + 1].atom_idxs.update(fragments[idx].shift_atoms)
					fragments[idx].shift_atoms = {}
				else:
					# 次のフラグメントが TER の場合、元に戻す
					fragments[idx].atom_idxs.update(fragments[idx].shift_atoms)
			else:
				# 末端のフラグメントの場合、元に戻す
				fragments[idx].atom_idxs.update(fragments[idx].shift_atoms)
		elif fragments[idx].type == "NucleicAcid":
			if idx < len(fragments) - 1:
				if fragments[idx + 1].type != "TER":
					# 次のフラグメントが TER でない場合
					fragments[idx].atom_idxs.update(fragments[idx + 1].shift_atoms)
					fragments[idx + 1].shift_atoms = {}
				else:
					# 次のフラグメントが TER の場合、元に戻す
					fragments[idx].atom_idxs.update(fragments[idx].shift_atoms)
			else:
				# 末端のフラグメントの場合、元に戻す
				fragments[idx].atom_idxs.update(fragments[idx].shift_atoms)
		elif fragments[idx].type == "Ion":
			ligands[fragments[idx].residue_name] = fragments[idx].charge
		elif fragments[idx].type == "Other":
			ligands[fragments[idx].residue_name] = 0

		if fragments[idx].charge != "ERR":
			total_charge += fragments[idx].charge
		new_fragments.append(fragments[idx])
		new_idx += 1

		if fragments[idx].type == "NucleicAcid" and len(fragments[idx].sep_atoms) != 0:
			# 塩基と分割する場合、塩基のフラグメントを追加
			fragment = FragmentData(args.flag_sep)
			fragment.atom_idxs = new_fragments[new_idx].sep_atoms
			new_fragments[new_idx].sep_atoms = {}
			fragment.connectivity = new_fragments[new_idx].sep_connectivity
			new_fragments[new_idx].sep_connectivity = [-1, -1]
			fragment.charge = -1
			fragment.type = "NucleicAcid"
			fragment.residue_name = new_fragments[new_idx].residue_name
			fragment.bda = new_fragments[new_idx].sep_bda
			new_fragments[new_idx].sep_bda = 0
			new_fragments.append(fragment)
			total_charge += fragment.charge
			new_idx += 1

	fragments = new_fragments
	del(new_fragments)

	ligand_list = []
	for residue, charge in ligands.items():
		ligand_list.append("{0}={1}".format(residue, charge))
	ligand = ",".join(ligand_list)

	# 出力
	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT)

	cnt_fragment = 0
	with open(args.OUTPUT, "w") as obj_output:
		obj_output.write("  FNo.  | Charge | BDA | Atoms of fragment\n")
		for idx, obj in enumerate(fragments):
			if obj.type != "TER":
				cnt_fragment += 1
				atom_order = sorted(list(map(lambda x : int(x), obj.atom_idxs.keys())))
				atom_order = map(lambda x : "{0:>8}".format(x), atom_order)
				atom_list = " ".join(list(atom_order))
				if obj.charge == "ERR":
					sys.stderr.write(" WARNING: Fragment No. {0} cannot determine charge and BDA, because unknown residue\n".format(cnt_fragment))
				obj_output.write("{fragment_number:>7} |  {charge:>3}   | {bda:^3} | {atom}\n".format(fragment_number = cnt_fragment, charge = obj.charge, bda = obj.bda, atom = atom_list))
		obj_output.write("\n")
		obj_output.write("<< connections (ex. \"Next_fragment_atom   Prev_fragment_atom\") >>\n")
		for obj in fragments:
			if obj.type != "TER":
				if -1 not in obj.connectivity:
					connect_list = " ".join(list(map(lambda x : "{0:>9}".format(x), obj.connectivity)))
					obj_output.write("{0}\n".format(connect_list))
		obj_output.write("\n")
		with open(template_file, "r") as obj_input:
			for line in obj_input:
				if "$title" in line:
					line = line.replace("$title", title)
				if "$total_atom" in line:
					line = line.replace("$total_atom", str(total_atom))
				if "$total_charge" in line:
					line = line.replace("$total_charge", str(total_charge))
				if "$in" in line:
					line = line.replace("$in", args.INPUT)
				if "$cpf" in line:
					line = line.replace("$cpf", title + ".cpf")
				if "$option_nuc" in line:
					line = line.replace("$option_nuc", option_nuc)
				if "$total_fragment" in line:
					line = line.replace("$total_fragment", str(cnt_fragment))
				if "$ligand" in line:
					line = line.replace("$ligand", ligand)
				obj_output.write(line)
