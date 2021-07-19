#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fred4 - fragment editor for mizuho ABINIT-MP
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import os
import re
from mods.func_prompt_io import *



# =============== common variables =============== #
# general
RE_WSP = re.compile(r"[\s\t]+")
RE_QUOTE_H = re.compile(r"^['\"]")
RE_QUOTE_T = re.compile(r"['\"]$")
RE_EMPTY = re.compile(r"^[\s\t]*$")
RE_NAMELIST_H = re.compile(r"^[\s\t]*\&")
RE_NAMELIST_T = re.compile(r"^[\s\t]*/$")
RE_DIGIT = re.compile(r"[\d\s]+")


# ネームリスト更新
RE_NATOM = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
RE_NF = re.compile(r"NF", re.IGNORECASE)
RE_CHARGE = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
RE_AUTOFRAG = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
RE_FRAGMENT = re.compile(r"\&FRAGMENT", re.IGNORECASE)
RE_COORD = re.compile(r"\{\.{3}\}")

ATOM_CHARGES = {
	"H": 1, "Li": 3, "Be": 4, "B": 5,
	"C": 6, "N": 7, "O": 8, "F": 9,
	"Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
	"S": 16, "Cl": 17, "K": 19, "Ca": 20,
	"Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
	"Br": 35
}


# =============== functions =============== #
# split_n
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


# check_charge
def check_charge(fragment_members, charges, pdb):
	atom_orders = []
	atom_types = []

	with open(pdb, "r") as obj_pdb:
		for line_val in obj_pdb:
			if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
				atom_orders.append(line_val[6:11].strip())
				atom = RE_DIGIT.sub("", line_val[12:14].strip())
				atom = RE_QUOTE_H.sub("", atom)
				atom = RE_QUOTE_T.sub("", atom)
				if atom in ["HO", "HH"]:
					atom = "H"
				atom_types.append(atom)
				if not atom in ATOM_CHARGES:
					sys.stderr.write("ERROR: Unknown atomtype (%s). Skipped...\n" % atom)

	flag_error = 0
	for i in range(len(fragment_members)):
		charge = 0
		electron_info = []
		for j in range(len(fragment_members[i])):
			try:
				charge += ATOM_CHARGES[atom_types[atom_orders.index(str(fragment_members[i][j]))]]
				electron_info.append("{0:>5} {1} = {2}\n".format(fragment_members[i][j], atom_types[atom_orders.index(str(fragment_members[i][j]))], ATOM_CHARGES[atom_types[atom_orders.index(str(fragment_members[i][j]))]]))
			except ValueError:
				sys.stderr.write("ERROR: %d is not in list. Check the atom order in fred and pdb file.\n" % fragment_members[i][j])
				sys.exit(1)
		charge += charges[i]
		if charge % 2 != 0:
			sys.stderr.write("ERROR: Invalid number of fragment electrons.\n       The number of electrons in fragment No. %d is %d.\n" % (i + 1, charge))
			for error_line in electron_info:
				sys.stderr.write("       " + error_line)
			sys.stderr.write("\n")
			flag_error = 1

	if flag_error == 1:
		sys.stderr.write("ERROR: The number of electrons for some fragments were not even number.\n       Prceeding? (y/N): ")
		sys.stderr.flush()
		user = sys.stdin.readline().rstrip("\r\n")
		if user.lower() != "y":
			sys.exit(0)
	else:
		sys.stderr.write("INFO: check_charge is ok.\n")


# convert_ajf
def convert_ajf(lists, width, n):
	lists = list(map(lambda data : str(data), lists))
	new_lists = [""]

	format_data = "%" + str(width) + "s"

	index = 0
	for i in range(len(lists)):
		new_lists[index] += format_data % lists[i]
		if i != 0 and (i + 1) % n == 0:
			index += 1
			new_lists.append("")

	if len(new_lists[len(new_lists) - 1]) == 0:
		# 指定された個数でちょうど終わる場合は、無駄な要素が追加されるので削除する
		del(new_lists[len(new_lists) - 1])

	return new_lists


# load_ajf
def load_ajf(file_input, file_reference):
	re_namelist_fragment_h = re.compile(r"^[\s\t]*\&FRAGMENT")
	re_namelist_fragment_t = re.compile(r"^[\s\t]*\/")

	flag_read = 0
	atom = 0
	atom_now = 0
	fragment = 0
	fragment_now = 0
	namelists = []
	fragment_atoms = []
	fragment_members = []
	fragment_members_tmp = []
	charges = []
	BDAs = []
	connections = []

	with open(file_input, "r") as obj_input:
		for line_idx, line_val in enumerate(obj_input, 1):
			if flag_read == 0:
				# ネームリスト取得
				line_val = line_val.strip()

				if len(line_val.strip()) == 0:
					# 空行はスキップ
					continue

				if re_namelist_fragment_h.search(line_val):
					# FRAGMENT ネームリストの始端
					flag_read = 1

				elif "ReadGeom" in line_val:
					if file_reference != None:
						# 参照 PDB が指定されていた場合
						check_file(file_reference)
					else:
						# 参照 PDB が指定されていない場合
						file_reference = RE_WSP.sub("", line_val)
						file_reference = file_reference.replace("ReadGeom=", "")
						file_reference = RE_QUOTE_H.sub("", file_reference)
						file_reference = RE_QUOTE_T.sub("", file_reference)

					check_file(file_reference)

					with open(file_reference, "r") as obj_pdb:
						for p_line in obj_pdb:
							if p_line.startswith("ATOM") or pline.startswith("HETATM"):
								atom += 1

				if re_namelist_fragment_t.search(line_val):
					line_val = "/\n"
				elif not RE_NAMELIST_H.search(line_val):
					line_val = "  " + line_val
				namelists.append(line_val)

			elif 0 < flag_read:
				# フラグメント情報取得
				line_val = line_val.rstrip("\r\n")

				if re_namelist_fragment_t.search(line_val):
					# FRAGMENT ネームリストの終端
					line_val = "{...}\n/\n"
					flag_read = 0
					namelists.append(line_val)

				else:
					# 8 文字ずつ区切る
					datas = list(map(lambda data : int(data), split_n(line_val, 8)))

					if flag_read == 1:
						# フラグメント構成原子取得
						if "0" in datas:
							# 構成原子 0 のフラグメントがある場合
							sys.stderr.write("ERROR: Invalid ajf file in {0}\n".format(line_idx))
							sys.stderr.write("       zero atoms in fragment found\n")
							sys.exit(1)

						fragment += len(datas)
						fragment_atoms.extend(datas)

						atom_now += sum(datas)
						if atom == atom_now:
							# PDB の総原子数と現在の原子数が一致した場合
							flag_read = 2
						elif atom < atom_now:
							# PDB の総原子数以上の原子を検出した場合
							sys.stderr.write("ERROR: The number of atoms in PDB file (%s) and ajf file mismatched (%d / %d)\n" % (file_reference, atom_now, atom))
							sys.stderr.write("       Maybe wrong PDB file was specified.\n")
							sys.stderr.write("       Please fix ReadGeom in ajf file, OR use -P option.\n")
							sys.exit(1)

					elif flag_read == 2:
						# 電荷情報取得
						fragment_now += len(datas)
						charges.extend(datas)

						if fragment == fragment_now:
							# フラグメント数が総フラグメント数と一致した場合
							flag_read = 3
							fragment_now = 0
						elif fragment < fragment_now:
							# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
							sys.stderr.write("ERROR: The number of fragments mismatched in charge section (%d / %d)\n" % (fragment_now, fragment))
							sys.exit(1)

					elif flag_read == 3:
						# BDA 取得
						fragment_now += len(datas)
						BDAs.extend(datas)

						if fragment == fragment_now:
							# フラグメント数が総フラグメント数と一致した場合
							flag_read = 4
							fragment_now = 0
						elif fragment < fragment_now:
							# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
							sys.stderr.write("ERROR: The number of fragments mismatched in BDA section\n" % (fragment_now, fragment))
							sys.exit(1)

					elif flag_read == 4:
						# フラグメント構成原子の原子順序番号取得
						fragment_members_tmp.extend(datas)

						if fragment_atoms[fragment_now] <= len(fragment_members_tmp):
							fragment_now += 1
							fragment_members.append(fragment_members_tmp)
							fragment_members_tmp = []

							if fragment == fragment_now:
								# フラグメント数が総フラグメント数と一致した場合
								flag_read = 5
								now_fragment = 0
							elif fragment < fragment_now:
								# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
								sys.stderr.write("ERROR: The number of fragments mismatched in fragment atom section\n" % (fragment_now, fragment))
								sys.exit(1)

					elif flag_read == 5:
						# 接続情報取得
						line_val = line_val.strip()
						connections.append(list(map(lambda data : int(data),RE_WSP.split(line_val))))
	return fragment_atoms, charges, BDAs, fragment_members, connections, namelists


# load_fred
def load_fred(file_input):
	flag_read = 0
	fragment = 0
	atom = 0
	flag_fragment = 0

	fragment_atoms = []
	charges = []
	BDAs = []
	fragment_members = []
	connections = []
	namelists = []

	re_connection = re.compile(r"<< connections", re.IGNORECASE)
	re_namelist_mark = re.compile(r"=+< namelist >=+")
	RE_NATOM = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
	RE_NF = re.compile(r"NF", re.IGNORECASE)
	RE_CHARGE = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
	RE_AUTOFRAG = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
	RE_FRAGMENT = re.compile(r"\&FRAGMENT", re.IGNORECASE)
	RE_COORD = re.compile(r"\{\.{3}\}")

	with open(file_input, "r") as obj_input:
		for line_idx, line_val in enumerate(obj_input, 1):
			line_val = line_val.strip()

			if line_idx == 1 or RE_EMPTY.search(line_val):
				# 1行目と空行は無視
				continue

			elif re_connection.search(line_val):
				# 接続情報フラグ
				flag_read = 1

			elif re_namelist_mark.search(line_val):
				# この行以降がネームリスト
				flag_read = 2

				# フラグメント情報の統計
				# フラグメントに含まれる原子数算出
				for i in range(len(fragment_members)):
					fragment_atoms.append(len(fragment_members[i]))

				# フラグメント数算出
				fragment = len(fragment_members)

				# 電荷の合計算出
				charge = sum(charges)

			elif flag_read == 0:
				# フラグメント情報
				datas = line_val.split("|")
				datas = list(map(lambda data : data.strip(), datas))
				for item in datas:
					if item == "":
						sys.stderr.write("ERROR: Invalid format in line {0}.\n".format(line_idx))
						sys.exit(1)

				charges.append(int(datas[1]))
				BDAs.append(int(datas[2]))
				tmp = list(map(lambda data : int(data), RE_WSP.split(datas[3])))
				fragment_members.append(tmp)
				atom += len(tmp)

			elif flag_read == 1:
				# 接続情報
				tmp = list(map(lambda data : int(data), RE_WSP.split(line_val)))
				connections.append(tmp)

			elif flag_read == 2:
				# ネームリスト
				if RE_NATOM.search(line_val):
					# 原子数の更新
					line_val = re.sub(r"\d+", str(atom), line_val)
				elif RE_NF.search(line_val):
					# フラグメント数の更新
					line = re.sub(r"\d+", str(fragment), line_val)
				elif RE_CHARGE.search(line):
					# 電荷の更新
					line_val = re.sub(r"-?\d+", str(charge), line_val)
				elif RE_AUTOFRAG.search(line_val):
					# autofrag の更新
					line_val = re.sub(r"=.+$", "='OFF'", line_val)

				if RE_NAMELIST_T.search(line_val):
					line = "/\n"
				elif not RE_NAMELIST_H.search(line_val):
					line_val = "  " + line_val

				namelists.append(line_val)
	return fragment_atoms, charges, BDAs, fragment_members, connections, namelists



# write_data (データの書き出し; ファイルが指定されていた場合はファイルに書き出し)
def write_data(contents, output):
	with open(output, "w") as obj_output:
		for line_val in contents:
			obj_output.write("{0}\n".format(line_val))



# =============== main =============== #
if __name__ == '__main__':
	try:
		parser = argparse.ArgumentParser(description="Fragment editor for mizuho ABINIT-MP", formatter_class=argparse.RawTextHelpFormatter)

		subparser = parser.add_subparsers(help="Sub-command")
		subparser.required = True

		parser_edit = subparser.add_parser("edit", help="Convert ajf to fred (ajf -> fred)")
		parser_edit.set_defaults(func="edit")
		parser_edit.add_argument("-i", dest="INPUT_FILE", metavar="INPUT", required=True, help="ajf file")
		parser_edit.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT", help="output file")
		parser_edit.add_argument("-p", "--pdb", metavar="PDB", help="reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
		parser_edit.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_rewrite = subparser.add_parser("rewrite", help="Rewrite fred fred -> fred")
		parser_rewrite.set_defaults(func="rewrite")
		parser_rewrite.add_argument("-i", dest="INPUT_FILE", metavar="INPUT", required=True, help="fred")
		parser_rewrite.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT", help="Output (Default: STDOUT)")
		parser_rewrite.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_output = subparser.add_parser("output", help="Convert fred to ajf (fred -> ajf)")
		parser_output.set_defaults(func="output")
		parser_output.add_argument("-i", dest="INPUT_FILE", metavar="INPUT", required=True, help="fred")
		parser_output.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT", help="Output (Default: STDOUT)")
		parser_output.add_argument("-p", "--pdb", metavar="PDB", help="Reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
		parser_output.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_autofrag = subparser.add_parser("autofrag", help="Auto fragmentation for PDB (pdb -> fred)")
		parser_autofrag.set_defaults(func="autofrag")
		parser_autofrag.add_argument("-i", dest="INPUT_FILE", metavar="INPUT", required=True, help="PDB")
		parser_autofrag.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT", help="Output (Default: STDOUT)")
		parser_autofrag.add_argument("-s", "--separate", action="store_true", help="Nucleotide is separates to base and sugar+phosphate")
		parser_autofrag.add_argument("-v", "--version", choices=["3", "5", "m"], help="ajf version: 3 = abinit-mp3, 5 = abinitmp5, m = mizuho")
		parser_autofrag.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_editfrag = subparser.add_parser("editfrag", help="Create new fred in which fragments were devided based on PDB and fred (pdb + fred -> fred)")
		parser_editfrag.set_defaults(func="editfrag")
		parser_editfrag.add_argument("-i", dest="INPUT_FILE", metavar="INPUT", required=True, help="PDB")
		parser_editfrag.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT", help="Output (Default: STDOUT)")
		parser_editfrag.add_argument("-f", metavar="fred", dest="fred", required=True, help="fred")
		parser_editfrag.add_argument("-b", metavar="pdb", dest="pdb", required=True, help="pdb")
		parser_editfrag.add_argument("-n", metavar="pdb", dest="pdb", required=True, nargs="+", help="pdb for each fragments")
		parser_editfrag.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		args = parser.parse_args()
	except TypeError:
		sys.stderr.write("ERROR: No sub-command (autofrag | edit | rewrite | output | editfrag)\n")
		sys.exit(1)

	check_exist(args.INPUT_FILE, 2)

	if args.func == "edit":
		# 編集ファイルに変換

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = load_ajf(args.input, args.pdb)

		# 整形
		flag_fragment = 0
		output = []
		output.append("  FNo.  | charge | BDA | atoms of fragment")
		for i in range(len(fragment_atoms)):
			output.append("{0:7} |{1:6}  |{2:3}  |{3}".format(i + 1, charges[i], BDAs[i], " ".join(map(lambda data : "{0:8d}".format(data), fragment_members[i]))))
		output.append("\n<< connections (ex. \"Next_fragment_atom    Prev_fragment_atom\")>>")
		for i in range(len(connections)):
			output.append("".join(map(lambda data : "{0:8d}".format(data), connections[i])))
		output.append("\n")
		output.append("===============< namelist >===============")
		for line_val in namelists:
			if RE_NATOM.search(line_val):
				# 原子数の更新
				atom = sum(fragment_atoms)
				line_val = re.sub(r"\d+", str(atom), line_val)
			elif RE_NF.search(line_val):
				# フラグメント数の更新
				fragment = len(fragment_members)
				line_val = re.sub(r"\d+", str(fragment), line_val)
			elif RE_CHARGE.search(line_val):
				# 電荷の更新
				charge = sum(charges)
				line_val = re.sub(r"-?\d+", str(charge), line_val)
			elif RE_AUTOFRAG.search(line_val):
				# autofrag の更新
				line_val = re.sub(r"=.+$", "='OFF'", line_val)

			output.append(line_val)

		# 書き出し
		write_data(output, args.OUTPUT_FILE)


	elif args.func == "rewrite":
		# リロード

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = load_fred(args.INPUT_FILE)

		# 整形
		output = []
		output.append("  FNo.  | charge | BDA | atoms of fragment")
		for i in range(len(fragment_atoms)):
			output.append("%7s |%6s  |%3s  |%s" % (i + 1, charges[i], BDAs[i], " ".join(map(lambda data : "%8d" % data, fragment_members[i]))))
		output.append("\n<< connections (ex. \"Next_fragment_atom    Prev_fragment_atom\")>>")
		for i in range(len(connections)):
			output.append("".join(map(lambda data : "%8d" % data, connections[i])))
		output.append("\n")
		output.append("===============< namelist >===============")
		for line_val in namelists:
			if RE_NATOM.search(line_val):
				# 原子数の更新
				atom = sum(fragment_atoms)
				line_val = re.sub(r"\d+", str(atom), line_val)
			elif RE_NF.search(line_val):
				# フラグメント数の更新
				fragment = len(fragment_members)
				line_val = re.sub(r"\d+", str(fragment), line_val)
			elif RE_CHARGE.search(line_val):
				# 電荷の更新
				charge = sum(charges)
				line_val = re.sub(r"-?\d+", str(charge), line_val)
			elif RE_AUTOFRAG.search(line):
				# autofrag の更新
				line_val = re.sub(r"=.+$", "='OFF'", line_val)
			output.append(line_val)

		# 書き出し
		write_data(output, args.OUTPUT_FILE)

	elif args.func == "output":
		# ajf ファイルに変換

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = load_fred(args.INPUT_FILE)

		file_reference = ""
		if args.pdb != None:
			check_file(args.pdb)
			file_reference = args.pdb
		else:
			for item in namelists:
				if "ReadGeom" in item:
					file_reference = item.strip()
					file_reference = RE_WSP.sub("", file_reference)
					file_reference = file_reference.replace("ReadGeom=", "")
					file_reference = RE_QUOTE_H.sub("", file_reference)
					file_reference = RE_QUOTE_T.sub("", file_reference)
					check_exist(file_reference, 2)
					break

		check_charge(fragment_members, charges, file_reference)

		# 整形
		output = []
		flag_fragment = 0
		for line_val in namelists:
			if RE_NATOM.search(line_val):
				# 原子数の更新
				atom = sum(fragment_atoms)
				line_val = re.sub(r"\d+", str(atom), line_val)
			elif RE_NF.search(line_val):
				# フラグメント数の更新
				fragment = len(fragment_members)
				line_val = re.sub(r"\d+", str(fragment), line_val)
			elif RE_CHARGE.search(line_val):
				# 電荷の更新
				charge = sum(charges)
				line_val = re.sub(r"-?\d+", str(charge), line_val)
			elif RE_AUTOFRAG.search(line_val):
				# autofrag の更新
				line_val = re.sub(r"=.+$", "='OFF'", line_val)
			elif RE_FRAGMENT.search(line_val):
				# フラグメント情報書き出し場所を検索
				flag_fragment = 1
			elif flag_fragment == 1 and RE_COORD.search(line_val):
				# フラグメント情報の書き出し
				for line_val in convert_ajf(fragment_atoms, 8, 10):
					output.append(line_val)

				for line_val in convert_ajf(charges, 8, 10):
					output.append(line_val)

				for line_val in convert_ajf(BDAs, 8, 10):
					output.append(line_val)

				for j in range(len(fragment_members)):
					for line_val in convert_ajf(fragment_members[j], 8, 10):
						output.append(line_val)

				for j in range(len(connections)):
					for line_val in convert_ajf(connections[j], 8, 10):
						output.append(line_val)
				flag_fragment = 0
				continue
			output.append(line_val)

		# 出力
		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		write_data(output, args.OUTPUT_FILE)

	elif args.func == "autofrag":
		pass


	elif args.func == "editfrag":
		pass
