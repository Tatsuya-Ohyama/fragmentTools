#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fred4 - fragment editor for mizuho ABINIT-MP
"""

import sys, os, re
import argparse

# =============== common variables =============== #
# general
re_wsp = re.compile(r"[\s\t]+")
re_quote_h = re.compile(r"^['\"]")
re_quote_t = re.compile(r"['\"]$")
re_empty = re.compile(r"^[\s\t]*$")
re_namelist_h = re.compile(r"^[\s\t]*\&")
re_namelist_t = re.compile(r"^[\s\t]*/$")
re_pdb_atom = re.compile(r"^((HETATM)|(ATOM))")

# ネームリスト更新
re_natom = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
re_nf = re.compile(r"NF", re.IGNORECASE)
re_charge = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
re_autofrag = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
re_fragment = re.compile(r"\&FRAGMENT", re.IGNORECASE)
re_coord = re.compile(r"\{\.{3}\}")


# =============== functions =============== #
# check_file
def check_file(file):
	if not os.path.exists(file):
		sys.stderr.write("ERROR: No such file (%s)\n" % file)
		sys.exit(1)
	else:
		return file


# check_overwrite
def check_overwrite(file):
	if os.path.exists(file):
		sys.stderr.write("WARN: %s exists. Overwrite it? (y/N): " % file)
		sys.stderr.flush()
		user = sys.stdin.readline().strip("\r\n")
		if not (user == "Y" or user == "y"):
			sys.exit(0)
		else:
			return file


# split_n
def split_n(line, length):
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
	atom_charges = {"H": 1, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "K": 19, "Ca": 20, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Br": 35}
	atom_orders = []
	atom_types = []
	re_digit = re.compile(r"[\d\s]+")

	with open(pdb, "r") as obj_pdb:
		for line in obj_pdb:
			if re_pdb_atom.search(line):
				atom_orders.append(line[6:11].strip())
				atom = re_digit.sub("", line[12:14].strip())
				if atom == "HO":
					atom = "H"
				atom_types.append(atom)
				if not atom in atom_charges:
					sys.stderr.write("ERROR: Unknown atomtype (%s)\n" % atom)
					sys.exit(1)
				# atom_chages.append(atom_charges[atom])

	flag_error = 0
	for i in range(len(fragment_members)):
		charge = 0
		for j in range(len(fragment_members[i])):
			try:
				charge += atom_charges[atom_types[atom_orders.index(str(fragment_members[i][j]))]]
			except ValueError:
				sys.stderr.write("ERROR: %d is not in list. Check the atom order in fred and pdb file.\n" % fragment_members[i][j])
				sys.exit(1)
		charge += charges[i]
		if charge % 2 != 0:
			sys.stderr.write("ERROR: Invalid number of fragment electrons.\n       charge of fragment No. %d is %d.\n\n" % (i + 1, charge))
			flag_error = 1

	if flag_error == 1:
		sys.stderr.write("ERROR: The number of electrons for some fragments were not even number.\n       Prceeding? (y/N): ")
		sys.stderr.flush()
		user = sys.stdin.readline().rstrip("\r\n")
		if not (user == "Y" or user == "y"):
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
		line_count = 0
		for line in obj_input:
			line_count += 1

			if flag_read == 0:
				# ネームリスト取得
				line = line.strip()

				if re_empty.search(line):
					# 空行はスキップ
					continue

				if re_namelist_fragment_h.search(line):
					# FRAGMENT ネームリストの始端
					flag_read = 1

				elif "ReadGeom" in line:
					if file_reference != None:
						# 参照 PDB が指定されていた場合
						check_file(file_reference)
					else:
						# 参照 PDB が指定されていない場合
						file_reference = re_wsp.sub("", line)
						file_reference = file_reference.replace("ReadGeom=", "")
						file_reference = re_quote_h.sub("", file_reference)
						file_reference = re_quote_t.sub("", file_reference)

					check_file(file_reference)

					with open(file_reference, "r") as obj_pdb:
						for p_line in obj_pdb:
							if re_pdb_atom.search(p_line):
								atom += 1

				if re_namelist_fragment_t.search(line):
					line = "/\n"
				elif not re_namelist_h.search(line):
					line = "  " + line
				namelists.append(line)

			elif 0 < flag_read:
				# フラグメント情報取得
				line = line.rstrip("\r\n")

				if re_namelist_fragment_t.search(line):
					# FRAGMENT ネームリストの終端
					line = "{...}\n/\n"
					flag_read = 0
					namelists.append(line)

				else:
					# 8 文字ずつ区切る
					datas = list(map(lambda data : int(data), split_n(line, 8)))

					if flag_read == 1:
						# フラグメント構成原子取得
						if "0" in datas:
							# 構成原子 0 のフラグメントがある場合
							sys.stderr.write("ERROR: Invalid ajf file in %d\n" % line_count)
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
						line = line.strip()
						connections.append(list(map(lambda data : int(data),re_wsp.split(line))))
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
	re_natom = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
	re_nf = re.compile(r"NF", re.IGNORECASE)
	re_charge = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
	re_autofrag = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
	re_fragment = re.compile(r"\&FRAGMENT", re.IGNORECASE)
	re_coord = re.compile(r"\{\.{3}\}")

	with open(file_input, "r") as obj_input:
		line_count = 0

		for line in obj_input:
			line_count += 1
			line = line.strip()

			if line_count == 1 or re_empty.search(line):
				# 1行目と空行は無視
				continue

			elif re_connection.search(line):
				# 接続情報フラグ
				flag_read = 1

			elif re_namelist_mark.search(line):
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
				datas = line.split("|")
				datas = list(map(lambda data : data.strip(), datas))
				charges.append(int(datas[1]))
				BDAs.append(int(datas[2]))
				tmp = list(map(lambda data : int(data), re_wsp.split(datas[3])))
				fragment_members.append(tmp)
				atom += len(tmp)

			elif flag_read == 1:
				# 接続情報
				tmp = list(map(lambda data : int(data), re_wsp.split(line)))
				connections.append(tmp)

			elif flag_read == 2:
				# ネームリスト
				if re_natom.search(line):
					# 原子数の更新
					line = re.sub(r"\d+", str(atom), line)
				elif re_nf.search(line):
					# フラグメント数の更新
					line = re.sub(r"\d+", str(fragment), line)
				elif re_charge.search(line):
					# 電荷の更新
					line = re.sub(r"-?\d+", str(charge), line)
				elif re_autofrag.search(line):
					# autofrag の更新
					line = re.sub(r"=.+$", "='OFF'", line)

				if re_namelist_t.search(line):
					line = "/\n"
				elif not re_namelist_h.search(line):
					line = "  " + line

				namelists.append(line)
	return fragment_atoms, charges, BDAs, fragment_members, connections, namelists



# write_data (データの書き出し; ファイルが指定されていた場合はファイルに書き出し)
def write_data(contents, output):
	if output == None:
		for line in contents:
			print(line)
	else:
		if args.flag_overwrite == False:
			check_overwrite(output)
		with open(output, "w") as obj_output:
			for line in contents:
				obj_output.write("%s\n" % line)


# =============== main =============== #
if __name__ == '__main__':
	try:
		parser = argparse.ArgumentParser(description = "Fragment editor for mizuho ABINIT-MP", formatter_class=argparse.RawTextHelpFormatter)
		parser.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
		parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")

		subparser = parser.add_subparsers(help = "Sub-command")
		subparser.required = True

		parser_edit = subparser.add_parser("edit", help = "Convert ajf to fred (ajf -> fred)")
		parser_edit.set_defaults(func = "edit")
		parser_edit.add_argument("-p", "--pdb", metavar = "PDB", help = "Reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
		parser_edit.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
		parser_edit.add_argument("input", help = "ajf file")

		parser_rewrite = subparser.add_parser("rewrite", help = "Rewrite fred fred -> fred")
		parser_rewrite.set_defaults(func = "rewrite")
		parser_rewrite.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
		parser_rewrite.add_argument("input", help = "fred")

		parser_output = subparser.add_parser("output", help = "Convert fred to ajf (fred -> ajf)")
		parser_output.set_defaults(func = "output")
		parser_output.add_argument("-p", "--pdb", metavar = "PDB", help = "Reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
		parser_output.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
		parser_output.add_argument("input", help = "fred")

		parser_autofrag = subparser.add_parser("autofrag", help = "Auto fragmentation for PDB (pdb -> fred)")
		parser_autofrag.set_defaults(func = "autofrag")
		parser_autofrag.add_argument("-s", "--separate", action = "store_true", help = "Nucleotide is separates to base and sugar+phosphate")
		parser_autofrag.add_argument("-v", "--version", choices = ["3", "5", "m"], help = "ajf version: 3 = abinit-mp3, 5 = abinitmp5, m = mizuho")
		parser_autofrag.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
		parser_autofrag.add_argument("input", help = "PDB")

		parser_editfrag = subparser.add_parser("editfrag", help = "Create new fred in which fragments were devided based on PDB and fred (pdb + fred -> fred)")
		parser_editfrag.set_defaults(func = "editfrag")
		parser_editfrag.add_argument("-f", metavar = "fred", dest = "fred", required = True, help = "fred")
		parser_editfrag.add_argument("-b", metavar = "pdb", dest = "pdb", required = True, help = "pdb")
		parser_editfrag.add_argument("-n", metavar = "pdb", dest = "pdb", required = True, nargs = "+", help = "pdb for each fragments")
		parser_editfrag.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
		parser_editfrag.add_argument("input", help = "PDB")

		args = parser.parse_args()
	except TypeError:
		sys.stderr.write("ERROR: No sub-command (autofrag | edit | rewrite | output | editfrag)\n")
		sys.exit(1)

	check_file(args.input)

	if args.func == "edit":
		# 編集ファイルに変換

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = load_ajf(args.input, args.pdb)

		# 整形
		flag_fragment = 0
		output = []
		output.append("  FNo.  | charge | BDA | atoms of fragment")
		for i in range(len(fragment_atoms)):
			output.append("%7s |%6s  |%3s  |%s" % (i + 1, charges[i], BDAs[i], " ".join(map(lambda data : "%8d" % data, fragment_members[i]))))
		output.append("\n<< connections (ex. \"Next_fragment_atom    Prev_fragment_atom\")>>")
		for i in range(len(connections)):
			output.append("".join(map(lambda data : "%8d" % data, connections[i])))
		output.append("\n")
		output.append("===============< namelist >===============")
		for line in namelists:
			if re_natom.search(line):
				# 原子数の更新
				atom = sum(fragment_atoms)
				line = re.sub(r"\d+", str(atom), line)
			elif re_nf.search(line):
				# フラグメント数の更新
				fragment = len(fragment_members)
				line = re.sub(r"\d+", str(fragment), line)
			elif re_charge.search(line):
				# 電荷の更新
				charge = sum(charges)
				line = re.sub(r"-?\d+", str(charge), line)
			elif re_autofrag.search(line):
				# autofrag の更新
				line = re.sub(r"=.+$", "='OFF'", line)

			output.append(line)

		# 書き出し
		write_data(output, args.output_path)


	elif args.func == "rewrite":
		# リロード

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = load_fred(args.input)

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
		for line in namelists:
			if re_natom.search(line):
				# 原子数の更新
				atom = sum(fragment_atoms)
				line = re.sub(r"\d+", str(atom), line)
			elif re_nf.search(line):
				# フラグメント数の更新
				fragment = len(fragment_members)
				line = re.sub(r"\d+", str(fragment), line)
			elif re_charge.search(line):
				# 電荷の更新
				charge = sum(charges)
				line = re.sub(r"-?\d+", str(charge), line)
			elif re_autofrag.search(line):
				# autofrag の更新
				line = re.sub(r"=.+$", "='OFF'", line)
			output.append(line)

		# 書き出し
		write_data(output, args.output_path)

	elif args.func == "output":
		# ajf ファイルに変換

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = load_fred(args.input)

		file_reference = ""
		if args.pdb != None:
			check_file(args.pdb)
			file_reference = args.pdb
		else:
			for item in namelists:
				if "ReadGeom" in item:
					file_reference = item.strip()
					file_reference = re_wsp.sub("", file_reference)
					file_reference = file_reference.replace("ReadGeom=", "")
					file_reference = re_quote_h.sub("", file_reference)
					file_reference = re_quote_t.sub("", file_reference)
					check_file(file_reference)
					break

		check_charge(fragment_members, charges, file_reference)

		# 整形
		output = []
		flag_fragment = 0
		for line in namelists:
			if re_natom.search(line):
				# 原子数の更新
				atom = sum(fragment_atoms)
				line = re.sub(r"\d+", str(atom), line)
			elif re_nf.search(line):
				# フラグメント数の更新
				fragment = len(fragment_members)
				line = re.sub(r"\d+", str(fragment), line)
			elif re_charge.search(line):
				# 電荷の更新
				charge = sum(charges)
				line = re.sub(r"-?\d+", str(charge), line)
			elif re_autofrag.search(line):
				# autofrag の更新
				line = re.sub(r"=.+$", "='OFF'", line)
			elif re_fragment.search(line):
				# フラグメント情報書き出し場所を検索
				flag_fragment = 1
			elif flag_fragment == 1 and re_coord.search(line):
				# フラグメント情報の書き出し
				for line in convert_ajf(fragment_atoms, 8, 10):
					output.append(line)

				for line in convert_ajf(charges, 8, 10):
					output.append(line)

				for line in convert_ajf(BDAs, 8, 10):
					output.append(line)

				for j in range(len(fragment_members)):
					for line in convert_ajf(fragment_members[j], 8, 10):
						output.append(line)

				for j in range(len(connections)):
					for line in convert_ajf(connections[j], 8, 10):
						output.append(line)
				flag_fragment = 0
				continue
			output.append(line)

		# 出力
		write_data(output, args.output_path)

	elif args.func == "autofrag":
		pass


	elif args.func == "editfrag":
		pass
