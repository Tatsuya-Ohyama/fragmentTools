#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fred4 - fragment editor for mizuho ABINIT-MP
"""

import sys, os, re
import argparse
import functools

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
		user = sys.stdin.readline()
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

	return new_lists


# write_data (データの書き出し; ファイルが指定されていた場合はファイルに書き出し)
def write_data(contents, output):
	if output == None:
		for line in contents:
			print(line)
	else:
		check_overwrite(output)
		with open(output, "w") as obj_output:
			for line in contents:
				obj_output.write("%s\n" % line)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "Fragment editor for mizuho ABINIT-MP", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("input", help = "Input file")
	group = parser.add_mutually_exclusive_group(required = True)
	group.add_argument("--edit", "-E", action="store_true", default = False, help = "ajf -> fred_text")
	group.add_argument("--reload", "-R", action="store_true", default = False, help = "fred_text -> fred_text")
	group.add_argument("--output", "-O", action="store_true", default = False, help = "fred_text -> ajf")
	parser.add_argument("-o", metavar = "OUTPUT", dest = "output_path", help = "Output (Default: STDOUT)")
	parser.add_argument("-p", "--pdb", metavar = "PDB", help = "Reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
	args = parser.parse_args()

	check_file(args.input)

	re_wsp = re.compile(r"[\s\t]+")
	re_quote_h = re.compile(r"^['\"]")
	re_quote_t = re.compile(r"['\"]$")
	re_empty = re.compile(r"^[\s\t]*$")
	re_namelist_h = re.compile(r"^[\s\t]*\&")
	re_namelist_t = re.compile(r"^[\s\t]*/$")

	if args.edit == True:
		# 編集ファイルに変換
		re_namelist_fragment_h = re.compile(r"^[\s\t]*\&FRAGMENT")
		re_namelist_fragment_t = re.compile(r"^[\s\t]*\/")
		re_pdb_atom = re.compile(r"^((HETATM)|(ATOM))")

		flag_read = 0
		namelists = []
		all_atom = 0
		now_atom = 0
		fragment = 0
		now_fragment = 0
		fragment_atoms = []
		reference_pdb = ""
		fragment_atom_numbers = []
		charges = []
		BDAs = []
		now_fragment_atoms = []
		connections = []

		with open(args.input, "r") as obj_input:
			line_count = 0
			for line in obj_input:
				line_count += 1
				line = line.rstrip("\r\n")

				if re_empty.search(line):
					continue

				elif re_namelist_fragment_h.search(line):
					# FRAGMENT ネームリストの始端
					flag_read = 1
					namelists.append(line)
					continue

				elif 0 < flag_read:
					if re_namelist_fragment_t.search(line):
						# FRAGMENT ネームリストの終端
						namelists.append("{...}\n/\n")
						flag_read = 0
						continue

					datas = list(map(lambda data : int(data), split_n(line, 8)))

					if flag_read == 1:
						# フラグメント構成原子取得
						if "0" in datas:
							# 構成原子 0 のフラグメントがある場合
							sys.stderr.write("ERROR: Invalid ajf file in %d\n" % line_count)
							sys.stderr.write("       zero atoms in fragment found\n")
							sys.exit(1)

						fragment += len(datas)
						fragment_atom_numbers.extend(datas)

						now_atom += functools.reduce(lambda x, y: x + y, datas)
						if all_atom == now_atom:
							# PDB の総原子数と現在の原子数が一致した場合
							flag_read = 2
						elif all_atom < now_atom:
							# PDB の総原子数以上の原子を検出した場合
							sys.stderr.write("ERROR: The number of atoms in PDB file (%s) and ajf file mismatched (%d / %d)\n" % (reference_pdb, now_atom, all_atom))
							sys.stderr.write("       Maybe wrong PDB file was specified.\n")
							sys.stderr.write("       Please fix ReadGeom in ajf file, OR use -P option.\n")
							sys.exit(1)

					elif flag_read == 2:
						# 電荷情報取得
						now_fragment += len(datas)
						charges.extend(datas)

						if fragment == now_fragment:
							# フラグメント数が総フラグメント数と一致した場合
							flag_read = 3
							now_fragment = 0
						elif fragment < now_fragment:
							# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
							sys.stderr.write("ERROR: The number of fragments mismatched in charge section (%d / %d)\n" % (now_fragment, fragment))
							sys.exit(1)

					elif flag_read == 3:
						# BDA 取得
						now_fragment += len(datas)
						BDAs.extend(datas)

						if fragment == now_fragment:
							# フラグメント数が総フラグメント数と一致した場合
							flag_read = 4
							now_fragment = 0
						elif fragment < now_fragment:
							# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
							sys.stderr.write("ERROR: The number of fragments mismatched in BDA section\n" % (now_fragment, fragment))
							sys.exit(1)

					elif flag_read == 4:
						# フラグメント構成原子の原子順序番号取得
						now_fragment_atoms.extend(datas)

						if fragment_atom_numbers[now_fragment] <= len(now_fragment_atoms):
							now_fragment += 1
							fragment_atoms.append(now_fragment_atoms)
							now_fragment_atoms = []

							if fragment == now_fragment:
								# フラグメント数が総フラグメント数と一致した場合
								flag_read = 5
								now_fragment = 0
							elif fragment < now_fragment:
								# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
								sys.stderr.write("ERROR: The number of fragments mismatched in fragment atom section\n" % (now_fragment, fragment))
								sys.exit(1)

					elif flag_read == 5:
						# 接続情報取得
						line = line.strip()
						connections.append(list(map(lambda data : int(data),re_wsp.split(line))))

				else:
					if "ReadGeom" in line:
						if args.pdb != None:
							# 参照 PDB が指定されていた場合
							if not os.path.exists(args.pdb):
								sys.stderr.write("ERROR: No such PDB file (%s)\n" % args.pdb)
								sys.exit(1)
							else:
								reference_pdb = args.pdb
						else:
							# 参照 PDB が指定されていない場合
							reference_pdb = line.strip()
							reference_pdb = re_wsp.sub("", reference_pdb)
							reference_pdb = reference_pdb.replace("ReadGeom=", "")
							reference_pdb = re_quote_h.sub("", reference_pdb)
							reference_pdb = re_quote_t.sub("", reference_pdb)

						check_file(reference_pdb)

						with open(reference_pdb, "r") as refPDB:
							for p_line in refPDB:
								if re_pdb_atom.search(p_line):
									all_atom += 1

					elif re_namelist_fragment_t.search(line):
						line = "/\n"

					namelists.append(line)

		# 整形
		output = []
		output.append("  FNo.  | charge | BDA | atoms of fragment")
		for i in range(len(fragment_atom_numbers)):
			output.append("%7s |%6s  |%3s  |%s" % (i + 1, charges[i], BDAs[i], " ".join(map(lambda data : "%8d" % data, fragment_atoms[i]))))
		output.append("\n<< connections (ex. \"Next_fragment_atom    Prev_fragment_atom\")>>")
		for i in range(len(connections)):
			output.append("".join(map(lambda data : "%8d" % data, connections[i])))
		output.append("\n")
		output.append("===============< namelist >===============")
		for line in namelists:
			output.append(line)

		write_data(output, args.output_path)


	elif args.reload == True:
		# リロード
		flag_read = 0
		flag_namelist = 0
		fragment = 0
		charge = 0
		atom = 0
		output = []

		re_connection = re.compile(r"<<\sconnections")
		re_readmark = re.compile(r"=+<\snamelist\s>=+")
		re_data = re.compile(r"^[\s\t]*\d+[\s\t]*\|[\s\t]*-?\d+[\s\t]*\|[\s\t]*\d+[\s\t]*\|")
		re_namelist_ajf = re.compile(r"\&CNTRL", re.IGNORECASE)
		re_fmocntrl = re.compile(r"\&FMOCNTRL", re.IGNORECASE)
		re_natom = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
		re_charge = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
		re_nf = re.compile(r"^[\s\t]*NF", re.IGNORECASE)

		with open(args.input) as obj_input:
			for line in obj_input:
				line = line.strip()

				if re_empty.search(line):
					continue

				elif re_connection.search(line):
					# 接続情報取得フラグ
					flag_read = 1
					line = "\n<< connections (ex. \"Next_fragment_atom    Prev_fragment_atom\") >>"

				elif re_readmark.search(line):
					# ネームリスト取得フラグ
					flag_read = 2
					line = "\n\n===============< namelist >==============="

				elif flag_read == 0:
					# フラグメント情報取得
					if re_data.search(line):
						fragment += 1
						datas = line.split("|")
						datas = list(map(lambda data : data.strip(), datas))
						charge += int(datas[1])
						atoms = re_wsp.split(datas[3])
						atom += len(atoms)
						atoms = map(lambda data : "%8s" % data, atoms)
						line = "%7s |%6s  |%3s  |%s" % (datas[0], datas[1], datas[2], " ".join(atoms))

				elif flag_read == 1:
					# 接続情報取得
					datas = re_wsp.split(line)
					line = "".join(list(map(lambda data : "%8s" % data, datas)))

				elif flag_read == 2:
					# ネームリスト取得
					if re_namelist_t.search(line):
						line = "/\n"
					elif not re_namelist_h.search(line):
						line = "  " + line

					if re_namelist_ajf.search(line):
						flag_namelist = 1

					elif re_fmocntrl.search(line):
						flag_namelist = 2

					if flag_namelist == 1:
						if re_natom.search(line):
							line = re.sub(r"\d+", str(atom), line)
							sys.stderr.write("Atom:     %d\n" % atom)
							atom = "-"

						elif re_charge.search(line):
							line = re.sub(r"-?\d+", str(charge), line)
							sys.stderr.write("Charge:   %d\n" % charge)
							charge = "-"

						elif re_namelist_t.search(line):
							line = ""
							if atom != "-":
								line += "  Natom=%d\n" % atom
								sys.stderr.write("Atom:     %d\n" % atom)
							elif charge != "-":
								line += "  Charge=%d\n" % charge
								sys.stderr.write("Charge:   %d\n" % charge)
							line += "/\n"
							flag_namelist = 0

					elif flag_namelist == 2:
						# NF の修正
						if re_nf.search(line):
							line = re.sub(r"\d+", str(fragment), line)
							sys.stderr.write("Fragment: %d\n" % fragment)
							fragment = "-"
						elif re_namelist_t.search(line):
							line = ""
							if fragment != "-":
								line += "  NF=%d" % fragment
								sys.stderr.write("Fragment: %d\n" % fragment)
							line += "/"
							flag_namelist = 0
						elif re_namelist_t.search(line):
							line = "/\n"

				output.append(line)
		write_data(output, args.output_path)


	elif args.output == True:
		# ajf ファイルに変換
		flag_read = 0
		charges = []
		fragment = 0
		BDAs = []
		fragment_atoms = []
		connections = []
		output = []
		atom = 0
		fragment_atom_numbers = []
		flag_fragment = 0

		re_connection = re.compile(r"<< connections", re.IGNORECASE)
		re_readmark = re.compile(r"=+<\snamelist\s>=+")
		re_natom = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
		re_nf = re.compile(r"NF", re.IGNORECASE)
		re_charge = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
		re_autofrag = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
		re_fragment = re.compile(r"\&FRAGMENT", re.IGNORECASE)
		re_coord = re.compile(r"\{\.{3}\}")

		with open(args.input) as obj_input:
			line_count = 0

			for line in obj_input:
				line_count += 1
				line = line.strip()

				if line_count == 1:
					continue

				elif re_connection.search(line):
					# 接続情報
					flag_read = 1
					continue

				elif re_readmark.search(line):
					# この行以降がネームリスト
					flag_read = 2

					# フラグメントに含まれる原子数を計算
					for i in range(len(fragment_atoms)):
						fragment_atom_numbers.append(len(fragment_atoms[i]))
					fragment = len(fragment_atoms)

					# 電荷の合計算出
					charge = functools.reduce(lambda x, y: x + y, charges)

					continue

				elif re_empty.search(line):
					# 空行は無視
					continue

				if flag_read == 0:
					# フラグメント情報
					datas = line.split("|")
					datas = list(map(lambda data : data.strip(), datas))
					charges.append(int(datas[1]))
					BDAs.append(int(datas[2]))
					tmp = list(map(lambda data : int(data), re_wsp.split(datas[3])))
					fragment_atoms.append(tmp)
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
					elif re_fragment.search(line):
						# フラグメント情報書き出し場所を検索
						flag_fragment = 1
					elif flag_fragment == 1 and re_coord.search(line):
						# フラグメント情報の書き出し
						for line in convert_ajf(fragment_atom_numbers, 8, 10):
							output.append(line)

						for line in convert_ajf(charges, 8, 10):
							output.append(line)

						for line in convert_ajf(BDAs, 8, 10):
							output.append(line)

						for j in range(len(fragment_atoms)):
							for line in convert_ajf(fragment_atoms[j], 8, 10):
								output.append(line)

						for j in range(len(connections)):
							for line in convert_ajf(connections[j], 8, 10):
								output.append(line)
						continue

					if re_namelist_h.search(line):
						flag_output = 1

					elif re_namelist_t.search(line):
						flag_output = 2


					if flag_output == 0:
						output.append("  " + line)
					elif flag_output == 1:
						output.append(line)
						flag_output = 0
					elif flag_output == 2:
						output.append("/\n")
						flag_output = 0


		# 出力
		write_data(output, args.output_path)


	else:
		sys.stderr.write("ERROR: Unknown action")
		sys.exit(1)
