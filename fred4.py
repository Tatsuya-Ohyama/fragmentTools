#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fred4 - fragment editor for mizuho ABINIT-MP
"""

import sys, os, re
import argparse
sys.path.append("modules")
import modules.basic

# =============== variables =============== #
re_wsp = re.compile(r"[\s\t]+")
re_quote_h = re.compile(r"^['\"]")
re_quote_t = re.compile(r"['\"]$")


# =============== functions =============== #
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


# write_data (データの書き出し; ファイルが指定されていた場合はファイルに書き出し)
def write_data(contents, output):
	if output is None:
		for line in contents:
			print(line)
	else:
		if args.flag_overwrite == False:
			modules.basic.check_overwrite(output)
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

	modules.basic.check_exist(args.input, 2)

	if args.func == "edit":
		# 編集ファイルに変換

		import modules.fred

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = modules.fred.load_ajf(args.input, args.pdb)

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

		import modules.fred

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = modules.fred.load_fred(args.input)

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

		import modules.fred

		# 読み込み
		(fragment_atoms, charges, BDAs, fragment_members, connections, namelists) = modules.fred.load_fred(args.input)

		file_reference = ""
		if args.pdb is not None:
			modules.basic.check_exist(args.pdb, 2)
			file_reference = args.pdb
		else:
			for item in namelists:
				if "ReadGeom" in item:
					file_reference = item.strip()
					file_reference = re_wsp.sub("", file_reference)
					file_reference = file_reference.replace("ReadGeom=", "")
					file_reference = re_quote_h.sub("", file_reference)
					file_reference = re_quote_t.sub("", file_reference)
					modules.basic.check_exist(file_reference, 2)
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
				for line in modules.fred.convert_ajf(fragment_atoms, 8, 10):
					output.append(line)

				for line in modules.fred.convert_ajf(charges, 8, 10):
					output.append(line)

				for line in modules.fred.convert_ajf(BDAs, 8, 10):
					output.append(line)

				for j in range(len(fragment_members)):
					for line in modules.fred.convert_ajf(fragment_members[j], 8, 10):
						output.append(line)

				for j in range(len(connections)):
					for line in modules.fred.convert_ajf(connections[j], 8, 10):
						output.append(line)
				flag_fragment = 0
				continue
			output.append(line)

		# 出力
		write_data(output, args.output_path)

	elif args.func == "autofrag":
		import modules.autofrag
		pass


	elif args.func == "editfrag":
		import modules.editfrag
		pass
