#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fred4 - fragment editor
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import os
import re

from mods.func_prompt_io import *
from mods.FileAJF import FileAJF
from mods.FileFred import FileFred
from mods.FragmentData import FragmentData



# =============== common variables =============== #
# general
RE_WSP = re.compile(r"[\s\t]+")
RE_QUOTE_H = re.compile(r"^['\"]")
RE_QUOTE_T = re.compile(r"['\"]$")
RE_EMPTY = re.compile(r"^[\s\t]*$")
RE_NAMELIST_H = re.compile(r"^[\s\t]*\&")
RE_NAMELIST_T = re.compile(r"^[\s\t]*/$")
RE_DIGIT = re.compile(r"[\d\s]+")

RE_NATOM = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
RE_NF = re.compile(r"NF", re.IGNORECASE)
RE_CHARGE = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
RE_AUTOFRAG = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
RE_FRAGMENT = re.compile(r"\&FRAGMENT", re.IGNORECASE)
RE_COORD = re.compile(r"\{\.{3}\}")

RE_NAMELIST_FRAGMENT_H = re.compile(r"^[\s\t]*\&FRAGMENT")
RE_NAMELIST_FRAGMENT_T = re.compile(r"^[\s\t]*\/")

ATOM_CHARGES = {
	"H": 1, "Li": 3, "Be": 4, "B": 5,
	"C": 6, "N": 7, "O": 8, "F": 9,
	"Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
	"S": 16, "Cl": 17, "K": 19, "Ca": 20,
	"Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
	"Br": 35
}


# =============== functions =============== #
def check_electrons(fragment_members, charges, pdb):
	"""
	Function to check system charge

	Args:
		fragment_members (list): fragment information
		charges (list): fragment charge
		pdb (str): .pdb file
	"""
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
		sys.stderr.write("INFO: check_electrons is ok.\n")


def convert_ajf(fragment_info, width, n):
	"""
	Function to convert fragment information to .ajf file

	Args:
		fragment_info (list): fragment information
		width (int): width
		n (int): number of elements

	Returns:
		list
	"""
	fragment_info = list(map(lambda data : str(data), fragment_info))
	new_fragment_info = [""]

	format_data = "{" + str(width) + "}"

	index = 0
	for i in range(len(fragment_info)):
		new_fragment_info[index] += format_data.format(fragment_info[i])
		if i != 0 and (i + 1) % n == 0:
			index += 1
			new_fragment_info.append("")

	if len(new_fragment_info[len(new_fragment_info) - 1]) == 0:
		# 指定された個数でちょうど終わる場合は、無駄な要素が追加されるので削除する
		del(new_fragment_info[len(new_fragment_info) - 1])

	return new_fragment_info




# =============== main =============== #
if __name__ == '__main__':
	try:
		parser = argparse.ArgumentParser(description="Fragment editor for mizuho ABINIT-MP", formatter_class=argparse.RawTextHelpFormatter)

		subparser = parser.add_subparsers(help="subcommand")
		subparser.required = True

		parser_edit = subparser.add_parser("edit", help="Convert ajf to fred (ajf -> fred)")
		parser_edit.set_defaults(func="edit")
		parser_edit.add_argument("-i", dest="INPUT_FILE", metavar="INPUT", required=True, help="ajf file")
		parser_edit.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT", help="output file")
		parser_edit.add_argument("-p", "--pdb", dest="PDB_FILE", metavar="PDB", help="reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
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

	obj_fred = FragmentData()
	if args.func == "edit":
		# edit mode (convert to edit-style)

		# read .ajf file
		obj_ajf = FileAJF()
		obj_ajf.read(args.INPUT_FILE)
		list_fragments, list_connections = obj_ajf.create_fragment_objects(args.PDB_FILE)

		# convert .fred file
		obj_fred = FileFred()
		obj_fred.set_n_atom(sum([len(v.atoms) for v in list_fragments]))
		obj_fred.set_charge(sum([v.charge for v in list_fragments]))
		obj_fred.set_fragments(list_fragments)
		obj_fred.set_connections(list_connections)
		obj_fred.set_other_info(obj_ajf.parameter_list(True))

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_fred.write(args.OUTPUT_FILE)


	elif args.func == "rewrite":
		# rewrite mode

		obj_fred = FileFred()
		obj_fred.read(args.INPUT_FILE)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_fred.write(args.OUTPUT_FILE)


	elif args.func == "output":
		# output mode (convert to ajf)

		obj_fred = FileFred()
		obj_fred.read(args.INPUT_FILE)

		obj_ajf = FileAJF()
		obj_ajf.set_parameters()

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_ajf.write(args.OUTPUT_FILE)

		check_electrons(fragment_members, charges, file_reference)

		# # 整形
		# output = []
		# flag_fragment = 0
		# for line_val in namelists:
		# 	if RE_NATOM.search(line_val):
		# 		# 原子数の更新
		# 		atom = sum(fragment_atoms)
		# 		line_val = re.sub(r"\d+", str(atom), line_val)
		# 	elif RE_NF.search(line_val):
		# 		# フラグメント数の更新
		# 		fragment = len(fragment_members)
		# 		line_val = re.sub(r"\d+", str(fragment), line_val)
		# 	elif RE_CHARGE.search(line_val):
		# 		# 電荷の更新
		# 		charge = sum(charges)
		# 		line_val = re.sub(r"-?\d+", str(charge), line_val)
		# 	elif RE_AUTOFRAG.search(line_val):
		# 		# autofrag の更新
		# 		line_val = re.sub(r"=.+$", "='OFF'", line_val)
		# 	elif RE_FRAGMENT.search(line_val):
		# 		# フラグメント情報書き出し場所を検索
		# 		flag_fragment = 1
		# 	elif flag_fragment == 1 and RE_COORD.search(line_val):
		# 		# フラグメント情報の書き出し
		# 		continue
		# 	output.append(line_val)


	elif args.func == "autofrag":
		pass


	elif args.func == "editfrag":
		pass
