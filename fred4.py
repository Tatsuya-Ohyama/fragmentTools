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
from mods.MoleculeInformation import MoleculeInformation
from mods.func_string import target_range



# =============== common variables =============== #
# general
RE_QUOTE_H = re.compile(r"^['\"]")
RE_QUOTE_T = re.compile(r"['\"]$")
RE_DIGIT = re.compile(r"[\d\s]+")
RE_CONNECT = re.compile(r'^\d+-\d+$')

ATOM_ELECTRON = {
	"H": 1, "Li": 3, "Be": 4, "B": 5,
	"C": 6, "N": 7, "O": 8, "F": 9,
	"Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
	"S": 16, "Cl": 17, "K": 19, "Ca": 20,
	"Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
	"Br": 35
}



# =============== functions =============== #
def check_electrons(list_obj_fragments, pdb_file):
	"""
	Function to check number of electrons

	Args:
		fragment_members (list): fragment information
		pdb_file (str): .pdb file
	"""
	element_table = {}

	# get atom_index and element
	with open(pdb_file, "r") as obj_input:
		for line_val in obj_input:
			if line_val.startswith("ATOM") or line_val.startswith("HETATM"):
				atom_idx = int(line_val[6:11].strip())
				element = RE_DIGIT.sub("", line_val[12:14].strip())
				element = RE_QUOTE_H.sub("", element)
				element = RE_QUOTE_T.sub("", element)
				if element in ["HO", "HH"]:
					element = "H"

				element_table[atom_idx] = element
				if not element in ATOM_ELECTRON:
					sys.stderr.write("ERROR: Unknown atomtype ({0}). Skipped...\n".format(element))

	flag_error = False
	for fragment_idx, obj_fragment in enumerate(list_obj_fragments, 1):
		electron_fragment = 0
		for atom_i in obj_fragment.atoms:
			try:
				electron_atom = ATOM_ELECTRON[element_table[atom_i]]
				electron_fragment += electron_atom
			except ValueError:
				sys.stderr.write("ERROR: atom index `{0}` is not found in list.\n".format(atom_i))
				sys.exit(1)

		electron_fragment += (-1 * obj_fragment.charge)
		if electron_fragment % 2 != 0:
			sys.stderr.write("ERROR: Invalid number of fragment electrons.\n")
			sys.stderr.write("       The number of electrons in fragment No. {0} is {1}.\n".format(fragment_idx, electron_fragment))
			flag_error = True

	if flag_error:
		sys.stderr.write("ERROR: The number of electrons for some fragments were not even number.\n")
		sys.stderr.write("       Prceeding? (y/N): ")
		sys.stderr.flush()
		user = sys.stdin.readline().strip()
		if user.lower() != "y":
			sys.exit(0)
	else:
		sys.stderr.write("INFO: check_electrons is ok.\n")



# =============== main =============== #
if __name__ == '__main__':
	try:
		parser = argparse.ArgumentParser(description="Fragment editor for mizuho ABINIT-MP", formatter_class=argparse.RawTextHelpFormatter)

		subparser = parser.add_subparsers(help="subcommand")
		subparser.required = True

		parser_edit = subparser.add_parser("edit", help="Convert ajf to fred (ajf -> fred)")
		parser_edit.set_defaults(func="edit")
		parser_edit.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.ajf", required=True, help="ajf file")
		parser_edit.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.fred", required=True, help="fred file")
		parser_edit.add_argument("-p", "--pdb", dest="PDB_FILE", metavar="REF.pdb", help="reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
		parser_edit.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_rewrite = subparser.add_parser("rewrite", help="Rewrite fred fred -> fred")
		parser_rewrite.set_defaults(func="rewrite")
		parser_rewrite.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.fred", required=True, help="fred file")
		parser_rewrite.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.fred", required=True, help="fred file")
		parser_rewrite.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_output = subparser.add_parser("output", help="Convert fred to ajf (fred -> ajf)")
		parser_output.set_defaults(func="output")
		parser_output.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.fred", required=True, help="fred file")
		parser_output.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.ajf", required=True, help="ajf file")
		parser_output.add_argument("-p", "--pdb", dest="PDB_FILE", metavar="REF.pdb", help="Reference PDB (if not specify, this program use ReadGeom PDB in ajf file)")
		parser_output.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_autofrag = subparser.add_parser("autofrag", help="Auto fragmentation for PDB (pdb -> fred)")
		parser_autofrag.set_defaults(func="autofrag")
		parser_autofrag.add_argument("-p", dest="INPUT_FILE", metavar="INPUT.pdb", required=True, help="structure file")
		parser_autofrag.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.fred", required=True, help="fred file")
		parser_autofrag.add_argument("-s", "--separate", action="store_true", help="Nucleotide is separates to base and sugar+phosphate")
		parser_autofrag.add_argument("-v", "--version", choices=["3", "5", "m"], help="ajf version: 3 = abinit-mp3, 5 = abinitmp5, m = mizuho")
		parser_autofrag.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		parser_editfrag = subparser.add_parser("editfrag", help="Create new fred in which fragments were devided based on PDB and fred (pdb + fred -> fred)")
		parser_editfrag.set_defaults(func="editfrag")
		parser_editfrag.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.fred", required=True, help="fred file")
		parser_editfrag.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.fred", help="fred file")
		parser_editfrag.add_argument("-p", dest="STRUCTURE_FILE", metavar="STRUCTURE.pdb", required=True, help="system structure file")
		parser_editfrag.add_argument("-n", dest="FRAGMENT_STRUCTURE_LIST", metavar="FRAGMENT_STRCUTURE.pdb", required=True, nargs="+", help="fragment structure file")
		parser_editfrag.add_argument("-t", dest="TARGET_FRAGMENT_LIST", metavar="TARGET_FRAGMENT", nargs="+", help="fragment indexes to which fragmentation is applied (separate with white space or specify by `1,2,5-10`)")
		parser_editfrag.add_argument("-m", dest="FLAG_MULTI", action="store_true", default=False, help="applying separation to other the same type residues")
		parser_editfrag.add_argument("-c", dest="CONNECTION_LIST", metavar = "ATOM1-ATOM2", nargs = "+", required=True, help = "connection list described by Ambermask (Ex: :EG@C9-:EG@C10  34-25)")
		parser_editfrag.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

		args = parser.parse_args()

	except TypeError:
		sys.stderr.write("ERROR: No sub-command (autofrag | edit | rewrite | output | editfrag)\n")
		sys.exit(1)


	check_exist(args.INPUT_FILE, 2)

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
		obj_fred.set_parameters(obj_ajf.parameters)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_fred.write(args.OUTPUT_FILE)


	elif args.func == "rewrite":
		# rewrite mode

		# read .fred file
		obj_fred = FileFred()
		obj_fred.read(args.INPUT_FILE)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_fred.write(args.OUTPUT_FILE)


	elif args.func == "output":
		# output mode (convert to ajf)

		# read .fred file
		obj_fred = FileFred()
		obj_fred.read(args.INPUT_FILE)

		# convert .ajf file
		obj_ajf = FileAJF()
		obj_ajf.set_parameters(obj_fred.complete_parameters)

		file_reference = obj_fred.parameters["&CNTRL"]["ReadGeom"]
		if args.PDB_FILE is not None:
			file_reference = args.PDB_FILE
		check_electrons(obj_fred.fragments, file_reference)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_ajf.write(args.OUTPUT_FILE)


	elif args.func == "autofrag":
		# autofrag mode
		pass


	elif args.func == "editfrag":
		# editfrag mode

		# read structure
		check_exist(args.STRUCTURE_FILE, 2)
		base_structure = MoleculeInformation(args.STRUCTURE_FILE)

		# read .fred
		obj_fred = FileFred().read(args.INPUT_FILE)

		# create target fragment list
		target_fragment_atoms = []
		if args.TARGET_FRAGMENT_LIST is not None:
			target_fragment = [[int(v)] if v.isdigit() else target_range(v, start=1) for v in args.TARGET_FRAGMENT_LIST]
			target_fragment = [v2 for v1 in target_fragment for v2 in v1]
			target_fragment = list(sorted(set(target_fragment)))
			target_fragment_atoms = [set(v.atoms) for i, v in enumerate(obj_fred.fragments, 1) if i in target_fragment]

		cnt_total = 0
		for structure_idx, fragment_structure_file in enumerate(args.FRAGMENT_STRUCTURE_LIST):
			# loop for new fragment structure

			# read fragment structure
			check_exist(fragment_structure_file, 2)
			fragment_structure = MoleculeInformation(fragment_structure_file)

			# change atom member in existing fragment
			cnt = 0
			for obj_fragment in fragment_structure.output_fragmentdata("object", base_structure, args.FLAG_MULTI):
				atoms = set(obj_fragment.atoms)
				if args.TARGET_FRAGMENT_LIST is None or any([len(v & atoms) for v in target_fragment_atoms]):
					obj_fred.add_fragment(obj_fragment)
				cnt += 1
			cnt_total += cnt
			sys.stderr.write("Replace {0} fragments with {1}\n".format(cnt, fragment_structure_file))

		# add connection
		for str_connection in args.CONNECTION_LIST:
			atom1, atom2 = str_connection.split("-", maxsplit=1)
			if RE_CONNECT.search(str_connection):
				# direct format of connection (atom indexes)
				obj_fred.add_connection([atom1, atom2])
			else:
				mask_list1 = base_structure.convert_number(atom1)
				mask_list2 = base_structure.convert_number(atom2)
				for mask1, mask2 in zip(mask_list1, mask_list2):
					obj_fred.add_connection([mask1, mask2])

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)
		obj_fred.write(args.OUTPUT_FILE)
