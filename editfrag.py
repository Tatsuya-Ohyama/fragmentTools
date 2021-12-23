#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
editfrag.py
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse

from mods.func_prompt_io import check_exist, check_overwrite
from mods.func_string import target_range
from mods.MoleculeInformation import MoleculeInformation
from mods.FileFred import FileFred



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Create new fred in which fragments were devided based on PDB and fred (pdb + fred -> fred)", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.fred", required=True, help="fred file")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.fred", help="fred file")
	parser.add_argument("-p", dest="STRUCTURE_FILE", metavar="STRUCTURE.pdb", required=True, help="system structure file")
	parser.add_argument("-n", dest="FRAGMENT_STRUCTURE_LIST", metavar="FRAGMENT_STRCUTURE.pdb", required=True, nargs="+", help="fragment structure file")
	parser.add_argument("-t", dest="TARGET_FRAGMENT_LIST", metavar="TARGET_FRAGMENT", nargs="+", help="fragment indexes to which fragmentation is applied (separate with white space or specify by `1,2,5-10`)")
	parser.add_argument("-m", dest="FLAG_MULTI", action="store_true", default=False, help="applying separation to other the same type residues")
	parser.add_argument("-c", dest="CONNECTION_LIST", metavar = "ATOM1-ATOM2", nargs = "+", required=True, help = "Ambermask or atom index for connected atom pair list (Ex: :EG@C9-:EG@C10  34-25)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite_forcibly")

	args = parser.parse_args()


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
	fragment_index = obj_fred.fragments[-1].fragment_index + 1
	for structure_idx, fragment_structure_file in enumerate(args.FRAGMENT_STRUCTURE_LIST):
		# loop for new fragment structure

		# read fragment structure
		check_exist(fragment_structure_file, 2)
		fragment_structure = MoleculeInformation(fragment_structure_file)

		# change atom member in existing fragment
		cnt = 0
		new_fragment = []
		for obj_fragment in fragment_structure.output_fragmentdata("object", base_structure, args.FLAG_MULTI):
			atoms = set(obj_fragment.atoms)
			if args.TARGET_FRAGMENT_LIST is None or any([len(v & atoms) for v in target_fragment_atoms]):
				obj_fragment.set_fragment_index(fragment_index)
				fragment_index += 1
				obj_fred.add_fragment(obj_fragment)
				new_fragment.append(obj_fragment.fragment_index)
			cnt += 1
		cnt_total += cnt
		sys.stderr.write("Replace {0} fragments with {1} => Fragment {2}\n".format(cnt, fragment_structure_file, ", ".join([str(v) for v in new_fragment])))

	# add connection
	for str_connection in args.CONNECTION_LIST:
		pair = str_connection.split("-", maxsplit=1)
		if all([v.isdigit() for v in pair]):
			pair = [int(v) for v in pair]
		else:
			list_atom1 = base_structure.convert_number(pair[0])
			list_atom2 = base_structure.convert_number(pair[1])
			if len(list_atom1) != 1 or len(list_atom2) != 1:
				sys.stderr.write("ERROR: An atomic group cannot be specified as a connection pair. `{0}`\n".format(str_connection))
				sys.exit()
			pair = [int(list_atom1[0]), int(list_atom2[0])]

		for obj_fragment in obj_fred.fragments:
			if pair[1] in obj_fragment.atoms:
				obj_fragment.append_connection(pair)
				break

	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)
	obj_fred.write(args.OUTPUT_FILE)
