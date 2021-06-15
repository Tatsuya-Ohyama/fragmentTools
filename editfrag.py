#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
editfrag.py
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from classes.func_prompt_io import check_exist, check_overwrite
from classes.MoleculeInformation import MoleculeInformation
from classes.FragmentData import FragmentData
from classes.FredData import FredData

# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "editfrag.py", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-b", dest = "base_structure", metavar = "WHOLE.pdb", required = True, help = "whole structure for ABINIT-MP")
	parser.add_argument("-n", dest = "add_structures", metavar = "ADD.pdb", required = True, nargs = "+", help = "fragment structure")
	parser.add_argument("-f", dest = "fred", metavar = "BASE.fred", required = True, help = "original fred file")
	parser.add_argument("-o", dest = "output_file", metavar = "NEW.fred", required = True, help = "output fred")
	parser.add_argument("-c", dest = "connection", metavar = "ATOM1-ATOM2", nargs = "+", help = "connection list (Ex: :EG@C9-:EG@C10  34-25)")
	parser.add_argument("-m", dest = "flag_multi", action = "store_true", default = False, help = "applying separation to other the same type residues")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	# 構造ファイルの読み込み
	check_exist(args.base_structure, 2)
	base_structure = MoleculeInformation(args.base_structure)

	# fred ファイルの読み込み
	check_exist(args.fred, 2)
	obj_fred = FredData(args.fred)

	# 構造の追加
	cnt_total = 0
	for new_structure_file in args.add_structures:
		# 追加構造ファイルで回す
		# 構造の読み込み
		check_exist(new_structure_file, 2)
		new_structure = MoleculeInformation(new_structure_file)

		# 構造の追加
		cnt = 0
		for record in new_structure.output_fragmentdata(base_structure, args.flag_multi):
			obj_fred.add_fragment(FragmentData(record))
			cnt += 1
		cnt_total += cnt
		sys.stderr.write("Replace {0} fragments with {1}\n".format(cnt, new_structure_file))

	if args.connection is not None:
		re_num = re.compile(r'^\d+-\d+$')
		for connection in args.connection:
			connections = connection.split("-", 2)
			if re_num.search(connection):
				obj_fred.add_connection(connections)
			else:
				mask_list1 = base_structure.convert_number(connections[0])
				mask_list2 = base_structure.convert_number(connections[1])
				for mask1, mask2 in zip(mask_list1, mask_list2):
					obj_fred.add_connection([mask1, mask2])

	elif cnt_total % len(args.add_structures) == 0:
		obj_fred.add_connection(int(cnt_total / len(args.add_structures)))
	else:
		obj_fred.add_connection(1)

	if args.flag_overwrite == False:
		check_overwrite(args.output_file)
	obj_fred.write_file(args.output_file)
