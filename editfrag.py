#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
editfrag.py
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
from basic_func import check_exist, check_overwrite
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
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.base_structure, 2)
	base_structure = MoleculeInformation(args.base_structure)

	check_exist(args.fred, 2)
	obj_fred = FredData(args.fred)

	for new_structure_file in args.add_structures:
		check_exist(new_structure_file, 2)
		new_structure = MoleculeInformation(new_structure_file)
		new_structure.update_index(base_structure)
		obj_fred.add_fragment(FragmentData(new_structure.output_fragmentdata()))
		sys.stderr.write("Added new fragment ({0}) to end of fragments\n".format(new_structure_file))

	if args.flag_overwrite == False:
		check_overwrite(args.output_file)
	obj_fred.write_file(args.output_file)
