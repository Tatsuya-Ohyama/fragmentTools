#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fragseparator.py

Program to write .pdb file for each fragment
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse

import parmed

from mods.func_prompt_io import *
from mods.FileFred import FileFred



# =============== functions =============== #


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="fragseparator", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-f", dest="FRED_FILE", metavar="INPUT.fred", help="Input fred file")
	parser.add_argument("-p", dest="PDB_FILE", metavar="INPUT.pdb", help="Input pdb file")
	parser.add_argument("-o", dest="OUTPUT_PREFIX", metavar="OUTPUT_PREFIX", help="Prefix for output")
	args = parser.parse_args()

	check_exist(args.FRED_FILE, 2)
	check_exist(args.PDB_FILE, 2)

	# fred ファイル読み込み
	obj_fred = FileFred().read(args.FRED_FILE)

	# PDB ファイル読み込み
	obj_mol = parmed.load_file(args.PDB_FILE)
	n_atoms = len(obj_mol.atoms)

	# 出力
	for obj_fragment in obj_fred.fragments:
		list_strip_flag = [False if i in obj_fragment.atoms else True for i in range(1, n_atoms+1)]
		obj_mol_output = obj_mol.copy(parmed.Structure)
		obj_mol_output.strip(list_strip_flag)

		output_file = "{}_{:03d}.pdb".format(args.OUTPUT_PREFIX, obj_fragment.index)
		obj_mol_output.write_pdb(output_file, renumber=True)
		sys.stderr.write("INFO: {0} was created.\n".format(output_file))
