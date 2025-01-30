#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import parmed


# =============== Function =============== #
def write_frag(output_prefix, list_obj_fragments, obj_mol):
	n_atoms = len(obj_mol.atoms)

	# 出力
	for obj_fragment in list_obj_fragments:
		list_strip_flag = [False if i in obj_fragment.atoms else True for i in range(1, n_atoms+1)]
		obj_mol_output = obj_mol.copy(parmed.Structure)
		obj_mol_output.strip(list_strip_flag)

		output_file = "{}_{:03d}.pdb".format(output_prefix, obj_fragment.index)
		obj_mol_output.write_pdb(output_file, renumber=True)
		sys.stderr.write("INFO: {0} was created.\n".format(output_file))
