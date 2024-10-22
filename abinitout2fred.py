#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Program to convert .out (.log) file for ABINIT-MP to .fred file
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse

from mods.func_prompt_io import check_exist, check_overwrite
from mods.FragmentData  import FragmentData
from mods.FileAJF import FileAJF
from mods.FileFred import FileFred



# =============== function =============== #
def read_out_file(obj_fred, out_file):
	"""
	Function to add fragment information to fred object

	Args:
		obj_fred (FileFred object): FileFred object
		out_file (str): .out (.log) file for ABINIT-MP

	Returns:
		FileFred object
	"""
	read_section = None
	column_pos = None
	list_fragments = []
	list_fragment_idx = []
	obj_fragment_current = None
	with open(out_file, "r") as obj_input:
		for line_val in obj_input:
			if "Seq." in line_val and  "Frag." in line_val:
				# Start Residue Level in Auto-fragmentation
				column_pos = {}
				column_pos["Frag."] = [line_val.index("Frag."), None]
				if "Residue" in line_val:
					# Amino acid
					read_section = "auto-frag residue amino"
					column_pos["Frag."][1] = line_val.index("Residue")
					column_pos["Charge"] = [line_val.index("Charge"), None]
					column_pos["N-term."] = [line_val.index("N-term."), line_val.index("C-Term.")]
					column_pos["C-term."] = [line_val.index("C-Term."), line_val.index("Charge")]

				elif "Base" in line_val:
					# Nucleotide
					read_section = "auto-frag residue nucleotide"
					column_pos["Frag."][1] = line_val.index("Base")
					column_pos["Charge"] = [line_val.index("Formal charge"), None]
					column_pos["Base"] = [line_val.index("Base"), line_val.index("5'-term.")]
					column_pos["5'-term."] = [line_val.index("5'-term."), line_val.index("3'-Term.")]
					column_pos["3'-term."] = [line_val.index("3'-Term."), line_val.index("Formal charge")]
				continue

			if "Frag." in line_val and "ATOM" in line_val:
				# Start Atom Level in Auto-fragmentation
				read_section = "auto-frag atom"
				column_pos = {}
				column_pos["Frag."] = [line_val.index("Frag."), line_val.index("Elec")]
				column_pos["Atom"] = [line_val.index("ATOM"), None]
				list_fragment_idx = [obj_fragment.fragment_index for obj_fragment in list_fragments]
				continue

			if "Frag." in line_val and "Bonded Atom" in line_val:
				# Start connection
				read_section = "connection"
				column_pos = {}
				column_pos["Frag."] = [line_val.index("Frag."), line_val.index("Bonded Atom")]
				column_pos["Bonded Atom"] = [line_val.index("Bonded Atom"), line_val.index("Proj.")]
				continue

			if len(line_val.strip()) == 0:
				# End of section (Blank line)
				read_section = None
				column_pos = None
				continue

			if read_section is None:
				continue

			if "auto-frag residue" in read_section:
				# Residue Level in Auto-fragmentation
				obj_fragment = FragmentData()
				obj_fragment.set_fragment_index(int(line_val[column_pos["Frag."][0]:column_pos["Frag."][1]].strip()))
				if read_section == "auto-frag residue amino":
					# Amino acid
					obj_fragment.set_charge(int(line_val[column_pos["Charge"][0]:column_pos["Charge"][1]].strip().split()[0]))

				elif read_section == "auto-frag residue nucleotide":
					# Nucleotide
					obj_fragment.set_charge(int(line_val[column_pos["Charge"][0]:column_pos["Charge"][1]].strip()))

				list_fragments.append(obj_fragment)
				continue

			if read_section == "auto-frag atom":
				# Atom Level in Auto-fragmentation
				if "Invalid" in line_val or "Fragment No." in line_val:
					sys.stderr.write(line_val)
					continue

				fragment_idx = line_val[column_pos["Frag."][0]:column_pos["Frag."][1]].strip()
				if len(fragment_idx) != 0:
					fragment_idx = int(fragment_idx)
					list_idx = list_fragment_idx.index(fragment_idx)
					obj_fragment_current = list_fragments[list_idx]
				atoms = [int(v) for v in line_val[column_pos["Atom"][0]:column_pos["Atom"][1]].strip().split()]
				obj_fragment_current.set_atoms(obj_fragment_current.atoms + atoms)
				continue

			if read_section == "connection":
				fragment_idx = int(line_val[column_pos["Frag."][0]:column_pos["Frag."][1]].strip())
				list_idx = list_fragment_idx.index(fragment_idx)
				obj_fragment_current = list_fragments[list_idx]
				connections = [int(v) for v in line_val[column_pos["Bonded Atom"][0]:column_pos["Bonded Atom"][1]].strip().split()]
				obj_fragment_current.append_connection(connections)

	obj_fred.set_fragments(list_fragments)
	return obj_fred


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Program to convert .out (.log) file for ABINIT-MP to .fred file", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest="INPUT_OUT_FILE", metavar="INPUT.out", required=True, help=".out (.log) file for ABINIT-MP")
	parser.add_argument("-a", dest="INPUT_AJF_FILE", metavar="INPUT.fred", required=True, help=".fred file")
	parser.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.fred", required=True, help=".fred file")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	args = parser.parse_args()

	check_exist(args.INPUT_OUT_FILE, 2)
	check_exist(args.INPUT_AJF_FILE, 2)

	obj_ajf = FileAJF()
	obj_ajf.read(args.INPUT_AJF_FILE)

	# convert .ajf file
	obj_fred = FileFred()
	obj_fred.set_parameters(obj_ajf.parameters)

	obj_fred = read_out_file(obj_fred, args.INPUT_OUT_FILE)

	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)
	obj_fred.write(args.OUTPUT_FILE)
