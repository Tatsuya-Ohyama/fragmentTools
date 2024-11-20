#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import parmed
from mods.FragmentData import FragmentData


# =============== Constant =============== #
RESIDUE_TYPES = {
	"AminoAcid": ["ACE", "ALA", "ARG", "ASN", "ASP", "CYS", "CYM", "GLN", "GLU", "GLH", "GLY", "HIS", "HIP", "HID", "HIE", "ILE", "LEU", "LYS", "MET", "NME", "PHE", "PRO", "SER", "SYM", "THR", "TRP", "TYR", "VAL"],
	"NucleicAcid": ["DA5", "DT5", "DG5", "DC5", "DA3", "DT3", "DG3", "DC3", "DA", "DT", "DG", "DC", "RA5", "RU5", "RG5", "RC5", "RA3", "RU3", "RG3", "RC3", "RA", "RU", "RG", "RC", "DNA_base"],
	"Water": ["SOL", "WAT", "HOH"],
	"Ion": ["Na", "Mg", "K", "Ca", "Cl", "Zn"]
}



# =============== Function =============== #
def fragmentation(structure_file, sep_amino="+amino", sep_nuc="+base"):
	"""
	Function to separate molecule into fragments

	Args:
		structure_file (str): structure file path
		sep_amino (str, optional): separation option for amino acids ("+amino", "/amino", "+peptide", or "/peptide" ;Default: "+amino")
		sep_nuc (str, optional): separation option for nucleic acids ("+base", "/base", or "/sugar"; Default: "+base")

	Returns:
		list: [[obj_atom, ...], ...]	# list of obj_atom list for each fragment
	"""
	# check arguments
	if sep_amino not in ["+amino", "/amino", "+peptide", "/peptide"]:
		sys.stderr.write("ERROR: `sep_amino` should be `+amino`, `/amino`, `+peptide`, or `/peptide`.\n")
		sys.exit(1)

	if sep_nuc not in ["+base", "/base", "/sugar"]:
		sys.stderr.write("ERROR: `sep_nuc` should be `+base`, `/base`, or `/sugar`.\n")
		sys.exit(1)

	obj_mol = parmed.load_file(structure_file)
	if obj_mol.atoms[-1].idx+1 != obj_mol.atoms[-1].number:
		sys.stderr.write("WARNING: The index of atoms in the PDB file does not match the index by number of atoms.\n         The output file is indexed by number of atoms.\n")

	list_obj_fragments = [[FragmentData().set_atoms(list(obj_residue.atoms))] for obj_residue in obj_mol.residues]
	list_atom_info = {obj_atom: None for obj_atom in obj_mol.atoms}
	list_obj_atom_shift = [[[]] for obj_residue in obj_mol.residues]
	for frag_i, obj_residue in enumerate(obj_mol.residues):
		res_type = None
		for res_type_ref, list_resnames in RESIDUE_TYPES.items():
			if obj_residue.name in list_resnames:
				res_type = res_type_ref
				break

		if res_type == "AminoAcid":
			# determine fragment charge
			charge = None
			if obj_residue.name in ["LYS", "ARG", "HIP", "SYM"]:
				charge = 1

			elif obj_residue.name in ["ASP", "GLU"]:
				charge = -1

			elif obj_residue.name in ["ACE", "ALA", "ASN", "CYS", "GLN", "GLY", "HIS", "HID", "HIE", "ILE", "LEU", "MET", "NME", "PHE", "PRO", "SER", "SYM", "THR", "TRP", "TYR", "VAL"]:
				charge = 0

			else:
				sys.stderr.write("WARNING: Unable to determine fragment charge by undefined amino acid residues (`{}`).\n".format(obj_residue.name))


			# determine terminal
			list_atoms = set(obj_residue.atoms)
			list_bond_partners = [set(obj_atom.bond_partners) - list_atoms for obj_atom in obj_residue.atoms]
			list_bond_partners = [v for v in list_bond_partners if len(v) != 0]
			term_type = None

			if len(list_bond_partners) == 0:
				# single residue
				term_type = "S"

			elif len(list_bond_partners) == 1:
				# bond with one residue -> terminal
				obj_atom_partner = list_bond_partners.pop().pop()
				if frag_i == 0:
					# N-terminal (first fragment)
					term_type = "N"
					obj_atom_N = [obj_atom for obj_atom in obj_residue.atoms if obj_atom.name == "N"][0]
					list_obj_atom_H = [obj_atom for obj_atom in obj_atom_N.bond_partners if obj_atom.element == 1]
					if len(list_obj_atom_H) == 3:
						list_obj_fragments[frag_i][0].add_charge(1)

				elif obj_atom_partner not in obj_mol.residues[frag_i-1].atoms:
					# N-terminal (not connected previous residue)
					term_type = "N"
					obj_atom_N = [obj_atom for obj_atom in obj_residue.atoms if obj_atom.name == "N"][0]
					list_obj_atom_H = [obj_atom for obj_atom in obj_atom_N.bond_partners if obj_atom.element == 1]
					if len(list_obj_atom_H) == 3:
						list_obj_fragments[frag_i][0].add_charge(1)

				else:
					# C-terminal
					term_type = "C"
					obj_atom_C = [obj_atom for obj_atom in obj_residue.atoms if obj_atom.name == "C"][0]
					list_obj_atom_O = [obj_atom for obj_atom in obj_atom_C.bond_partners if obj_atom.element == 8]
					list_bond_partners_O = [set(obj_atom.bond_partners) - set([obj_atom_C]) for obj_atom in list_obj_atom_O]
					list_bond_partners_O = [v for v in list_bond_partners_O if len(v) != 0]
					if len(list_bond_partners_O) == 0:
						list_obj_fragments[frag_i][0].add_charge(-1)

			if sep_amino.startswith("+"):
				# prepare fragment
				list_obj_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
				list_obj_fragments[frag_i][0].add_charge(charge)

			elif sep_amino.startswith("/"):
				# prepare fragment
				list_obj_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
				list_obj_fragments[frag_i][0].add_charge(0)
				list_obj_fragments[frag_i].append(FragmentData())	# add side chain fragment
				list_obj_atom_shift[frag_i].append([])
				list_obj_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_amino, "Sidechain"))
				list_obj_fragments[frag_i][1].add_charge(charge)

				# separate into side chain
				list_main_chain = []
				list_side_chain = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if obj_residue.name in ["GLY", "PRO"]:
						# do not separate side chain for GLY
						list_main_chain.append(obj_atom)

					elif obj_atom.name in ["N", "H", "C", "O", "CA", "HA"]:
						list_main_chain.append(obj_atom)

					elif term_type == "N":
						# N-terminal and atoms connected with "N" -> main chain
						for obj_atom_partner in obj_atom.bond_partners:
							if obj_atom_partner.name == "N" and obj_residue.number == obj_atom_partner.residue.number:
								# intra N-H in N-terminal -> main chain
								list_main_chain.append(obj_atom)
								break

						else:
							list_side_chain.append(obj_atom)

					elif term_type == "C":
						# C-terminal and atoms connected with "C" -> main chain
						for obj_atom_partner in obj_atom.bond_partners:
							if obj_atom_partner.name == "C" and obj_residue.number == obj_atom_partner.residue.number:
								# intra C-O in C-terminal (OXT) -> main chain
								list_main_chain.append(obj_atom)
								break

						else:
							list_side_chain.append(obj_atom)

					else:
						list_side_chain.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_main_chain)
				list_obj_fragments[frag_i][1].set_atoms(list_side_chain)


			if sep_amino in ["+amino", "/amino"]:
				# shift atoms in main chain fragments
				list_remain_atoms = []
				list_shift_atoms = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if obj_atom.name in ["C", "O"] and term_type != "C":
						# shift to next fragment
						list_shift_atoms.append(obj_atom)

					else:
						list_remain_atoms.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_remain_atoms)
				if len(list_shift_atoms) != 0:
					if frag_i+1 < len(list_obj_atom_shift):
						list_obj_atom_shift[frag_i+1][0].extend(list_shift_atoms)

					else:
						# end of shift atom
						sys.stderr.write("ERROR: Algorithm error in amino acid fragmentation. Please report to developer.\n")
						sys.exit(1)
			continue


		if res_type == "NucleicAcid":
			# determine terminal
			list_atoms = set(obj_residue.atoms)
			list_bond_partners = [set(obj_atom.bond_partners) - list_atoms for obj_atom in obj_residue.atoms]
			list_bond_partners = [v for v in list_bond_partners if len(v) != 0]
			term_type = None

			if len(list_bond_partners) == 0:
				# single residue
				term_type = "S"

			elif len(list_bond_partners) == 1:
				# bond with one residue -> terminal
				obj_atom_partner = list_bond_partners.pop().pop()
				if frag_i == 0:
					# 5'-terminal (first fragment)
					term_type = "5"

				elif obj_atom_partner not in obj_mol.residues[frag_i-1].atoms:
					# N-terminal (not connected previous residue)
					term_type = "5"

				else:
					# C-terminal
					term_type = "3"

			if sep_nuc == "+base":
				# prepare fragment
				list_obj_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
				list_obj_fragments[frag_i][0].add_charge(0)

				# shift atoms in main chains
				list_remain_atoms = []
				list_shift_atoms = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if term_type in [None, "3"] and obj_atom.name in ["P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2", "1H5'", "2H5'"]:
						list_shift_atoms.append(obj_atom)

					else:
						list_remain_atoms.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_remain_atoms)
				if len(list_shift_atoms) != 0:
					if 0 < frag_i:
						list_obj_atom_shift[frag_i-1][0].extend(list_shift_atoms)
						if len([obj_atom for obj_atom in list_shift_atoms if obj_atom.name == "P"]) != 0:
							list_obj_fragments[frag_i-1][0].add_charge(-1)

					else:
						# end of shift atom
						sys.stderr.write("ERROR: Algorithm error in nucleic acid fragmentation (+base). Please report to developer.\n")
						sys.exit(1)

			elif sep_nuc == "/base":
				# prepare fragment
				list_obj_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
				list_obj_fragments[frag_i][0].add_charge(0)
				list_obj_fragments[frag_i].append(FragmentData())	# add base fragment
				list_obj_atom_shift[frag_i].append([])
				list_obj_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_nuc, "Base"))
				list_obj_fragments[frag_i][1].add_charge(0)

				# separate into backbone and base
				list_backbone = []
				list_base = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if "'" in obj_atom.name or obj_atom.name in ["P", "O1P", "OP1", "O2P", "OP2", "HO5", "HO3", "H5T", "H3T"]:
						# main chain
						list_backbone.append(obj_atom)

					else:
						# base
						list_base.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_backbone)
				list_obj_fragments[frag_i][1].set_atoms(list_base)

				# shift atoms in main chains
				list_remain_atoms = []
				list_shift_atoms = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if term_type in [None, "3"] and obj_atom.name in ["P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2", "1H5'", "2H5'"]:
						list_shift_atoms.append(obj_atom)

					else:
						list_remain_atoms.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_remain_atoms)
				if len(list_shift_atoms) != 0:
					if 0 < frag_i:
						list_obj_atom_shift[frag_i-1][0].extend(list_shift_atoms)
						if len([obj_atom for obj_atom in list_shift_atoms if obj_atom.name == "P"]) != 0:
							list_obj_fragments[frag_i-1][0].add_charge(-1)

					else:
						# end of shift atom
						sys.stderr.write("ERROR: Algorithm error in nucleic acid fragmentation (/base). Please report to developer.\n")
						sys.exit(1)

			elif sep_nuc == "/sugar":
				# prepare fragment
				list_obj_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
				list_obj_fragments[frag_i][0].add_charge(0)
				list_obj_fragments[frag_i].append(FragmentData())	# add sugar fragment
				list_obj_atom_shift[frag_i].append([])
				list_obj_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_nuc, "Sugar"))
				list_obj_fragments[frag_i][1].add_charge(0)
				list_obj_fragments[frag_i].append(FragmentData())	# add base fragment
				list_obj_atom_shift[frag_i].append([])
				list_obj_fragments[frag_i][2].set_type("{}{}-{}".format(res_type, sep_nuc, "Base"))
				list_obj_fragments[frag_i][2].add_charge(0)

				# separate into backbone, sugar, and base
				list_backbone = []
				list_sugar = []
				list_base = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if obj_atom.name in ["P", "O1P", "OP1", "O2P", "OP2", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2"]:
						# main chain (phosphate)
						list_backbone.append(obj_atom)

					elif "'" in obj_atom.name or obj_atom.name in ["H5T", "HO5", "H3T", "HO3"]:
						# sugar
						list_sugar.append(obj_atom)

					else:
						# base
						list_base.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_backbone)
				list_obj_fragments[frag_i][1].set_atoms(list_sugar)
				list_obj_fragments[frag_i][2].set_atoms(list_base)

				# shift atoms in main chains
				list_remain_atoms = []
				list_shift_atoms = []
				for obj_atom in list_obj_fragments[frag_i][0].atoms:
					if term_type in [None, "3"] and obj_atom.name in ["P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2", "1H5'", "2H5'"]:
						list_shift_atoms.append(obj_atom)

					else:
						list_remain_atoms.append(obj_atom)

				list_obj_fragments[frag_i][0].set_atoms(list_remain_atoms)
				if len(list_shift_atoms) != 0:
					if 0 < frag_i:
						list_obj_atom_shift[frag_i-1][1].extend(list_shift_atoms)
						if len([obj_atom for obj_atom in list_shift_atoms if obj_atom.name == "P"]) != 0:
							list_obj_fragments[frag_i-1][1].add_charge(-1)

					else:
						# end of shift atom
						sys.stderr.write("ERROR: Algorithm error in nucleic acid fragmentation (/sugar). Please report to developer.\n")
						sys.exit(1)
			continue


		if res_type == "Water":
			list_obj_fragments[frag_i][0].set_type("{}".format(res_type))
			list_obj_fragments[frag_i][0].set_atoms([obj_atom for obj_atom in obj_residue.atoms])
			list_obj_fragments[frag_i][0].add_charge(0)
			list_atom_info[obj_atom] = list_obj_fragments[frag_i][0]
			continue

		if res_type == "Ion":
			list_obj_fragments[frag_i][0].set_type("{}".format(res_type))
			list_obj_fragments[frag_i][0].set_atoms([obj_atom for obj_atom in obj_residue.atoms])
			list_atom_info[obj_atom] = list_obj_fragments[frag_i][0]
			if obj_atom.element in [3, 11, 19, 37]:
				list_obj_fragments[frag_i][0].add_charge(1)

			elif obj_atom.element in [4, 12, 20, 29, 30, 38, 56]:
				list_obj_fragments[frag_i][0].add_charge(2)

			elif obj_atom.element in [9, 17, 35, 53]:
				list_obj_fragments[frag_i][0].add_charge(-1)

			continue

		# Ligand
		list_obj_fragments[frag_i][0].set_type("{}".format("Ligand"))
		list_obj_fragments[frag_i][0].set_atoms([obj_atom for obj_atom in obj_residue.atoms])
		list_atom_info[obj_atom] = list_obj_fragments[frag_i][0]


	# update shift atoms
	for list_obj_fragments1, list_obj_atom_shift1 in zip(list_obj_fragments, list_obj_atom_shift):
		for obj_fragment, list_atom in zip(list_obj_fragments1, list_obj_atom_shift1):
			obj_fragment.set_atoms(obj_fragment.atoms + list_atom)

	# flatten fragment list
	list_obj_fragments = [v2 for v1 in list_obj_fragments for v2 in v1 if len(v2.atoms) != 0]

	# update fragment information
	list_obj_fragments = [obj_fragment.set_index(frag_idx).set_bda(0).sort_atoms() for frag_idx, obj_fragment in enumerate(list_obj_fragments, 1)]

	# update atom information
	for obj_fragment in list_obj_fragments:
		for obj_atom in obj_fragment.atoms:
			list_atom_info[obj_atom] = obj_fragment


	# make connection and bda information
	# BDA = atom in negative charged fragment; BAA = atom in positive charged fragment
	for obj_fragment in list_obj_fragments:
		list_obj_atom_member = set(obj_fragment.atoms)
		for obj_atom in obj_fragment.atoms:
			# make bonded atom list in other fragment
			list_obj_atom_partners = set(obj_atom.bond_partners) - list_obj_atom_member
			for obj_atom_partner in list_obj_atom_partners:
				obj_fragment_partner = list_atom_info[obj_atom_partner]

				# left: BDA (+1) / right: BAA (-1)
				is_bda = False
				if obj_fragment.type.endswith("Backbone") and obj_fragment_partner.type.endswith("Backbone") and obj_fragment.index < obj_fragment_partner.index:
					print("Backbone - Backbone", obj_fragment.index, obj_fragment_partner.index)
					is_bda = True

				elif obj_fragment.type.endswith("Sidechain") and obj_fragment_partner.type.endswith("Backbone"):
					is_bda = True

				elif obj_fragment.type.endswith("Backbone") and obj_fragment_partner.type.endswith("Base"):
					print("Backbone - Base", obj_fragment.index, obj_fragment_partner.index)
					is_bda = True

				elif obj_fragment.type.endswith("Backbone") and obj_fragment_partner.type.endswith("Sugar") and obj_fragment.index < obj_fragment_partner.index:
					print("Backbone - Sugar", obj_fragment.index, obj_fragment_partner.index)
					is_bda = True

				elif obj_fragment.type.endswith("Sugar") and obj_fragment_partner.type.endswith("Backbone") and obj_fragment.index < obj_fragment_partner.index:
					print("Sugar - Backbone", obj_fragment.index, obj_fragment_partner.index)
					is_bda = True

				elif obj_fragment.type.endswith("Sugar") and obj_fragment_partner.type.endswith("Base"):
					print("Sugar - Base", obj_fragment.index, obj_fragment_partner.index)
					is_bda = True

				if is_bda:
					obj_fragment_partner.append_connection([obj_atom, obj_atom_partner])
					obj_fragment_partner.add_bda(1)
					obj_fragment.add_charge(1)
					obj_fragment_partner.add_charge(-1)

	return list_obj_fragments
