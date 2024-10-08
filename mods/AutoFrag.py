#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
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
	list_fragments = [[] for _ in obj_mol.residues]
	list_atom_info = {obj_atom: None for obj_atom in obj_mol.atoms}
	for frag_i, obj_residue in enumerate(obj_mol.residues):
		res_type = None
		for res_type_ref, list_resnames in RESIDUE_TYPES.items():
			if obj_residue.name in list_resnames:
				res_type = res_type_ref
				break

		if res_type == "AminoAcid":
			# separate by fragment
			if sep_amino in ["+amino", "/amino"]:
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


				# create fragment object
				if sep_amino == "+amino":
					if len(list_fragments[frag_i]) == 0:
						list_fragments[frag_i].append(FragmentData())
						list_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
					list_fragments[frag_i][0].add_charge(charge)

				elif sep_amino == "/amino":
					if len(list_fragments[frag_i]) == 0:
						list_fragments[frag_i].append(FragmentData())	# for main chain
						list_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
						list_fragments[frag_i][0].add_charge(0)

						list_fragments[frag_i].append(FragmentData())	# for side chain
						list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_amino, "Sidechain"))

					elif len(list_fragments[frag_i]) == 1:
						list_fragments[frag_i].append(FragmentData())	# for side chain
						list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_amino, "Sidechain"))

					list_fragments[frag_i][1].add_charge(charge)

				# move main chain
				for obj_atom in obj_residue.atoms:
					if obj_atom.name in ["N", "H", "CA"]:
						# main chain
						list_fragments[frag_i][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]

					elif obj_atom.name in ["C", "O"]:
						# move to next fragment
						if len(list_fragments[frag_i+1]) == 0:
							list_fragments[frag_i+1].append(FragmentData())
							list_fragments[frag_i+1][0].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
						list_fragments[frag_i+1][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i+1][0]

					elif sep_amino == "+amino":
						# side chain atom with main chain
						list_fragments[frag_i][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]

					elif sep_amino == "/amino":
						# side chain atom without main chain
						list_fragments[frag_i][1].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][1]

			# separate by residue
			if sep_amino in ["+peptide", "/peptide"]:
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


				# create fragment object
				if sep_amino == "+peptide":
					if len(list_fragments[frag_i]) == 0:
						list_fragments[frag_i].append(FragmentData())
						list_fragments[frag_i].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
					list_fragments[frag_i][0].add_charge(charge)

				elif sep_amino == "/peptide":
					if len(list_fragments[frag_i]) == 0:
						list_fragments[frag_i].append(FragmentData())	# for main chain
						list_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_amino, "Backbone"))
						list_fragments[frag_i][0].add_charge(0)

						list_fragments[frag_i].append(FragmentData())	# for side chain
						list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_amino, "Sidechain"))
						list_fragments[frag_i][1].add_charge(charge)

					elif len(list_fragments[frag_i]) == 1:
						list_fragments[frag_i].append(FragmentData())	# for side chain
						list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_amino, "Sidechain"))

					list_fragments[frag_i][1].add_charge(charge)

				for obj_atom in obj_residue.atoms:
					if obj_atom.name in ["C", "O", "N", "H", "CA"]:
						# main chain
						list_fragments[frag_i][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]

					elif sep_amino == "+amino":
						# side chain atom with main chain
						list_fragments[frag_i][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]

					elif sep_amino == "/amino":
						list_fragments[frag_i][1].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][1]

			continue

		if res_type == "NucleicAcid":
			# determine charge
			if sep_nuc == "+base":
				list_fragments[frag_i].append(FragmentData())
				list_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
				list_fragments[frag_i][0].add_charge(0)

			elif sep_nuc == "/base":
				list_fragments[frag_i].append(FragmentData())
				list_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
				list_fragments[frag_i][0].add_charge(0)
				list_fragments[frag_i].append(FragmentData())
				list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_nuc, "Base"))
				list_fragments[frag_i][1].add_charge(0)

			elif sep_nuc == "/sugar":
				if len(list_fragments[frag_i]) == 0:
					list_fragments[frag_i].append(FragmentData())
					list_fragments[frag_i][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
					list_fragments[frag_i][0].add_charge(0)
					list_fragments[frag_i].append(FragmentData())
					list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_nuc, "Sugar"))
					list_fragments[frag_i][1].add_charge(0)
					list_fragments[frag_i].append(FragmentData())
					list_fragments[frag_i][2].set_type("{}{}-{}".format(res_type, sep_nuc, "Base"))
					list_fragments[frag_i][2].add_charge(0)

				elif len(list_fragments[frag_i]) == 1:
					list_fragments[frag_i].append(FragmentData())
					list_fragments[frag_i][1].set_type("{}{}-{}".format(res_type, sep_nuc, "Sugar"))
					list_fragments[frag_i][1].add_charge(0)
					list_fragments[frag_i].append(FragmentData())
					list_fragments[frag_i][2].set_type("{}{}-{}".format(res_type, sep_nuc, "Base"))
					list_fragments[frag_i][2].add_charge(0)

			for obj_atom in obj_residue.atoms:
				if sep_nuc == "+base":
					# nucleotide
					if obj_atom.name in ["P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2"]:
						# move to next
						if frag_i == 0:
							list_fragments[frag_i][0].append_atom(obj_atom)
							if obj_atom.name == "P":
								list_fragments[frag_i][0].add_charge(-1)
							list_atom_info[obj_atom] = list_fragments[frag_i][0]

						else:
							list_fragments[frag_i-1][0].append_atom(obj_atom)
							if obj_atom.name == "P":
								list_fragments[frag_i-1][0].add_charge(-1)
							list_atom_info[obj_atom] = list_fragments[frag_i-1][0]

					else:
						list_fragments[frag_i][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]


				elif sep_nuc == "/base":
					if obj_atom.name in ["P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2"]:
						# move to next
						if frag_i == 0:
							list_fragments[frag_i][0].append_atom(obj_atom)
							if obj_atom.name == "P":
								list_fragments[frag_i][0].add_charge(-1)
							list_atom_info[obj_atom] = list_fragments[frag_i][0]

						else:
							list_fragments[frag_i-1][0].append_atom(obj_atom)
							if obj_atom.name == "P":
								list_fragments[frag_i-1][0].add_charge(-1)
							list_atom_info[obj_atom] = list_fragments[frag_i-1][0]


					elif "'" in obj_atom.name or obj_atom.name in ["H5T", "H3T"]:
						# backbone
						list_fragments[frag_i][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]

					else:
						# base
						list_fragments[frag_i][1].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][1]

				elif sep_nuc == "/sugar":
					if not obj_atom.residue.name.endswith("5") and obj_atom.name in ["H5T", "P", "O1P", "OP1", "O2P", "OP2", "O5'", "C5'", "H5'", "H5'1", "H5''", "H5'2"]:
						# main chain (phosphate)
						list_fragments[frag_i][0].append_atom(obj_atom)
						if obj_atom.name == "P":
							list_fragments[frag_i][0].add_charge(-1)
						list_atom_info[obj_atom] = list_fragments[frag_i][0]

					elif obj_atom.name == "O3'" and not obj_atom.residue.name.endswith("3"):
						if len(list_fragments[frag_i+1]) == 0:
							list_fragments[frag_i+1].append(FragmentData())
							list_fragments[frag_i+1][0].set_type("{}{}-{}".format(res_type, sep_nuc, "Backbone"))
						list_fragments[frag_i+1][0].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i+1][0]

					elif "'" in obj_atom.name or obj_atom.name in ["H3T", "H5T"]:
						# sugar
						list_fragments[frag_i][1].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][1]

					else:
						# base
						list_fragments[frag_i][2].append_atom(obj_atom)
						list_atom_info[obj_atom] = list_fragments[frag_i][2]

			continue

		if res_type == "Water":
			list_fragments[frag_i].append(FragmentData())
			list_fragments[frag_i][0].set_type("{}".format(res_type))
			list_fragments[frag_i][0].set_atoms([obj_atom for obj_atom in obj_residue.atoms])
			list_fragments[frag_i][0].add_charge(0)
			list_atom_info[obj_atom] = list_fragments[frag_i][0]
			continue

		if res_type == "Ion":
			list_fragments[frag_i].append(FragmentData())
			list_fragments[frag_i][0].set_type("{}".format(res_type))
			list_fragments[frag_i][0].set_atoms([obj_atom for obj_atom in obj_residue.atoms])
			list_atom_info[obj_atom] = list_fragments[frag_i][0]
			if obj_atom.element in [3, 11, 19, 37]:
				list_fragments[frag_i][0].add_charge(1)

			elif obj_atom.element in [4, 12, 20, 29, 30, 38, 56]:
				list_fragments[frag_i][0].add_charge(2)

			elif obj_atom.element in [9, 17, 35, 53]:
				list_fragments[frag_i][0].add_charge(-1)

			continue

		# Ligand
		list_fragments[frag_i].append(FragmentData())
		list_fragments[frag_i][0].set_type("{}".format("Ligand"))
		list_fragments[frag_i][0].set_atoms([obj_atom for obj_atom in obj_residue.atoms])
		list_atom_info[obj_atom] = list_fragments[frag_i][0]

	# flatten fragment list
	list_fragments = [v2 for v1 in list_fragments for v2 in v1 if len(v2.atoms) != 0]
	list_fragments = [obj_fragment.set_index(frag_idx) for frag_idx, obj_fragment in enumerate(list_fragments, 1)]

	# make connection and bda information
	for obj_fragment in list_fragments:
		list_connections = []
		list_obj_atom_member = set(obj_fragment.atoms)
		for obj_atom in obj_fragment.atoms:
			list_obj_atom_partners = set(obj_atom.bond_partners)
			list_obj_atom_other = list_obj_atom_partners - list_obj_atom_member
			list_connections.extend([[obj_atom_partner, obj_atom] for obj_atom_partner in list_obj_atom_other])

		n_connection = len(list_connections)
		if obj_fragment.type in ["AminoAcid+amino-Backbone", "AminoAcid+peptide-Backbone", "NucleicAcid+base-Backbone", "NucleicAcid/sugar-Backbone"]:
			# chain type
			if n_connection == 1:
				# terminal
				if list_connections[0][0].idx > list_connections[0][1].idx:
					# N-terminal, 5'-terminal, or initial fragment
					obj_fragment.set_connections([])
					obj_fragment.set_bda(0)

				else:
					# C-terminal, 3'-terminal, or end fragment
					obj_fragment.set_connections(list_connections)
					obj_fragment.set_bda(1)
					obj_fragment.add_charge(-1)
					obj_fragment_partner = list_atom_info[list_connections[0][0]]
					obj_fragment_partner.add_charge(1)

			elif n_connection == 2:
				for obj_atom_partner, obj_atom in list_connections:
					if obj_atom.idx < obj_atom_partner.idx:
						continue

					obj_fragment.set_connections([[obj_atom_partner, obj_atom]])
					obj_fragment.add_charge(-1)
					obj_fragment_partner = list_atom_info[obj_atom_partner]
					obj_fragment_partner.add_charge(1)
				obj_fragment.set_bda(1)

			else:
				sys.stderr.write("ERROR: undefined bda type in fragment #{}.\n".format(obj_fragment.index))
				sys.exit(1)

		elif obj_fragment.type in ["AminoAcid/amino-Backbone", "AminoAcid/peptide-Backbone", "NucleicAcid/base-Backbone", "NucleicAcid/sugar-Sugar"]:
			# chain with branch type
			if n_connection == 2:
				# terminal
				for obj_atom_partner, obj_atom in list_connections:
					if not (list_atom_info[obj_atom_partner].type.endswith("Backbone") or list_atom_info[obj_atom_partner].type.endswith("Sugar")):
						# skip sidechain
						continue

					# backbone
					if obj_atom.idx < obj_atom_partner.idx:
						# N-terminal, 5'-terminal, or initial fragment
						obj_fragment.set_connections([])
						obj_fragment.set_bda(0)

					else:
						# C-terminal, 3'-terminal, or end fragment
						obj_fragment.set_connections([[obj_atom_partner, obj_atom]])
						obj_fragment.set_bda(1)
						obj_fragment.add_charge(-1)
						obj_fragment_partner = list_atom_info[obj_atom_partner]
						obj_fragment_partner.add_charge(1)

			elif n_connection == 3:
				for obj_atom_partner, obj_atom in list_connections:
					obj_fragment_partner = list_atom_info[obj_atom_partner]
					if obj_fragment_partner.type.endswith("Backbone") and obj_atom.idx > obj_atom_partner.idx:
						obj_fragment.append_connection([obj_atom_partner, obj_atom])
						obj_fragment.add_charge(-1)
						obj_fragment_partner.add_charge(1)
				obj_fragment.set_bda(1)

			else:
				sys.stderr.write("ERROR: undefined bda type in fragment #{}.\n".format(obj_fragment.index))
				sys.exit(1)


		elif obj_fragment.type in ["AminoAcid/amino-Sidechain", "AminoAcid/peptide-Sidechain", "NucleicAcid/base-Base", "NucleicAcid/sugar-Base"]:
			# branch type
			if n_connection == 1:
				obj_fragment.set_connections(list_connections)
				obj_fragment.set_bda(1)
				obj_fragment.add_charge(-1)
				obj_fragment_partner = list_atom_info[list_connections[0][0]]
				obj_fragment_partner.add_charge(1)

			else:
				sys.stderr.write("ERROR: undefined bda type in fragment #{}.\n".format(obj_fragment.index))
				sys.exit(1)

		elif obj_fragment.type in ["Water", "Ion"]:
			obj_fragment.set_connections([])
			obj_fragment.set_bda(0)

		elif obj_fragment.type == "Ligand":
			sys.stderr.write("WARNING: Unable to determine charge and BDA of fragment #{}, due to the unknown residue.\n".format(obj_fragment.index))

		else:
			sys.stderr.write("ERROR: undefined bda type in fragment #{}.\n".format(obj_fragment.index))

			sys.exit(1)

	return list_fragments
