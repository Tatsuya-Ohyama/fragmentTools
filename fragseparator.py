#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fragseparator
fred ファイルを元に、各フラグの PDB ファイルを生成するプログラム
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import basic_func
import argparse

# =============== functions =============== #


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "fragseparator", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-f", dest = "fred", metavar = "INPUT.fred", help = "Input fred file")
	parser.add_argument("-p", dest = "pdb", metavar = "INPUT.pdb", help = "Input pdb file")
	parser.add_argument("-o", dest = "output", metavar = "OUTPUT", help = "Prefix for output")
	args = parser.parse_args()

	basic_func.check_exist(args.fred, 2)
	basic_func.check_exist(args.pdb, 2)

	ref_atoms = []
	get_atoms = []
	flag_read = False

	# fred ファイル読み込み
	with open(args.fred, "r") as obj_input:
		re_wsp = re.compile(r"[\s\t]+")
		re_empty = re.compile(r"^[\s\t]*\n$")
		for line in obj_input:
			if "<< connections" in line:
				# 接続情報になったら終了
				break

			if re_empty.search(line):
				# 空行はスキップ
				continue

			if flag_read:
				# 読み込みフラグが True の場合
				datas = line.split("|")
				datas = re_wsp.split(datas[-1])
				tmp_atoms = list(map(lambda x : x.strip(), datas))
				ref_atoms.append(tmp_atoms)
				get_atoms.append([])
			else:
				# 1 行目はスルーする
				flag_read = True

	# PDB ファイル読み込み
	with open(args.pdb, "r") as obj_input:
		re_atom = re.compile(r"^(?:(?:ATOM)|(?:HETATM))")
		for line in obj_input:
			if re_atom.search(line):
				atom_idx = line[6:11].strip()
				idx = 0
				for refs in ref_atoms:
					if atom_idx in refs:
						get_atoms[idx].append(line)
						break
					idx += 1

	# 出力
	for idx, atoms in enumerate(get_atoms):
		output = "{prefix}_{idx:03d}.pdb".format(prefix = args.output, idx = idx + 1)
		with open(output, "w") as obj_output:
			for line in atoms:
				obj_output.write(line)
			obj_output.write("END\n")

		sys.stderr.write("INFO: {0} was created.\n".format(output))
