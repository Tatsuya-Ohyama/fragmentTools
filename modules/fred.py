#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
py_module_fred (fred の機能をモジュール化)
"""

import sys, re

# =============== variables =============== #
re_quote_h = re.compile(r"^['\"]")
re_quote_t = re.compile(r"['\"]$")
re_namelist_h = re.compile(r"^[\s\t]*\&")
re_namelist_t = re.compile(r"^[\s\t]*/$")
re_pdb_atom = re.compile(r"^((HETATM)|(ATOM))")
re_empty = re.compile(r"^[\s\t]*$")
re_wsp = re.compile(r"[\s\t]+")

# ネームリスト更新
re_natom = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
re_nf = re.compile(r"NF", re.IGNORECASE)
re_charge = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
re_autofrag = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
re_fragment = re.compile(r"\&FRAGMENT", re.IGNORECASE)
re_coord = re.compile(r"\{\.{3}\}")


# =============== functions =============== #
# convert_ajf (ajf ファイルの変換)
def convert_ajf(lists, width, n):
	lists = list(map(lambda data : str(data), lists))
	new_lists = [""]

	format_data = "%" + str(width) + "s"

	index = 0
	for i in range(len(lists)):
		new_lists[index] += format_data % lists[i]
		if i != 0 and (i + 1) % n == 0:
			index += 1
			new_lists.append("")

	if len(new_lists[len(new_lists) - 1]) == 0:
		# 指定された個数でちょうど終わる場合は、無駄な要素が追加されるので削除する
		del(new_lists[len(new_lists) - 1])

	return new_lists


# load_ajf (ajf ファイルの読み込み)
def load_ajf(file_input, file_reference):
	re_namelist_fragment_h = re.compile(r"^[\s\t]*\&FRAGMENT")
	re_namelist_fragment_t = re.compile(r"^[\s\t]*\/")

	flag_read = 0
	atom = 0
	atom_now = 0
	fragment = 0
	fragment_now = 0
	namelists = []
	fragment_atoms = []
	fragment_members = []
	fragment_members_tmp = []
	charges = []
	BDAs = []
	connections = []

	with open(file_input, "r") as obj_input:
		line_count = 0
		for line in obj_input:
			line_count += 1

			if flag_read == 0:
				# ネームリスト取得
				line = line.strip()

				if re_empty.search(line):
					# 空行はスキップ
					continue

				if re_namelist_fragment_h.search(line):
					# FRAGMENT ネームリストの始端
					flag_read = 1

				elif "ReadGeom" in line:
					if file_reference != None:
						# 参照 PDB が指定されていた場合
						basic.check_exist(file_reference, 2)
					else:
						# 参照 PDB が指定されていない場合
						file_reference = re_wsp.sub("", line)
						file_reference = file_reference.replace("ReadGeom=", "")
						file_reference = re_quote_h.sub("", file_reference)
						file_reference = re_quote_t.sub("", file_reference)

					basic.check_exist(file_reference, 2)

					with open(file_reference, "r") as obj_pdb:
						for p_line in obj_pdb:
							if re_pdb_atom.search(p_line):
								atom += 1

				if re_namelist_fragment_t.search(line):
					line = "/\n"
				elif not re_namelist_h.search(line):
					line = "  " + line
				namelists.append(line)

			elif 0 < flag_read:
				# フラグメント情報取得
				line = line.rstrip("\r\n")

				if re_namelist_fragment_t.search(line):
					# FRAGMENT ネームリストの終端
					line = "{...}\n/\n"
					flag_read = 0
					namelists.append(line)

				else:
					# 8 文字ずつ区切る
					datas = list(map(lambda data : int(data), basic.split_n(line, 8)))

					if flag_read == 1:
						# フラグメント構成原子取得
						if "0" in datas:
							# 構成原子 0 のフラグメントがある場合
							sys.stderr.write("ERROR: Invalid ajf file in %d\n" % line_count)
							sys.stderr.write("       zero atoms in fragment found\n")
							sys.exit(1)

						fragment += len(datas)
						fragment_atoms.extend(datas)

						atom_now += sum(datas)
						if atom == atom_now:
							# PDB の総原子数と現在の原子数が一致した場合
							flag_read = 2
						elif atom < atom_now:
							# PDB の総原子数以上の原子を検出した場合
							sys.stderr.write("ERROR: The number of atoms in PDB file (%s) and ajf file mismatched (%d / %d)\n" % (file_reference, atom_now, atom))
							sys.stderr.write("       Maybe wrong PDB file was specified.\n")
							sys.stderr.write("       Please fix ReadGeom in ajf file, OR use -P option.\n")
							sys.exit(1)

					elif flag_read == 2:
						# 電荷情報取得
						fragment_now += len(datas)
						charges.extend(datas)

						if fragment == fragment_now:
							# フラグメント数が総フラグメント数と一致した場合
							flag_read = 3
							fragment_now = 0
						elif fragment < fragment_now:
							# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
							sys.stderr.write("ERROR: The number of fragments mismatched in charge section (%d / %d)\n" % (fragment_now, fragment))
							sys.exit(1)

					elif flag_read == 3:
						# BDA 取得
						fragment_now += len(datas)
						BDAs.extend(datas)

						if fragment == fragment_now:
							# フラグメント数が総フラグメント数と一致した場合
							flag_read = 4
							fragment_now = 0
						elif fragment < fragment_now:
							# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
							sys.stderr.write("ERROR: The number of fragments mismatched in BDA section\n" % (fragment_now, fragment))
							sys.exit(1)

					elif flag_read == 4:
						# フラグメント構成原子の原子順序番号取得
						fragment_members_tmp.extend(datas)

						if fragment_atoms[fragment_now] <= len(fragment_members_tmp):
							fragment_now += 1
							fragment_members.append(fragment_members_tmp)
							fragment_members_tmp = []

							if fragment == fragment_now:
								# フラグメント数が総フラグメント数と一致した場合
								flag_read = 5
								now_fragment = 0
							elif fragment < fragment_now:
								# フラグメント数が総フラグメント数以上のフラグメントを検出した場合
								sys.stderr.write("ERROR: The number of fragments mismatched in fragment atom section\n" % (fragment_now, fragment))
								sys.exit(1)

					elif flag_read == 5:
						# 接続情報取得
						line = line.strip()
						connections.append(list(map(lambda data : int(data),re_wsp.split(line))))
	return fragment_atoms, charges, BDAs, fragment_members, connections, namelists


# load_fred (fred ファイルの読み込み)
def load_fred(file_input):
	flag_read = 0
	fragment = 0
	atom = 0
	flag_fragment = 0

	fragment_atoms = []
	charges = []
	BDAs = []
	fragment_members = []
	connections = []
	namelists = []

	re_connection = re.compile(r"<< connections", re.IGNORECASE)
	re_namelist_mark = re.compile(r"=+< namelist >=+")
	re_natom = re.compile(r"^[\s\t]*Natom", re.IGNORECASE)
	re_nf = re.compile(r"NF", re.IGNORECASE)
	re_charge = re.compile(r"^[\s\t]*Charge", re.IGNORECASE)
	re_autofrag = re.compile(r"^[\s\t]*AutoFrag", re.IGNORECASE)
	re_fragment = re.compile(r"\&FRAGMENT", re.IGNORECASE)
	re_coord = re.compile(r"\{\.{3}\}")

	with open(file_input, "r") as obj_input:
		line_count = 0

		for line in obj_input:
			line_count += 1
			line = line.strip()

			if line_count == 1 or re_empty.search(line):
				# 1行目と空行は無視
				continue

			elif re_connection.search(line):
				# 接続情報フラグ
				flag_read = 1

			elif re_namelist_mark.search(line):
				# この行以降がネームリスト
				flag_read = 2

				# フラグメント情報の統計
				# フラグメントに含まれる原子数算出
				for i in range(len(fragment_members)):
					fragment_atoms.append(len(fragment_members[i]))

				# フラグメント数算出
				fragment = len(fragment_members)

				# 電荷の合計算出
				charge = sum(charges)

			elif flag_read == 0:
				# フラグメント情報
				datas = line.split("|")
				datas = list(map(lambda data : data.strip(), datas))
				charges.append(int(datas[1]))
				BDAs.append(int(datas[2]))
				tmp = list(map(lambda data : int(data), re_wsp.split(datas[3])))
				fragment_members.append(tmp)
				atom += len(tmp)

			elif flag_read == 1:
				# 接続情報
				tmp = list(map(lambda data : int(data), re_wsp.split(line)))
				connections.append(tmp)

			elif flag_read == 2:
				# ネームリスト
				if re_natom.search(line):
					# 原子数の更新
					line = re.sub(r"\d+", str(atom), line)
				elif re_nf.search(line):
					# フラグメント数の更新
					line = re.sub(r"\d+", str(fragment), line)
				elif re_charge.search(line):
					# 電荷の更新
					line = re.sub(r"-?\d+", str(charge), line)
				elif re_autofrag.search(line):
					# autofrag の更新
					line = re.sub(r"=.+$", "='OFF'", line)

				if re_namelist_t.search(line):
					line = "/\n"
				elif not re_namelist_h.search(line):
					line = "  " + line

				namelists.append(line)
	return fragment_atoms, charges, BDAs, fragment_members, connections, namelists
