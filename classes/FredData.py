#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
from classes.FragmentData import FragmentData

# =============== variables =============== #
re_fragment = re.compile(r"^[\s\t]*\d+[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|[\s\t]*(?:(?:ERR)|(?:-?\d+))[\s\t]*|(?:[\s\t]*\d+)+")
re_connection = re.compile(r"^(?:[\s\t]*\d+){2}[\s\t]*$")


# =============== classes =============== #
class FredData:
	""" fred ファイルのクラス """
	def __init__(self, input_file):
		self._charge = 0
		self._atom = 0
		self._fragments = []
		self._connection = []
		self._others = []

		self._load_file(input_file)


	def _load_file(self, input_file):
		""" fred データを読み込むメソッド """
		with open(input_file, "r") as obj_input:
			cnt_line = 0
			flag_read = 0
			for line in obj_input:
				cnt_line += 1
				line = line.rstrip("\r\n")

				if cnt_line == 2:
					flag_read = 1

				elif "connections" in line:
					flag_read = 2

				elif "namelist" in line:
					flag_read = 3

				if flag_read == 1 and re_fragment.search(line):
					fragment = FragmentData(line)
					if fragment.get_charge() != "ERR":
						self._charge += fragment.get_charge()
					self._fragments.append(fragment)
					self._atom += len(fragment.get_atoms())

				elif flag_read == 2 and re_connection.search(line):
					self._connection.append([int(x) for x in line.strip().split()])

				elif flag_read == 3:
					self._others.append(line)


	def add_fragment(self, obj_fragment):
		""" フラグメントを追加するメソッド """
		for idx, fragment in enumerate(self._fragments):
			fragment.duplicate_atom(obj_fragment.get_atoms())
			if len(fragment.get_atoms()) == 0:
				self._fragments.pop(idx)

		self._atom += len(obj_fragment.get_atoms())
		self._fragments.append(obj_fragment)
		return self


	def add_connection(self, connection):
		""" フラグメントの接続情報を増やすメソッド """
		if isinstance(connection, int):
			for i in range(connection):
				self._connection.append(["*", "*"])
		elif isinstance(connection, list):
			self._connection.append(connection)
		return self


	def get_fragment_number(self):
		""" フラグメント数を返すメソッド """
		return len(self._fragments)


	def get_last_fragment_index(self):
		""" 最後のフラグメント番号を返すメソッド """
		return self._fragments[-1].get_fragment_index()


	def write_file(self, output_file):
		""" ファイルに出力するメソッド """
		re_NF = re.compile(r"NF=[\s\t]*-?\d+")
		tmp_fragments = sorted([[obj_fragment.get_min_idx(), obj_fragment] for obj_fragment in self._fragments], key = lambda x : x[0])
		self._fragments = [obj_fragment[1].update_fragment_index(idx + 1) for idx, obj_fragment in enumerate(tmp_fragments)]
		self._connection = sorted(self._connection, key = lambda x : x[0])

		with open(output_file, "w") as obj_output:
			obj_output.write("  FNo.  | Charge | BDA | Atoms of fragment\n")
			for idx, fragment in enumerate(self._fragments):
				fragment.update_fragment_index(idx + 1)
				obj_output.write("{0:>7} |{1:>6}  |{2:>3}  |{3}\n".format(fragment.get_fragment_index(), fragment.get_charge(), fragment.get_bda(), " ".join(["{0:>8}".format(x) for x in fragment.get_atoms()])))
			obj_output.write("\n")

			obj_output.write("<< connections (ex. \"Next_fragment_atom   Prev_fragment_atom\") >>\n")
			for connection in self._connection:
				obj_output.write("{0}\n".format(" ".join(["{0:>9}".format(x) for x in connection])))
			obj_output.write("\n")

			for line in self._others:
				redb_NF = re_NF.search(line)
				if redb_NF:
					line = line[: redb_NF.start()] + "NF={0}".format(len(self._fragments)) + line[redb_NF.end() :]
				obj_output.write(line + "\n")


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
