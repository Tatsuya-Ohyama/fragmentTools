#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

# =============== classes =============== #
class FragmentData:
	""" フラグメントデータのクラス """
	def __init__(self, line):
		line = line.strip("\r\n")
		datas = [x.strip() for x in line.split("|")]

		re_int = re.compile(r"^-?\d+$")
		self._idx = int(datas[0])
		self._charge = 0
		if re_int.search(datas[1]):
			self._charge = int(datas[1])
		else:
			self._charge = "ERR"

		self._bda = 0
		if re_int.search(datas[2]):
			self._bda = int(datas[2])
		else:
			self._bda = "ERR"

		self._atoms = sorted([int(x) for x in datas[3].split()])

	def update_fragment_index(self, index):
		""" フラグメント番号を更新するメソッド """
		self._idx = index

	def duplicate_atom(self, atoms):
		""" フラグメント構成原子の重複を削除するメソッド """
		self._atoms = sorted(list(set(self._atoms) - set(atoms)))

	def get_fragment_index(self):
		""" フラグメント番号を返すメソッド """
		return self._idx

	def get_charge(self):
		""" 電荷を返すメソッド """
		return self._charge

	def get_bda(self):
		""" BDA を返すメソッド """
		return self._bda

	def get_atoms(self):
		""" 構成原子を返すメソッド """
		return self._atoms


# =============== main =============== #
# if __name__ == '__main__':
# 	main()
