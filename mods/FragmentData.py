#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re



# =============== variables =============== #
RE_INT = re.compile(r"^-?\d+$")



# =============== classes =============== #
class FragmentData:
	""" FragmentData class """
	def __init__(self, line):
		# member
		self._idx = None
		self._charge = 0
		self._bda = 0
		self._atoms = []

		# initiation
		line = line.strip("\r\n")
		datas = [x.strip() for x in line.split("|")]
		self._idx = int(datas[0])

		if RE_INT.search(datas[1]):
			self._charge = int(datas[1])
		else:
			self._charge = "ERR"

		self._bda = 0
		if RE_INT.search(datas[2]):
			self._bda = int(datas[2])
		else:
			self._bda = "ERR"

		self._atoms = sorted([int(x) for x in datas[3].split()])


	@property
	def fragment_index(self):
		return self._idx

	@property
	def charge(self):
		return self._charge

	@property
	def bda(self):
		return self._bda

	@property
	def atoms(self):
		return self._atoms

	@property
	def min_index(self):
		return min(self._atoms)


	def update_fragment_index(self, index):
		"""
		Method to update fragment index

		Args:
			index (int): fragment index

		Returns:
			self
		"""
		self._idx = index
		return self


	def delete_duplicate_atom(self, atoms):
		"""
		Method to delete duplicate atoms

		Args:
			atoms (list): atom list

		Returns:
			self
		"""
		self._atoms = sorted(list(set(self._atoms) - set(atoms)))
		return self
