#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re



# =============== variables =============== #
RE_INT = re.compile(r"^-?\d+$")



# =============== classes =============== #
class FragmentData:
	""" FragmentData class """
	def __init__(self):
		# member
		self._idx = None
		self._charge = 0
		self._bda = 0
		self._atoms = []


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
		print(self._atoms)
		return min(self._atoms)


	def set_fragment_index(self, fragment_index):
		"""
		Method to set fragment index

		Args:
			fragment_index (int): fragment index

		Returns:
			self
		"""
		self._idx = fragment_index
		return self


	def set_charge(self, charge):
		"""
		Method to set fragment charge

		Args:
			charge (int): fragment charge

		Returns:
			self
		"""
		self._charge = charge
		return self


	def set_bda(self, bda):
		"""
		Method to set fragment BDA

		Args:
			bda (int): fragment BDA

		Returns:
			self
		"""
		self._bda = bda
		return self


	def set_atoms(self, atoms):
		"""
		Method to set fragment atoms

		Args:
			atoms (list): fragment atoms

		Returns:
			self
		"""
		self._atoms = atoms
		return self


	def create_from_string(self, line_val):
		"""
		Method to create fragment object from string

		Args:
			line_val (str): string information (.fred format)

		Returns:
			self
		"""
		line_val = line_val.strip()
		datas = [v.strip() for v in line_val.split("|")]
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
		return self
