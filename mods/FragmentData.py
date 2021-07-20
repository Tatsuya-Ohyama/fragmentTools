#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys



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
