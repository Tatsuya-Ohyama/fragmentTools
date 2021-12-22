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
		self._bda = None
		self._atoms = []
		self._connections = []


	@property
	def fragment_index(self):
		return self._idx

	@property
	def charge(self):
		return self._charge

	@property
	def bda(self):
		if self._bda is not None and len(self._connections) != 0:
			if not (self._bda == len([v for v in self._connections if v[1] in self._atoms])):
				print("BDA: conflict (Fragment {0})".format(self._idx))

		if self._bda is not None:
			return self._bda

		elif len(self._connections) != 0:
			return len([v for v in self._connections if v[1] in self._atoms])

		else:
			return 0

	@property
	def connections(self):
		return self._connections

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


	def set_connections(self, connections):
		"""
		Method to set connection list

		Args:
			connections (list): connection list

		Returns:
			self
		"""
		self._connections = connections
		return self


	def append_connection(self, connection):
		"""
		Method to append connection information

		Args:
			connection (list): connection information

		Returns:
			self
		"""
		self._connections.append(connection)
		return self


	def get_connections(self):
		"""
		Method to output connection list

		Returns:
			list: connection list
		"""
		return [v for v in self._connections if v[1] in self._atoms]
