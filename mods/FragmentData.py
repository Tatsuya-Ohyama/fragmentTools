#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import parmed


# =============== classes =============== #
class FragmentData:
	""" FragmentData class """
	def __init__(self):
		# member
		self._idx = None
		self._type = None
		self._charge = None
		self._bda = None
		self._atoms = []
		self._connections = []


	@property
	def index(self):
		return self._idx

	@property
	def type(self):
		return self._type

	@property
	def charge(self):
		return self._charge

	@property
	def bda(self):
		if self._bda is not None and len(self._connections) != 0:
			if not (self._bda == len([v for v in self._connections if v[1] in self._atoms])):
				sys.stderr.write("WARNING: BDA in Fragment {0} conflicts between given BDA and connection.\n".format(self._idx))
				sys.stderr.write("         BDA: {0} / connection: {1} ({2})\n".format(self._bda, len(self._connections), ", ".join(["{0[0]}-{0[1]}".format(v) for v in self._connections])))

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


	def set_index(self, fragment_index):
		"""
		Method to set fragment index

		Args:
			index (int): fragment index
		"""
		self._idx = fragment_index
		return self


	def set_type(self, type):
		"""
		Method to set fragment type

		Args:
			type (str): fragment type
		"""
		self._type = type
		return self


	def set_charge(self, charge):
		"""
		Method to set fragment charge

		Args:
			charge (int): fragment charge
		"""
		self._charge = charge
		return self


	def add_charge(self, charge):
		"""
		Method to add fragment charge

		Args:
			charge (int): fragment charge
		"""
		if self._charge is None:
			self._charge = 0
		self._charge += charge
		return self


	def set_bda(self, bda):
		"""
		Method to set fragment BDA

		Args:
			bda (int): fragment BDA
		"""
		self._bda = bda
		return self


	def add_bda(self, bda):
		"""
		Method to add fragment BDA

		Args:
			bda (int): fragment BDA
		"""
		if self._bda is None:
			self._bda = 0
		self._bda += bda
		return self


	def append_atom(self, atom):
		"""
		Method to append fragment atom

		Args:
			atom (int): atom index
		"""
		self._atoms.append(atom)
		return self


	def set_atoms(self, list_atoms):
		"""
		Method to set fragment atoms

		Args:
			list_atoms (list): fragment atoms [atom_index(int), ...]
		"""
		self._atoms = list_atoms
		return self


	def set_connections(self, connections):
		"""
		Method to set connection list

		Args:
			connections (list): connection list
		"""
		self._connections = connections
		return self


	def append_connection(self, connection):
		"""
		Method to append/extend connection information

		Args:
			connection (list): connection information
		"""
		if type(connection) == list and type(connection[0]) == list:
			self._connections.extend(connection)

		else:
			self._connections.append(connection)

		return self


	def get_atoms(self):
		"""
		Method to output atom list

		Returns:
			list: [atom_idx(int)]
		"""
		if isinstance(self._atoms[0], parmed.topologyobjects.Atom):
			return [obj_atom.idx+1 for obj_atom in self._atoms]

		else:
			return self._atoms


	def get_charge(self):
		"""
		Methdo to output charge

		Returns:
			int: charge
		"""
		return self.charge


	def get_bda(self):
		"""
		Method to output bda

		Returns:
			int: bda
		"""
		return self._bda


	def get_connections(self):
		"""
		Method to output connection list

		Returns:
			list: connection list
		"""
		if isinstance(self._atoms[0], parmed.topologyobjects.Atom):
			return [[obj_atom.idx+1 for obj_atom in connection] for connection in self._connections]

		else:
			return [v for v in self._connections if v[1] in self._atoms]
