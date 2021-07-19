#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)



# =============== constant =============== #
INDENT = "  "



# =============== class =============== #
class FileAJF:
	""" .ajf file class """
	def __init__(self, ajf_file):
		# member
		self._parameters = {}

		# init
		self.read_ajf(ajf_file)


	@property
	def parameters(self):
		return self._parameters


	def read_ajf(self, input_file):
		"""
		function to read .ajf file

		Args:
			ajf_file (str): .ajf file (input)

		Returns:
			self
		"""
		self._parameters = {"LIST_ORDER": []}
		is_list_type = False
		with open(input_file, "r") as obj_input:
			group_name = None
			for line_val in obj_input:
				line_orig = line_val
				line_val = line_val.strip()
				if len(line_val.strip()) == 0:
					continue

				if line_val.startswith("&"):
					# start of name list
					group_name = line_val
					self._parameters["LIST_ORDER"].append(group_name)

					if group_name in ["&XYZ", "&FRAGMENT", "&FRAGPAIR"]:
						is_list_type = True
						self._parameters[group_name] = []
					else:
						is_list_type = False
						self._parameters[group_name] = {}

				elif line_val.startswith("/"):
					# end of name list
					group_name = None
					is_list_type = False
					continue

				elif is_list_type:
					# in name list of list type (order is important)
					self._parameters[group_name].append(line_orig)

				else:
					# in namelist of non-list type (order is not important)
					tmp_values = [v.strip() for v in line_val.split("=", maxsplit = 1)]
					self._parameters[group_name][tmp_values[0]] = tmp_values[1]

		return self


	def set_parameters(self, parameters):
		"""
		Method to set parameters

		Args:
			parameters (dict): parameters

		Returns:
			self
		"""
		self._parameters = parameters
		return self


	def write_ajf(self, output_file):
		"""
		Method to write .ajf file

		Args:
			output_file (str): .ajf file

		Returns:
			self
		"""
		with open(output_file, "w") as obj_output:
			for group_name in self._parameters["LIST_ORDER"]:
				obj_output.write("{0}\n".format(group_name))
				if isinstance(self._parameters[group_name], dict):
					for parameter_name, parameter_value in self._parameters[group_name].items():
						obj_output.write("{0}{1}={2}\n".format(
							INDENT,
							parameter_name,
							parameter_value
						))
				else:
					for line_val in self._parameters[group_name]:
						obj_output.write(line_val)
				obj_output.write("/\n\n")
