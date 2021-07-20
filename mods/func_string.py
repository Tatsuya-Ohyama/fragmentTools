#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
string function module
"""

import sys
import random



# =============== split_n =============== #
def split_n(line, length):
	"""
	Function to split string into n characters

	Args:
		line(str): target string
		length: number of characters to split

	Returns:
		string_list(list)
	"""
	line = line.rstrip("\r\n")
	datas = []
	pos = 0
	while pos + length <= len(line):
		datas.append(line[pos : pos + length])
		pos += length

	if pos != len(line):
		datas.append(line[pos:len(line)])

	datas = list(map(lambda data:data.strip(), datas))

	return datas


# =============== str2bool =============== #
def str2bool(string):
	"""
	Function to convert from string to bool

	Args:
		string(str): target string

	Returns:
		bool
	"""
	if type(string) == str:
		string = string.lower()
		if string == "true":
			return True
		else:
			return False
	else:
		sys.stderr.write("ERROR: function error in str2bool().\n")
		sys.exit(1)


# =============== summarized_range =============== #
def summarized_range(seq, sign_sep=",", sign_range="-"):
	"""
	Function to make int list a range specifier

	Args:
		seq(list): int list to make range specifier
		sign_sep: Symbol to separate multiple range specifiers (Default: ",")
		sign_range: Symbol indicating range (Default: "-")

	Returns:
		str
	"""
	new_seq = []

	for i in sorted(seq):
		if len(new_seq) == 0:
			# 初回登録
			new_seq.append([i])

		elif new_seq[-1][-1] + 1 != i:
			# 不連続の場合
			new_seq.append([i])

		else:
			# 連続の場合
			new_seq[-1].append(i)

	return sign_sep.join(["{0}{1}{2}".format(v[0], sign_range, v[-1]) if len(v) != 1 else str(v[0]) for v in new_seq])


# =============== target_range =============== #
def target_range(range_string, start=-1, end=-1):
	"""
	Function to convert range specifier to list

	Args:
		range_string(str): range specifier
		start(int): Start value when start value is not specified (Default: -1)
		end(int): End value when end value is not specified (Default: -1)

	Returns:
		range_list(list)
	"""
	# split by ","
	range_list = range_string.split(",")

	range_list_new = []
	for value in range_list:
		if "-" in value:
			# if range is fully specified
			split_datas = value.split("-", 2)
			if split_datas[0] == "":
				# if only end value is specified
				split_datas[0] = start
			if split_datas[1] == "":
				# if only start value is specified
				split_datas[1] = end
			range_list_new.append([int(split_datas[0]), int(split_datas[1]) + 1])
		else:
			# single specified
			if value == "":
				sys.stderr.write("ERROR: empty value in target_range() of basic.py\n")
				sys.exit(1)
			value = int(value)
			range_list_new.append([value, value + 1])

	range_list = []
	for range_value in range_list_new:
		for value in range(range_value[0], range_value[1]):
			range_list.append(value)

	return range_list


# =============== get_random_string =============== #
def get_random_string(char_string, n):
	"""
	Function to generate random string

	Args:
		char_string (str): characters compose of random string
		n (int): length of random string

	Returns:
		str: random string
	"""
	version = sys.version_info.major + sys.version_info.minor / 10
	char_string = [str(x) for x in char_string]
	if version >= 3.6:
		return "".join(random.choices(char_string, k=n))
	else:
		return "".join([random.choice(char_string) for i in range(n)])
