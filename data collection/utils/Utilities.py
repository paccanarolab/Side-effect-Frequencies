import operator
import numpy as np

__author__ = 'diegogaleano'
__email__  = 'Diego.Galeano.2014@rhul.live.ac.uk'
__date__  = '19-10-2016'

class MyUtilities(object):
	"""
	Implementation of sub-structural analysis algorithms.
	"""

	def __init__(self, directory_file=''):

		self.directory_file = directory_file
		
	def DoubleDicttoMatrix(self, mydict):
		'''
			Create a numpy matrix.
				rows: key1 of the dictionary.
				columns: key2 of the dictionary.
				
			returns:
				matrix, rows_keys, column_keys
		'''
		
		rows = list(mydict.keys())
		columns = set()
		for r in rows: 
			for v in mydict[r]:
				columns.add(v)
		
		columns = list(columns)
		
		N1 = len(rows)
		N2 = len(columns)		
		M = np.zeros(shape=(N1,N2))
		
		for idx1, r in enumerate(rows):							
			for idx2,c in enumerate(columns):				
				if c in mydict[r]:					
					M[idx1][idx2] = mydict[r][c][0]
				
		return M, rows, columns
		
	def DoubleDicttoMatrixConstrained(self, mydict, my_row_list, my_col_list):
		'''
			Create a numpy matrix.
				rows: key1 of the dictionary.
				columns: key2 of the dictionary.
				
			returns:
				matrix, rows_keys, column_keys
		'''
		
		N1 = len(my_row_list)
		N2 = len(my_col_list)		
		M = np.zeros(shape=(N1,N2))
		
		for idx1, r in enumerate(my_row_list):							
			for idx2,c in enumerate(my_col_list):				
				if r in mydict and c in mydict[r]:					
					M[idx1][idx2] = mydict[r][c][0]
				
		return M
		
	def DoubleDicttoMatrixConstrainedValue(self, mydict, my_row_list, my_col_list, value= 1):
		'''
			Create a numpy matrix.
				rows: key1 of the dictionary.
				columns: key2 of the dictionary.
				
			returns:
				matrix, rows_keys, column_keys
		'''
		
		N1 = len(my_row_list)
		N2 = len(my_col_list)		
		M = np.zeros(shape=(N1,N2))
		
		for idx1, r in enumerate(my_row_list):							
			for idx2,c in enumerate(my_col_list):				
				if r in mydict and c in mydict[r]:					
					M[idx1][idx2] = value
				
		return M