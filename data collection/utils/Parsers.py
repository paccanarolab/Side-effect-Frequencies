import csv
import os
import pickle
import numpy, scipy.io
from collections import defaultdict
__author__ = 'diegogaleano'
__email__  = 'Diego.Galeano.2014@rhul.live.ac.uk'
__date__  = '19-10-2016'

class EasyParsers(object):
	"""
	Broad class that parser different types of files.
	"""
	def __init__(self, directory_file=''):

		self.directory_file = directory_file

	def get_data_directory(self):
		'''
		This function set the directory_file were the data folder is.
		:param filename: name of the file in data.
		'''
		my_dir = os.getcwd()
		#code_dir = os.path.abspath(os.path.join(my_dir, os.pardir))
		#main_dir = os.path.abspath(os.path.join(code_dir, os.pardir))
		#print data_dir
		data_dir = my_dir + '/data/databases/'

		return data_dir

	def set_file_directory(self, filename):
		'''
		This function set the data directory as files directory
		:param filename:
		:return:
		'''
		self.directory_file = filename

	def get_results_directory(self):
		'''
		This function set the directory_file were the results folder is.
		:param filename: name of the file in results.
		'''
		my_dir = os.getcwd()
		code_dir = os.path.abspath(os.path.join(my_dir, os.pardir))
		main_dir = os.path.abspath(os.path.join(code_dir, os.pardir))
		result_dir = main_dir + '/data/results/'

		return result_dir
	def get_images_directory(self):
		'''
		This function set the directory_file were the results folder is.
		:param filename: name of the file in results.
		'''
		my_dir = os.getcwd()
		code_dir = os.path.abspath(os.path.join(my_dir, os.pardir))
		main_dir = os.path.abspath(os.path.join(code_dir, os.pardir))
		result_dir = main_dir + '/data/images/'

		return result_dir
		
	def parse_csv(self, filename, quote = '|'):
		'''
		Parse the csv file in plain way.
		:return: list
		'''
		data = list()
		with open(filename, 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar=quote)
			for row in spamreader:
				data.append(row)

		return data

	def parse_tsv(self, filename):
		'''
		Parse the tsv file in plain way.
		:return: list
		'''
		data = list()
		with open(filename, 'rb') as tsvfile:
			spamreader = csv.reader(tsvfile, delimiter='\t')
			for row in spamreader:
				data.append(row)
		return data

	def parse_csv_header(self, filename, column_id = 0):
		'''
		Parse the csv file.
		:param columnID: which column should be use as first key.
		:return: dictionary, keys = headers from first line.
		'''
		my_dict = defaultdict()
		header = list()
		# Plain parsing
		data = self.parse_csv(filename)

		# We get the others headers.
		for h in data[0]:
			header.append(h)

		for d in data[1::]:
			if d[column_id] not in my_dict:
				my_dict[d[column_id]] = defaultdict(list)
			
			for idx,h in enumerate(header):
				if idx != column_id:
					my_dict[d[column_id]][h].append(d[idx].replace('"',''))

		return my_dict
		
	def parse_tsv_header(self, filename, column_id = 0):
		'''
		Parse the csv file.
		:param columnID: which column should be use as first key.
		:return: dictionary, keys = headers from first line.
		'''
		my_dict = defaultdict()
		header = list()
		# Plain parsing
		data = self.parse_tsv(filename)

		# We get the others headers.
		for h in data[0]:
			header.append(h)

		for d in data[1::]:
			if d[column_id] not in my_dict:
				my_dict[d[column_id]] = defaultdict(list)
			
			for idx,h in enumerate(header):
				if idx != column_id:
					my_dict[d[column_id]][h].append(d[idx].replace('"',''))

		return my_dict

	def save_pickle(self, directory,variable_name, variable):
		'''
		Save the data in pickle format.
		:param variable_name:
		:param variable:
		:return:
		'''
		output = open(directory + variable_name, 'wb')
		pickle.dump(variable, output)
		output.close()

	def save_matlab(self, directory,variable_name, variable):
		'''
		Save the data in .mat format.
		:param variable_name:
		:param variable:
		:return:
		'''
		scipy.io.savemat(directory + variable_name,mdict ={variable_name: variable})

	def read_pickle(self,directory, variable_name):
		'''
		Read pickle files
		:param variable_name:
		:return: data.
		'''
		pkl_file = open(directory + variable_name, 'rb')
		data = pickle.load(pkl_file)
		return data

	def features_dictionary_to_npmatrix(self, my_dict, list_order, list_features, list_fcfp):
		'''
		Convert a dictionary into a np matrix according to the list that states the order of the indexes.
		Only will take the features indicated in list_features.
		:param my_dict:
		:param list_order:
		:param list_features
		:return:
		'''
		Nfeatures = len(list_features)
		Ndata = len(list_order)
		Nfcfp = len(list_fcfp)
		MatrixF1 = numpy.zeros(shape=(Ndata, Nfeatures))
		MatrixF2 = numpy.zeros(shape=(Ndata, Nfcfp))

		for idx, k in enumerate(list_order):
			for idy, f in enumerate(list_features):
				MatrixF1[idx][idy] = my_dict[k][f]

			for idy, f in enumerate(list_fcfp):
				if f in my_dict[k]['fcfp']:
					MatrixF2[idx][idy] = my_dict[k]['fcfp'][f]

		return MatrixF1, MatrixF2


	 

	#if __name__ == '__main__':
	# Example
	# csv_parser = EasyParsers(' ')
	# csv_parser.default_data_directory('BroadScreeningMoleculesFeatures.csv')
	# data = csv_parser.parse_csv_header()

	#for k,v in data.iteritems():
	#print k,v

	#print len(data)

