import csv
import os,sys
# Add utils to the path
sys.path.insert(0, os.getcwd() + '/utils/')

import Parsers as parsers
import pickle
import numpy, scipy.io
from xml.etree.ElementTree import iterparse
from collections import defaultdict
from pprint import pprint
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles
from collections import Counter

__author__ = 'diegogaleano'
__email__  = 'Diego.Galeano.2014@rhul.live.ac.uk'
__date__  = '04-03-2017'

class SIDERParser(object):
	'''
	 Here we parser the files from SIDER: database for side-effects
	'''
	def __init__(self):
		
		# Get the default directories
		self.easyparser = parsers.EasyParsers()
		self.data_directory = self.easyparser.get_data_directory()
		self.result_directory = self.easyparser.get_results_directory()
	
	def parser_meddra_all_se(self, my_preferred_term = 'PT'):
		'''
		Description of the file in the README:
		1 & 2: STITCH compound ids (flat/stereo, see above)
		3: UMLS concept id as it was found on the label
		4: MedDRA concept type (LLT = lowest level term, PT = preferred term; in a few cases the term is neither LLT nor PT)
		5: UMLS concept id for MedDRA term
		6: side effect name

		All side effects found on the labels are given as LLT. Additionally, the PT is shown. There is at least one
		PT for every LLT, but sometimes the PT is the same as the LLT. LLTs are sometimes too detailed, and therefore
		you might want to filter for PT. E.g. for this term:

		PT      C0235431        Blood creatinine increased

		there are several LLT (leftmost number = count in SIDER 4.1)

		149     C0151578        LLT     C0151578        Creatinine increased
		100     C0235431        LLT     C0235431        Blood creatinine increased
		93      C0700225        LLT     C0700225        Serum creatinine increased
		2       C0858118        LLT     C0858118        Plasma creatinine increased

		All of these LLTs are equivalent for most purposes and to the same PT. 

		344     PT      C0235431        Blood creatinine increased

		The mapping was performed by extracting the LLT-->PT relations from UMLS. 

		'''
		
		
		# get all the side-effects
		selist = self.easyparser.parse_tsv(self.data_directory + 'meddra_all_se.tsv')
		
		# get the stitch CID (PubChem ID)
		drugs_se = defaultdict(set)
		
		for row in selist:
		    # the stereo ID is equivalent to PubChem ID
			drug_id = int(row[1][3::])
			
			# term type: LLT or PT
			termtype = row[3]		
			
			if my_preferred_term == termtype:
				# side-effect
				se = row[5].lower().strip()
				
				drugs_se[str(drug_id)].add(se)
		
		return drugs_se
		
	def parser_meddra_freq(self, my_preferred_term = 'PT'):
		'''
		Description of the file in the README:
		
		This file contains the frequencies of side effects as extracted from the labels. Format:

		1 & 2: STITCH compound ids (flat/stereo, see above)
		3: UMLS concept id as it was found on the label
		4: "placebo" if the info comes from placebo administration, "" otherwise
		5: a description of the frequency: for example "postmarketing", "rare", "infrequent", "frequent", "common", or an exact
		   percentage
		6: a lower bound on the frequency
		7: an upper bound on the frequency
		8-10: MedDRA information as for meddra_all_se.tsv.gz

		The bounds are ranges like 0.01 to 1 for "frequent". If the exact frequency is known, then the lower bound
		matches the upper bound. Due to the nature of the data, there can be more than one frequency for the same label,
		e.g. from different clinical trials or for different levels of severeness.

		'''
		
		# get all the side-effects
		selist = self.easyparser.parse_tsv(self.data_directory + 'meddra_freq.tsv')
		
		# get the stitch CID (PubChem ID)
		drugs_se = defaultdict()
		
		for row in selist:
		    # the stereo ID is equivalent to PubChem ID
			drug_id = int(row[1][3::])
			
			# placebo, otherwise ''
			placebo = False
			if 'placebo' in row:
				placebo = True
			
			# frequency of the side-effect 
			frequency = row[4]

			# term type: LLT or PT
			termtype = row[7]

			# side-effect
			se = row[9].lower().strip()
			
			if my_preferred_term == termtype:
				if drug_id not in drugs_se:
					drugs_se[drug_id] = dict()
				
				if se not in drugs_se[drug_id]:
					drugs_se[drug_id][se] = defaultdict(list)
					
				drugs_se[drug_id][se]['placebo'].append(placebo)
				drugs_se[drug_id][se]['frequency'].append(frequency)
				
		return drugs_se
		
	def __flattoStereoID(self):
		'''
			dict:
				key: flat id, value: stereo id (PubChemID)
		'''
		# get all the side-effects
		selist = self.easyparser.parse_tsv(self.data_directory + 'meddra_all_se.tsv')
		
		# get the stitch CID (PubChem ID)
		drugID = dict()
		
		for row in selist:
		    # the stereo ID is equivalent to PubChem ID
			stereo_id = int(row[1][3::])
			
			
			flat_id = int(row[0][4::])
			
			drugID[flat_id] = stereo_id
		
		return drugID
		
	def parser_all_indications(self, my_preferred_term = 'PT'):
		'''
			0: STITCH compound id (flat, see above)
			1: UMLS concept id as it was found on the label
			2: method of detection: NLP_indication / NLP_precondition / text_mention
			3: concept name
			4: MedDRA concept type (LLT = lowest level term, PT = preferred term; in a few cases the term is neither LLT nor PT)
			5: UMLS concept id for MedDRA term
			6: MedDRA concept name

			All side effects found on the labels are given as LLT. Additionally, the PT is shown. There is at least one
			PT for every LLT, but sometimes the PT is the same as the LLT.
		'''
		# get all the side-effects
		indicationlist = self.easyparser.parse_tsv(self.data_directory + 'meddra_all_indications.tsv')
		
		# get the stitch CID 
		drugs_indication = defaultdict(set)
		
		# get flat to stereo ID
		IDConv = self.__flattoStereoID()
		
		for row in indicationlist:
		    # flat drug ID 
			flat_id = int(row[0][4::])
			
			# stereo ID
			if flat_id in IDConv:
				stereo_id = IDConv[flat_id]
				
				#indication
				myindication = row[6].lower().strip()
				
				# term type: LLT or PT
				termtype = row[4]
				
				# only if it was text mention
				if my_preferred_term == termtype:
					if 'text_mention' in row[2].strip():
						drugs_indication[str(stereo_id)].add(myindication)
				
		return drugs_indication
	
	def se_freq_breakdown(self, drugs_se):
		'''
		   REQUIREMENT: input should be output of parser_meddra_freq(...)
		   
		   We create a dictionary of drug-se pairs and we divide the info into:
		   1- exact_freq: exact frequency of the side effect. These are percentages, i.e., 5%.
		   2- placebo_exact_freq: frequency for the placebo (if any). There are always percentages,i.e. 3%
		   3- placebo_range_freq: frequency for the placebo (if any). There are always percentages,i.e. 3%
		   4- range_freq: frequency range. Sometimes it is provide frequency range 1-5%, 1 to 5% 
		   5- label_freq: label (common, very common, rare, ...).
		'''

		drug_se_pair = dict()
		   
		for drug_id, side_effects in drugs_se.iteritems():    
			drugID = str(drug_id) + "|"
			
			for se, data in side_effects.iteritems():
			    # the pair ID
				
				pairID = drugID + se		
				
				drug_se_pair[pairID] =  defaultdict(list)          
			 
				isPlacebo = data['placebo']
				
				# For this particular drug and side-effect pair:
				for Idx, fq in enumerate(data['frequency']):      
				    
					# CASE: Frequency is not placebo
					if isPlacebo[Idx] == False: 
					
						# subCASE: Frequency is an exact number
						try:
							exact_freq = float(fq.strip('%'))
							drug_se_pair[pairID]['exact_freq'].append(exact_freq)

						except:
							
							
							# Fix for cases in which we have FREQUENCIES like 1-9%. We append the average.                      
							fq = fq.replace('%','').replace(' ','').replace('-',';').replace('to',';').replace('<','') 

							# subCASE: is a range of frequency
							if ';' in fq: # it is a range 							   
								n1 = float(fq.split(';')[0])
								n2 = float(fq.split(';')[1])
								
								# We append both data in the range.
								drug_se_pair[pairID]['range_freq'].append(n1)
								drug_se_pair[pairID]['range_freq'].append(n2)

							else: # subCASE: is a label frequency
								fq = self.__NormalizeFrequencyLabels(fq) # normalization.
								drug_se_pair[pairID]['label_freq'].append(fq)

					else: # placebo frequencies processing
						
						# subCASE: Frequency is an exact number
						try:
							exact_freq = float(fq.strip('%'))
							drug_se_pair[pairID]['placebo_exact_freq'].append(exact_freq)
						except:
							# Fix for cases in which we have FREQUENCIES like 1-9%. We append the average.
							fq = fq.replace('%','').replace(' ','').replace('-',';').replace('to',';').replace('<','')                   

							if ';' in fq: # placebo is a range.
								n1 = float(fq.split(';')[0])
								n2 = float(fq.split(';')[1])
								drug_se_pair[pairID]['placebo_range_freq'].append(n1)
								drug_se_pair[pairID]['placebo_range_freq'].append(n2)

		return drug_se_pair  
		
	def VennCounterFreqType(self, drug_se_pair):
		'''
		 REQUIREMENT: input should be output of drug_se_pair(...)
		This function provides statistics about the type of data in the side-effects frequencies.
		   REQUIREMENT: input should be output of se_freq_breakdown(...)
		   This function counts all the possible intersection between the sets A, B and C.
		   A = ['exact_freq', 'range_freq']
		   B = ['label_freq']
		   C = ['placebo_exact_freq', 'placebo_range_freq']
		   
		   The function returns a dictionary.
		'''
		# Counter to be return
		set_counter = dict()

		# useful definitions
		typeFrequencies = ['exact_freq', 'range_freq', 'label_freq', 'placebo_exact_freq', 'placebo_range_freq']		
		A = ['exact_freq', 'range_freq']
		B = ['label_freq']
		C = ['placebo_exact_freq', 'placebo_range_freq']
		sets = ['A', 'B', 'C']
		operations = ['A-all','B-all', 'C-all', 'A & B', 'A & C', 'B & C', 'A & B & C']
		
		# possible set operations
		for v in operations:
			set_counter[v] = 0

		for eachpair, fiels in drug_se_pair.iteritems():
		   			
			if 'exact_freq' in fiels or 'range_freq' in fiels: # set A
				ban = 0 
			   
				if 'label_freq' in fiels: # set B
					set_counter['A & B'] += 1             
					
					ban = 1     
					if 'placebo_exact_freq' in fiels or 'placebo_range_freq' in fiels: # set C
						set_counter['A & B & C'] += 1 


				if 'placebo_exact_freq' in fiels or 'placebo_range_freq' in fiels: # set C
					set_counter['A & C'] += 1 
					ban = 1


				if ban == 0:
					set_counter['A-all'] += 1 

			if 'label_freq' in fiels: # set B.
			   
				ban = 0  
				if 'placebo_exact_freq' in fiels or 'placebo_range_freq' in fiels: # set C
					ban = 1
					set_counter['B & C'] += 1 					

				if 'exact_freq' in fiels or 'range_freq' in fiels: # set A
					ban = 1

				if ban == 0:
					set_counter['B-all'] += 1

			if 'placebo_exact_freq' in fiels or 'placebo_range_freq' in fiels: # set C
				ban = 0
			   
				if 'exact_freq' in fiels or 'range_freq' in fiels or 'label_freq' in fiels: # set A or B
					ban = 1
				#else:
					#print drugSEpair

				if ban == 0:
					set_counter['C-all'] += 1 

		return set_counter
	
	def plotVennDiagramFreqType(self, set_counter, directory = ''):
		'''	
		 REQUIREMENT: input should be output of VennCounterFreqType(...)
		'''
		#%matplotlib inline
		fig = plt.figure(figsize=(10,15))
		
		
		v = venn3(subsets=(set_counter['A-all'], set_counter['B-all'], set_counter['A & B'], set_counter['C-all'],
		set_counter['A & C'], set_counter['B & C'], set_counter['A & B & C']),
		set_labels = ('Frequency values (exact or range) \n (Set A) ', 'Frequency label\n (Set B)', 'Placebo frequency \n (Set C)'))
		#v.get_patch_by_id('100').set_alpha(1.0)

		# Change the color of the parts
		#v.get_patch_by_id('100').set_color('gray')
		#v.get_patch_by_id('110').set_color('gray')
		#v.get_patch_by_id('101').set_color('gray')
		#v.get_patch_by_id('111').set_color('gray')
		#v.get_patch_by_id('010').set_color('gray')
		#v.get_patch_by_id('011').set_color('gray')
		#v.get_patch_by_id('001').set_color('gray')

		# Replace the numbers by text in each part.
		v.get_label_by_id('100').set_text('A- all\n' + str(set_counter['A-all']))
		v.get_label_by_id('110').set_text('A & B\n' + str(set_counter['A & B'] ))
		v.get_label_by_id('101').set_text('A & C\n' + str(set_counter['A & C']))
		v.get_label_by_id('111').set_text('A & B & C\n'+ str(set_counter['A & B & C']))
		v.get_label_by_id('010').set_text('B-all\n' + str(set_counter['B-all'] ))
		v.get_label_by_id('011').set_text('B & C\n' + str(set_counter['B & C']))
		v.get_label_by_id('001').set_text('C-all\n' + str(set_counter['C-all']))

		# Percentage of the total number of pairs
		# v.get_label_by_id('100').set_text('10823\n 18.8 %')
		# v.get_label_by_id('110').set_text('5518\n 9.6 %')
		# v.get_label_by_id('101').set_text('10651\n 18.5 %')
		# v.get_label_by_id('111').set_text('2848\n 4.95 %')
		# v.get_label_by_id('010').set_text('33060\n 57.5 %')
		# v.get_label_by_id('011').set_text('2923\n 5.08 %')
		# v.get_label_by_id('001').set_text('209\n 0.36 %')

		plt.title("Venn diagram - Amount of Frequency types for drug-side effect pair")

		#plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-100,-80),
		#             ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
		#             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
		#if len(directory)>0: 
			#plt.savefig(directory + 'venndiagram.svg')
		plt.show()
		return fig


		
	def preprocessingToFrequencyLabels(self,drug_se_pair):
		'''
			This function pre-process and filter the drug_se_pair, according to different criterias:
			
				Set A: exact and range side-effect frequency.
				Set B: frequency labels.
				Set C: placebo frequency.
			
			General Rule: all the pairs that are left and have frequency of the drug associated with it, will be converted
				a equivalent label according to WHO (World Health Organization recommendations).
				http://www.who.int/medicines/areas/quality_safety/safety_efficacy/trainingcourses/definitions.pdf
				
				1) C-all pairs deleted, because they have placebo but not the drug exact frequency.
				2) A & C, we deleted those pairs for which the median frequency of the placebo is bigger than the median
					frequency of the drug. 
				3) A-all we computed the median.
				4) A & B, compute the median of the frequency. Keep the labels from B.
				5) A & B & C, we can obtain later. same criteria but intersect.
				6) B & C, deleted the placebo frequency, keep the labels.
				7) B-all, keep the labels.
		'''
		operations = ["A-[B U C]","B-[A U C]", "C-[A U B]", "A & B", "A & C", "B & C", "A & B & C"]
		drug_se_fingerprint = dict()
		placeboGroup = defaultdict(list)
		for v in operations:
			drug_se_fingerprint[v] = dict()
			
		for drugSEpair, fiels in drug_se_pair.iteritems():
		
			# We compute the median frequencies if we can.
			A = False
			B = False
			C = False
			
			if 'exact_freq' in fiels or 'range_freq' in fiels: # set A
				A = True            
				freq = list()
				if 'exact_freq' in fiels:
					for v in fiels['exact_freq']:
						freq.append(v)
				if 'range_freq' in fiels:
					for v in fiels['range_freq']:
						freq.append(v) 
						
				freq.sort()
				MedDrug = self.__median(freq)
			
			if 'placebo_exact_freq' in fiels or 'placebo_range_freq' in fiels: # set C
				C = True  
				freqPlacebo = list()
				if 'placebo_exact_freq' in fiels:
					for v in fiels['placebo_exact_freq']:
						freqPlacebo.append(v)
				if 'placebo_range_freq' in fiels:
					for v in fiels['placebo_range_freq']:
						freqPlacebo.append(v) 
						
				freqPlacebo.sort()
				MedPlacebo = self.__median(freqPlacebo)
			
			if 'label_freq' in fiels: # set B 
				B = True
				
			# We can ask about the cases.
			if A & B:            
				if drugSEpair not in drug_se_fingerprint["A & B"]:
					drug_se_fingerprint["A & B"][drugSEpair] = defaultdict(list)
					
				
				drug_se_fingerprint["A & B"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(self.__MappingFreqtoLabel(MedDrug)))    
				
				for v in fiels['label_freq']:                   
					drug_se_fingerprint["A & B"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(v))
						
			if A & C: 
				if (MedDrug > MedPlacebo):
					if drugSEpair not in drug_se_fingerprint["A & C"]:
						drug_se_fingerprint["A & C"][drugSEpair] = defaultdict(list)
						
					
					drug_se_fingerprint["A & C"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(self.__MappingFreqtoLabel(MedDrug))) 
				else:
					placeboGroup[drugSEpair].append(self.__SideEffectFingerprint(self.__MappingFreqtoLabel(MedDrug)))
					placeboGroup[drugSEpair].append(self.__SideEffectFingerprint(self.__MappingFreqtoLabel(MedPlacebo)))
					
			if A & ~B & ~C:
				if drugSEpair not in drug_se_fingerprint["A-[B U C]"]:
					drug_se_fingerprint["A-[B U C]"][drugSEpair] = defaultdict(list)
					
				#drug_se_fingerprint["A-[B U C]"][drugSEpair]['MedianFreq'].append(MedDrug)
				drug_se_fingerprint["A-[B U C]"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(self.__MappingFreqtoLabel(MedDrug)))   
				
			if B & C:
				if drugSEpair not in drug_se_fingerprint["B & C"]:
					drug_se_fingerprint["B & C"][drugSEpair] = defaultdict(list)
					
				# we only grap the labels from B.
				for v in fiels['label_freq']:                   
					drug_se_fingerprint["B & C"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(v))
		   
			if C & ~A & ~B: # we are not interested in this case, but save it to see it.
				if drugSEpair not in drug_se_fingerprint["C-[A U B]"]:
					drug_se_fingerprint["C-[A U B]"][drugSEpair] = defaultdict(list)
					
				drug_se_fingerprint["C-[A U B]"][drugSEpair]['MedPlacebo'].append(MedPlacebo)
			
			if A & B & C:
				if drugSEpair not in drug_se_fingerprint["A & B & C"]:
					drug_se_fingerprint["A & B & C"][drugSEpair] = defaultdict(list)
					
				if (MedDrug > MedPlacebo):
						
					#drug_se_fingerprint["A & B & C"][drugSEpair]['MedianFreq'].append(MedDrug) 
					drug_se_fingerprint["A & B & C"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(self.__MappingFreqtoLabel(MedDrug)))   
					
				for v in fiels['label_freq']:                   
					drug_se_fingerprint["A & B & C"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(v))
						
			if B & ~A & ~C:
				if drugSEpair not in drug_se_fingerprint["B-[A U C]"]:
					drug_se_fingerprint["B-[A U C]"][drugSEpair] = defaultdict(list)
						
				for v in fiels['label_freq']:                   
					drug_se_fingerprint["B-[A U C]"][drugSEpair]['LabelFreq'].append(self.__SideEffectFingerprint(v))

		
		return drug_se_fingerprint, placeboGroup

	def UnifySets(self,drug_se_fingerprint):
		'''
			We take only the operations we want to join
		'''
		filter_pairs = defaultdict(set)
		operations = ["A-[B U C]","B-[A U C]", "A & B", "A & C", "B & C"]

		# For each operation
		for op in operations:
			# For each drug-se pair
			for pair in drug_se_fingerprint[op]:
				# Each pair suppose to be unique
				if pair not in filter_pairs:
					for freq in drug_se_fingerprint[op][pair]['LabelFreq']:
						if freq > 0: # Filter postmarketing side-effects (CAREFUL!!!).
							filter_pairs[pair].add(freq)
		return filter_pairs
	
	def UnifySetsPostmarketing(self, drug_se_fingerprint, unique_se_with_freq):
		'''
			We only return post-marketing side-effects that have no known frequency
		'''
		filter_pairs = defaultdict(set)
		operations = ["A-[B U C]","B-[A U C]", "A & B", "A & C", "B & C"]

		# For each operation
		for op in operations:
			# For each drug-se pair
			for pair in drug_se_fingerprint[op]:
				# Each pair suppose to be unique
				if pair not in filter_pairs:
					listfreq = drug_se_fingerprint[op][pair]['LabelFreq']
					if -1 in listfreq and len(listfreq) == 1: # only postmarketing						
						drugID, se = pair.split('|')
						if se in unique_se_with_freq:
							filter_pairs[drugID].add(se)
							
		return filter_pairs
	
	def SideEffectswithNoFrequency(self, all_se, se_with_freq, pmktg_se, unique_se_with_freq):
		'''
			return side-effects with no frequency information or postmarketing
		'''
		
		drugs_no_freq = defaultdict(list)
		
		for drug, listse in all_se.iteritems():
			for se in listse:
				# the drug should be in my list of interest
				if drug in se_with_freq and se in unique_se_with_freq:
					# only if the se does not have frequency nor postmarketing
					if se not in se_with_freq[drug]:
						if drug not in pmktg_se or se not in pmktg_se[drug]:
							drugs_no_freq[drug].append(se)
				
		return drugs_no_freq
		
	def removepairsInconsistencyFrequency(self,filter_pairs):
		'''
			In this method we count how many pairs have different labels.			
			UPDATE: we will average those that are different.
		'''
		
		count = dict()
		
		DrugSEFilter = defaultdict(list)
		AllFrequencies = list()
		
		for pair,values in filter_pairs.iteritems():
			nofreq = len(values)
			if nofreq not in count:
				count[nofreq] = 1
			else:
				count[nofreq] +=1
			
			#if nofreq > 3:
				#print pair, nofreq,values, '\n'
			
			# We only keep those of len = 1
			#if nofreq == 1:
				#for freq in values:
					#DrugSEFilter[pair].append(freq)
					#AllFrequencies.append(freq)
			DrugSEFilter[pair].append(np.mean(list(values)))
			AllFrequencies.append(np.mean(list(values)))
			#if nofreq > 3:
				#print DrugSEFilter[pair], '\n'	
				
		#print count  
		
		return DrugSEFilter, Counter(AllFrequencies)
	
	def finalDrugSElist(self, DrugSEFilter):
		'''
			DICT: 
				key: pubchem ID of the drug.
				value: DICT:
					key: side-effect
					value: frequency of the side-effect for the drug.
		'''
		drug_se_profile = dict() 
		unique_side_effects = set()

		for drugID , frequency in DrugSEFilter.iteritems():
			
			data = drugID.split('|')
			pubchemID = data[0]
			side_effect = data[1].lower().strip()
			unique_side_effects.add(side_effect)			
			
			if pubchemID not in drug_se_profile:
				drug_se_profile[pubchemID] = dict()
				
			drug_se_profile[pubchemID][side_effect] = frequency

		unique_side_effects = list(unique_side_effects)
	
		return drug_se_profile, unique_side_effects
	
	def __median(self, lst):
		'''
			This function allows to compute the median for a list of numbers.
		'''
		return numpy.median(numpy.array(lst))
		
	def __SideEffectFingerprint(self,label):
		'''
			This function receives as input the frequency of the side-effect and returns the corresponding value for it, 
			according to the WHO (Collaborating Centre for Drug Statistics Methodology).
						very common >= 10% 
				 1% <=  common or frequent < 10%  
			   0.1% <=  uncommon or infrequent < 1%
			  0.01% <=  rare     < 0.1% 
						very rare < 0.01%
					   
		'''
		value = -1 # postmarketing

		if label == 'veryfrequent':
			value = 5
		elif label == 'frequent':
			value = 4
		elif label == 'infrequent':
			value = 3
		elif label == 'rare':
			value = 2
		elif label == 'veryrare':
			value = 1
			   
		return value
		
	def __MappingFreqtoLabel(self,fq):
		'''
			This function receives as input the frequency of the side-effect and returns the corresponding label for it, 
			according to the WHO (Collaborating Centre for Drug Statistics Methodology).
						very common >= 10% 
				 1% <=  common or frequent < 10%  
			   0.1% <=  uncommon or infrequent < 1%
			  0.01% <=  rare     < 0.1% 
						very rare < 0.01%
					   
		'''
		label = 'error'
		if fq >= 10:
			label = 'veryfrequent'
		elif fq >= 1:
			label = 'frequent'
		elif fq >= 0.1:
			label = 'infrequent'
		elif fq >= 0.01:
			label = 'rare'
		elif fq < 0.01:
			label = 'veryrare'
			   
		return label

	def __NormalizeFrequencyLabels(self, fqLabel):

		if fqLabel == 'common':       
			fqLabel = 'frequent'
		elif fqLabel == 'uncommon':
			fqLabel = 'infrequent'   
		elif fqLabel == 'verycommon':
			fqLabel = 'veryfrequent'  
			
		return fqLabel  
		
		
if __name__ == '__main__':

	# Get the default directories
	parser = parsers.EasyParsers()
	data_directory = parser.get_data_directory()
	result_directory = parser.get_results_directory()
	
	# 
