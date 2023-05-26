#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 2019

@author: mduncans@umich.edu
"""
# Replace the following two variables with paths to bbmerge and flash.
BBMerge_path = '~/bbmap/bbmerge.sh'
FLASH_path = '~/FLASH-1.2.11/flash'

import os
import re
import sys
import glob
import tqdm
import itertools
import numpy as np
import pandas as pd
from functools import partial
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from Bio.Seq import translate, Seq
from CDR_identifier import cdr_finder
from multiprocessing import get_context

#NGS_sample classes are used to analyze the raw .fastq files from the core to generate seq:freq data.
#R1.fastq + R2.fastq --> merged.fastq --> merged.fasta --> seq: freq
class NGS_sample(object):
	#NGS_sample class analyzes the .fastq.gz output files from a MiSeq NextGen Sequencing Run. The output files will be used with NGS_round_data.
	def __init__(self, sequence_length, starting_sequence, output_path, qtrim, cdr_regions = None, mutation_locations = None, use_cdr_finder = False, cdr_finder_kwargs = None, sub_string = None):
		'''
		Input:
			sequence_length: <class 'int'> the length of the DNA sequence
			starting_sequence: <class 'str'> a string of the starting amino acids in the sequence.
			output_path: <class 'str'> the directory where the output will be saved.
			qtrim: <class 'int'> qtrim parameter for bbmerge
			cdr_regions: <class 'list'> List containing tuples of <class 'int'> representing start and end positions of cdrs.
			mutation_locations: <class 'list'> containing <class 'int'> representing positions of mutations within cdrs.
			use_cdr_finder: <class 'bool'> Specifies whether cdr_regions will be used to get cdr locations or CDR_identifier.cdr_finder
			cdr_finder_kwargs: <class 'dict'> keyword dictionary containing arguments to CDR_identifier.cdr_finder:
				{'vh_loc':[x,x+n], # or None specifies location of variable domain
				'vl_loc':None, # or [y, y+m]
				'definition':'kabat', 
				'which_CDRHs':[True, True, True], 
				'which_CDRLs':[False, False, False]}
			sub_string: <class 'list'> representing the start and end of the sequence to consider for analysis.
		'''
		if cdr_finder_kwargs is not None:
			from CDR_identifier import cdr_finder, split_sequences  

		if cdr_regions == None and mutation_locations == None and use_cdr_finder == False:
			print('You must supply either cdr_regions and mutation_regions or set use_cdr_finder to True')
			quit()
		if use_cdr_finder:
			if not cdr_finder_kwargs:
				print('must provide cdr finder arguments.')
				quit()
		self.sequence_length = sequence_length
		if output_path.endswith('/'):
			self.output_path = output_path
		else:
			self.output_path = output_path + '/'
		self.starting_sequence = starting_sequence.upper()
		self.qtrim = qtrim
		self.cdr_regions = cdr_regions
		self.mutation_locations = mutation_locations
		self.cdr_finder_kwargs = cdr_finder_kwargs
		self.use_cdr_finder = use_cdr_finder
		self.sub_string = sub_string

	def run_flash(self):
		'''
		runs flash on unzipped .fastq files to form a single merged .fastq file. 
		
		creates the following files containing merged sequences, unmerged sequences, and histogram data of merged read length
			'merged_6436-MS-1_barcode-barcode_SX.extendedFrags.fastq'
			'merged_6436-MS-1_barcode-barcode_SX.notCombined_1.fastq'
			'merged_6436-MS-1_barcode-barcode_SX.notCombined_2.fastq'
			'merged_6436-MS-1_barcode-barcode_SX.hist'
			'merged_6436-MS-1_barcode-barcode_SX.histogram'
		Finds and uses two read files in directory
			'6436-MS-1_barcode-barcode_SX_R1_001.fastq'
			'6436-MS-1_barcode-barcode_SX_R2_001.fastq'

		'''
		m = glob.glob('*merged*.fastq')
		if len(m) == 0:
			f = glob.glob('*.fastq')
			os.system(f"{FLASH_path} {f[0]} {f[1]} -M 301 -o merged_{f[0].split('_R')[0]}")

	def run_bbmerge(self):
		'''
		runs bbmerge on unzipped .fastq files to form one .fastq file for the sample
		default value for qtrim is currently 10. Should be altered based on FastQC results.

		qtrim sets the quality trimming metric for bbmerge

		creates two files containing merged sequences and unmerged sequences
			'Sample_x_merged.fastq'
			'Sample_x_unmerged.fastq'
		Finds and uses two read files in directory
			'Sample_x_R1_001.fastq'
			'Sample_x_R2_001.fastq'

		'''
		m = glob.glob('*merged.fastq')
		if len(m) == 0:
			f = glob.glob('*.fastq')
			os.system(f"{BBMerge_path} in1={f[0]} in2={f[1]} out={re.sub('_R[1,2]_001.fastq', '_merged.fastq', f[0])} -qtrim=15") #outu={re.sub('_R[1,2]_001.fastq', '_unmerged.fastq', f[0])}
			print('{} and {} successfully merged.'.format(f[0], f[1]))
		else:
			for fqm in m:
				os.system('rm {}'.format(fqm))
			f = glob.glob('*.fastq')
			os.system(f"{BBMerge_path} in1={f[0]} in2={f[1]} out={re.sub('_R[1,2]_001.fastq', '_merged.fastq', f[0])} -qtrim=15") #outu={re.sub('_R[1,2]_001.fastq', '_unmerged.fastq', f[0])} 
			print('{} and {} successfully merged.'.format(f[0], f[1]))

	def process_fasta_line(self, line):
		line = line.split('\n')[0]
		seq_num = 0
		dna_lines = 0

		seq_dict = {}
		cdr_dict = {}
		if self.mutation_locations is not None:
			mutation_dict = {}
		
		if not line.startswith('>'):
			if self.sub_string is not None:
				seq = line[self.sub_string[0]:self.sub_string[1]]
			else:
				seq = line
			dna_lines += 1
			if len(seq) >= self.sequence_length and 'N' not in seq and len(seq) % 3 == 0:
				aa_seq = translate(seq)
				if not aa_seq.startswith(self.starting_sequence):
					seq_class = Seq(seq)
					aa_seq = str(translate(str(seq_class.reverse_complement())))
				if '*' not in aa_seq:
					seq_num += 1
					
					try:
						seq_dict[aa_seq] += 1
					except:
						seq_dict[aa_seq] = 1
					
					cdr_str = ''
					if self.cdr_regions is None:
						if self.cdr_finder_kwargs['vh_loc'] is not None:
							VH = aa_seq[self.cdr_finder_kwargs['vh_loc'][0]:self.cdr_finder_kwargs['vh_loc'][1]]
						else:
							VH = None
						if self.cdr_finder_kwargs['vl_loc'] is not None:
							VL = aa_seq[self.cdr_finder_kwargs['vl_loc'][0]:self.cdr_finder_kwargs['vl_loc'][1]]
						else:
							VL = None
						try:
							cdr_finder_dict = cdr_finder(VH_seq = VH, VL_seq = VL, **self.cdr_finder_kwargs)
						except:
							return [seq_dict, {}, seq_num, dna_lines]
						cdr_str = ''
						for cdr, cdr_set in cdr_finder_dict.items():
							cdr_str += cdr_set[0] + '_'
					else:
						for cdr in self.cdr_regions:
							cdr_str += f'{aa_seq[cdr[0]:cdr[1]]}_'
					cdr_str = cdr_str[0:len(cdr_str) - 1]
					
					if self.mutation_locations is not None:
						mutation_str = ''
						for mutation in self.mutation_locations:
							mutation_str += cdr_str[mutation]

						try:
							mutation_dict[mutation_str] += 1
						except:
							mutation_dict[mutation_str] = 1

					try:
						cdr_dict[cdr_str] += 1
					except:
						cdr_dict[cdr_str] = 1
		if self.mutation_locations is not None:
			return [seq_dict, cdr_dict, mutation_dict, seq_num, dna_lines]
		else:
			return [seq_dict, cdr_dict, seq_num, dna_lines]

	def condense_dictionaries(self, dict_list):
		condensed_dict = {}
		for d in dict_list:
			for key, value in d.items():
				try:
					condensed_dict[key] += value
				except:
					condensed_dict[key] = value
		return condensed_dict

	def get_freq_fasta_mp(self):
		try:
			print('Analyzing fasta file now.')
			for fa in glob.glob('*.fasta'):
				lines = open(fa)
				if os.cpu_count() > 1:
					p2 = get_context('fork').Pool(os.cpu_count()-1)
				else:
					p2 = get_context('fork').Pool()
					
				results = p2.map(self.process_fasta_line, lines)
				p2.close()
				p2.join()

				s_dicts = [results[i][0] for i in range(len(results))]
				c_dicts = [results[i][1] for i in range(len(results))]
				if self.mutation_locations is not None:
					m_dicts = [results[i][2] for i in range(len(results))]

				if os.cpu_count() > 1:
					p = get_context('fork').Pool(os.cpu_count()-1)
				else:
					p = get_context('fork').Pool()
				if self.mutation_locations is not None:
					seq_dict, cdr_dict, mutation_dict = p.map(self.condense_dictionaries, [s_dicts, c_dicts, m_dicts])
				else:
					seq_dict, cdr_dict = p.map(self.condense_dictionaries, [s_dicts, c_dicts])

				p.close()
				p.join()

				seq_num = 0
				dna_lines = 0
				if self.mutation_locations is not None:
					for result in results:
						seq_num += result[3]
						dna_lines += result[4]
				else:
					for result in results:
						seq_num += result[2]
						dna_lines += result[3]

				print(f'Creating csv file now. {seq_num} correct reads found from {dna_lines} total sequences')

				count_outfile = open(f"{self.output_path}Count_{''.join(str(int(fa.split('_merged')[0].split('_')[2].split('S')[1]) + 1000))}.csv", 'w')
				count_outfile.write(f'{seq_num}')
				count_outfile.close()

				seq_outfile = open('{}Sequences_{}.csv'.format(self.output_path, ''.join(str(int(fa.split('_merged')[0].split('_')[2].split('S')[1]) + 1000))), 'w')
				for sequence in sorted(seq_dict, key=seq_dict.get, reverse=True):
					seq_outfile.write('{}, {}\n'.format(sequence, seq_dict[sequence]/float(seq_num)))
					seq_dict.pop(sequence)
				seq_outfile.close()

				cdr_outfile = open('{}CDR_{}.csv'.format(self.output_path, ''.join(str(int(fa.split('_merged')[0].split('_')[2].split('S')[1]) + 1000))), 'w')
				for cdr in sorted(cdr_dict, key=cdr_dict.get, reverse=True):
					cdr_outfile.write('{}, {}\n'.format(cdr, cdr_dict[cdr]/float(seq_num)))
					cdr_dict.pop(cdr)
				cdr_outfile.close()

				if self.mutation_locations is not None:
					mutation_outfile = open('{}Mutation_{}.csv'.format(self.output_path, ''.join(str(int(fa.split('_merged')[0].split('_')[2].split('S')[1]) + 1000))), 'w')
					
					for mutation in sorted(mutation_dict, key=mutation_dict.get, reverse=True):
						mutation_outfile.write('{}, {}\n'.format(mutation, mutation_dict[mutation]/float(seq_num)))
						mutation_dict.pop(mutation)

					mutation_outfile.close()
				print('CSV File successfully created.')
		except:
			try:
				p2.close()
				p2.join()
				return
			except:
				p.close()
				p.join()
				return

	def analyze_sample_mp(self, sample_input):
		'''
		input:
			sample_input: <class 'int'> sample number ids.
			sample_input: <class 'tuple'> sample id and y/n for fasta recreation.
		'''
		if type(sample_input) == tuple:
			sample = sample_input[0]
			y_n = sample_input[1]
		else:
			sample = sample_input
			y_n = 'y'
		print('Analyzing Sample {}'.format(sample))
		os.chdir('Sample_{}'.format(sample))

		print('unzipping .fastq.gz files')
		unzip()
		
		#self.run_flash()
		self.run_bbmerge()

		print('converting to fasta')
		convert(y_n)

		os.chdir('../')

	def analyze_all_samples_mp(self, samples):
		'''
		input:
			samples: <class 'list'> list containg sample number ids.
		'''
		if os.cpu_count() > 1:
			p1 = get_context('fork').Pool(os.cpu_count()-1)
		else:
			p1 = get_context('fork').Pool()
		p1.map(self.analyze_sample_mp, samples)
		p1.close()
		p1.join()

		for sample in samples:
			print(f'Calculating frequencies of Sample {sample}')
			os.chdir('Sample_{}'.format(sample))
			self.get_freq_fasta_mp()
			os.chdir('../')

class NGS_phased_sample(NGS_sample):
	#To be used with samples prepared with amplicon phasing. Output will be used with NGS_round_data.
	def __init__(self, sequence_length, starting_sequence, output_path, qtrim, cdr_regions = None, mutation_locations = None, use_cdr_finder = False, cdr_finder_kwargs = None, sub_string = None):
		'''
		Input:
			sequence_length: <class 'int'> the length of the DNA sequence
			starting_sequence: <class 'str'> a string of the starting amino acids in the sequence.
			output_path: <class 'str'> the directory where the output will be saved.
			qtrim: <class 'int'> qtrim parameter for bbmerge
			cdr_regions: <class 'list'> List containing tuples of <class 'int'> representing start and end positions of cdrs.
			mutation_locations: <class 'list'> containing <class 'int'> representing positions of mutations within cdrs.
			use_cdr_finder: <class 'bool'> Specifies whether cdr_regions will be used to get cdr locations or CDR_identifier.cdr_finder
			cdr_finder_kwargs: <class 'dict'> keyword dictionary containing arguments to CDR_identifier.cdr_finder:
				{'vh_loc':[x,x+n], # or None specifies location of variable domain
				'vl_loc':None, # or [y, y+m]
				'definition':'kabat', 
				'which_CDRHs':[True, True, True], 
				'which_CDRLs':[False, False, False]}
			sub_string: <class 'dict'> {sample_1: [a, b], sample_2: [c, d], ...}
		'''
		if cdr_finder_kwargs is not None:
			from CDR_identifier import cdr_finder, split_sequences  

		if cdr_regions == None and mutation_locations == None and use_cdr_finder == False:
			print('You must supply either cdr_regions and mutation_regions or set use_cdr_finder to True')
			quit()
		if use_cdr_finder:
			if not cdr_finder_kwargs:
				print('must provide cdr finder arguments.')
				quit()
		self.sequence_length = sequence_length
		if output_path.endswith('/'):
			self.output_path = output_path
		else:
			self.output_path = output_path + '/'
		self.starting_sequence = starting_sequence.upper()
		self.qtrim = qtrim
		self.cdr_regions = cdr_regions
		self.mutation_locations = mutation_locations
		self.cdr_finder_kwargs = cdr_finder_kwargs
		self.use_cdr_finder = use_cdr_finder
		self.sub_string = sub_string

	def process_phased_fasta_line(self, sample, line):
		seq_num = 0
		dna_lines = 0

		seq_dict = {}
		cdr_dict = {}
		if self.mutation_locations is not None:
			mutation_dict = {}
		
		#only look at sequence lines
		if not line.startswith('>'):
			seq = line.split('\n')[0]
			dna_lines += 1
			
			#Check if sequence is ~correct length
			if len(seq) >= self.sequence_length and len(seq) % 3 == 0:
				
				#check if first few amino acids are correct or have ~one mutation causes incorrect start
				if hamming_distance(translate(seq[self.sub_string[sample][0]:self.sub_string[sample][1]])[0:len(self.starting_sequence)], self.starting_sequence) <= 1:
					seq_num += 1
					aa_seq = translate(seq[self.sub_string[sample][0]:self.sub_string[sample][1]])
					try:
						seq_dict[aa_seq] += 1
					except:
						seq_dict[aa_seq] = 1

				#check is seq is shifted upstream 1 bp
				elif hamming_distance(translate(seq[self.sub_string[sample][0] - 1:self.sub_string[sample][1]])[0:len(self.starting_sequence)], self.starting_sequence) <= 1:
					seq_num += 1
					aa_seq = translate(seq[self.sub_string[sample][0] - 1:self.sub_string[sample][1]])
					try:
						seq_dict[aa_seq] += 1
					except:
						seq_dict[aa_seq] = 1
				
				#check if seq is shifted downstream 1 bp
				elif hamming_distance(translate(seq[self.sub_string[sample][0] + 1:self.sub_string[sample][1]])[0:len(self.starting_sequence)], self.starting_sequence) <= 1:
					seq_num += 1
					aa_seq = translate(seq[self.sub_string[sample][0] + 1:self.sub_string[sample][1]])
					try:
						seq_dict[aa_seq] += 1
					except:
						seq_dict[aa_seq] = 1
				
				#Use reverse complement if every above check failed.
				else:
					seq = Seq(seq).reverse_complement()

					#check if first few amino acids are correct or have ~one mutation causes incorrect start
					if hamming_distance(translate(seq[self.sub_string[sample][0]:self.sub_string[sample][1]])[0:len(self.starting_sequence)], self.starting_sequence) <= 1:
						seq_num += 1
						aa_seq = translate(seq[self.sub_string[sample][0]:self.sub_string[sample][1]])
						try:
							seq_dict[aa_seq] += 1
						except:
							seq_dict[aa_seq] = 1

					#check is seq is shifted upstream 1 bp
					elif hamming_distance(translate(seq[self.sub_string[sample][0] - 1:self.sub_string[sample][1]])[0:len(self.starting_sequence)], self.starting_sequence) <= 1:
						seq_num += 1
						aa_seq = translate(seq[self.sub_string[sample][0] - 1:self.sub_string[sample][1]])
						try:
							seq_dict[aa_seq] += 1
						except:
							seq_dict[aa_seq] = 1
					
					#check if seq is shifted downstream 1 bp
					elif hamming_distance(translate(seq[self.sub_string[sample][0] + 1:self.sub_string[sample][1]])[0:len(self.starting_sequence)], self.starting_sequence) <= 1:
						seq_num += 1
						aa_seq = translate(seq[self.sub_string[sample][0] + 1:self.sub_string[sample][1]])
						try:
							seq_dict[aa_seq] += 1
						except:
							seq_dict[aa_seq] = 1					
				
				#Convert aa_seq to cdr_str and mut_str
				if 'aa_seq' in locals():
					try:
						cdr_str = ''
						if self.cdr_regions is None:
							if self.cdr_finder_kwargs['vh_loc'] is not None:
								VH = aa_seq[self.cdr_finder_kwargs['vh_loc'][0]:self.cdr_finder_kwargs['vh_loc'][1]]
							else:
								VH = None
							if self.cdr_finder_kwargs['vl_loc'] is not None:
								VL = aa_seq[self.cdr_finder_kwargs['vl_loc'][0]:self.cdr_finder_kwargs['vl_loc'][1]]
							else:
								VL = None
							try:
								cdr_finder_dict = cdr_finder(VH_seq = VH, VL_seq = VL, **self.cdr_finder_kwargs)
							except:
								return [seq_dict, {}, seq_num, dna_lines]
							cdr_str = ''
							for cdr, cdr_set in cdr_finder_dict.items():
								cdr_str += cdr_set[0] + '_'
						
						else:
							for cdr in self.cdr_regions:
								cdr_str += f'{aa_seq[cdr[0]:cdr[1]]}_'
						cdr_str = cdr_str[0:len(cdr_str) - 1]
						
						try:
							cdr_dict[cdr_str] += 1
						except:
							cdr_dict[cdr_str] = 1
					except Exception as e:
						#print('Errored out in CDR_STR')
						#print(f'Errored out due to {e}\n \t')
						exc_type, exc_obj, exc_tb = sys.exc_info()
						#print(exc_type, exc_obj, exc_tb)
						
					if self.mutation_locations is not None:
						try:
							mutation_str = ''
							for mutation in self.mutation_locations:
								mutation_str += cdr_str[mutation]

							try:
								mutation_dict[mutation_str] += 1
							except:
								mutation_dict[mutation_str] = 1
						except Exception as e:	
							#print(f'Errored out on MUT_STR, {aa_seq}, {cdr_str}, {mutation}')
							#print(f'Errored out due to {e}\n \t')
							exc_type, exc_obj, exc_tb = sys.exc_info()
							#print(exc_type, exc_obj, exc_tb)

		if self.mutation_locations is not None:
			return [seq_dict, cdr_dict, mutation_dict, seq_num, dna_lines]
		else:
			return [seq_dict, cdr_dict, seq_num, dna_lines]
	
	def get_freq_fasta_phased_mp(self, sample):
		try:
			print('Analyzing fasta file now.')
			for fa in glob.glob('*.fasta'):
				lines = open(fa)
				if os.cpu_count() > 1:
					p2 = get_context('fork').Pool(os.cpu_count()-1)
				else:
					p2 = get_context('fork').Pool()
				
				func = partial(self.process_phased_fasta_line, sample)
				results = p2.map(func, lines)
				p2.close()
				p2.join()
				print('FASTA processing done.')
				s_dicts = [results[i][0] for i in range(len(results))]
				c_dicts = [results[i][1] for i in range(len(results))]

				if self.mutation_locations is not None:
					m_dicts = [results[i][2] for i in range(len(results))]
					
				print('Condensing output now.')
				if os.cpu_count() > 1:
					p = get_context('fork').Pool(os.cpu_count()-1)
				else:
					p = get_context('fork').Pool()

				if self.mutation_locations is not None:
					seq_dict, cdr_dict, mutation_dict = p.map(self.condense_dictionaries, [s_dicts, c_dicts, m_dicts])
					print('Dictionaries concatenated.')
				else:
					seq_dict, cdr_dict = p.map(self.condense_dictionaries, [s_dicts, c_dicts])

				p.close()
				p.join()

				seq_num = 0
				dna_lines = 0
				if self.mutation_locations is not None:
					for result in results:
						seq_num += result[3]
						dna_lines += result[4]
				else:
					for result in results:
						seq_num += result[2]
						dna_lines += result[3]

				print(f'Creating csv file now. {seq_num} correct reads found from {dna_lines} total sequences')

				count_outfile = open(f"{self.output_path}Count_{sample + 1000}.csv", 'w')
				count_outfile.write(f'{seq_num}')
				count_outfile.close()

				seq_outfile = open(f"{self.output_path}Sequences_{sample + 1000}.csv", 'w')
				for sequence in sorted(seq_dict, key=seq_dict.get, reverse=True):
					seq_outfile.write('{}, {}\n'.format(sequence, seq_dict[sequence]/float(seq_num)))
					seq_dict.pop(sequence)
				seq_outfile.close()

				cdr_outfile = open(f"{self.output_path}CDR_{sample + 1000}.csv", 'w')
				for cdr in sorted(cdr_dict, key=cdr_dict.get, reverse=True):
					cdr_outfile.write('{}, {}\n'.format(cdr, cdr_dict[cdr]/float(seq_num)))
					cdr_dict.pop(cdr)
				cdr_outfile.close()

				if self.mutation_locations is not None:
					mutation_outfile = open(f"{self.output_path}Mutation_{sample + 1000}.csv", 'w')
					for mutation in sorted(mutation_dict, key=mutation_dict.get, reverse=True):
						mutation_outfile.write('{}, {}\n'.format(mutation, mutation_dict[mutation]/float(seq_num)))
						mutation_dict.pop(mutation)
					mutation_outfile.close()
				print('CSV File successfully created.')
		
		except Exception as e:
			#print(f'Errored out due to {e}\n \t')
			exc_type, exc_obj, exc_tb = sys.exc_info()
			#print(exc_type, exc_obj, exc_tb)
			try:
				p2.close()
				p2.join()
				return
			except:
				p.close()
				p.join()
				return

	def analyze_all_samples_phased_mp(self, samples):
		'''
		input:
			samples: <class 'list'> list containg sample number ids or tuple of sample id and y/n for fasta.
		'''
		if os.cpu_count() > 1:
			p1 = get_context('fork').Pool(os.cpu_count()-1)
		else:
			p1 = get_context('fork').Pool()
		p1.map(self.analyze_sample_mp, samples)
		p1.close()
		p1.join()

		for sample in samples:
			if type(sample) == tuple:
				sample = sample[0]
			print(f'Calculating frequencies of Sample {sample}')
			os.chdir('Sample_{}'.format(sample))
			self.get_freq_fasta_phased_mp(sample)
			os.chdir('../')

#NGS_round_data Classes. used to contain sequence:frequency data for each sample.
class NGS_round_data(object):
	#to be used for analysis of anlyzed ngs data. 
	def __init__(self, Round, sequence_type, samples, sample_of_interest, path, wild_type = None, mutations_dict = None):
		'''
		A class containing all frequency data from NGS sequencing experiments and to analyze the results.

		NGS data should be stored in same directory as this file. It should be stored in a sub directory Rx/ where x is the round.
		Inside of Rx/ there should be several files.
				Sequences_XXX.csv
				Mutation_XXX.csv
				CDR_XXX.csv
				RX_sort_counts.csv
				these files will be created upon analyzing the NGS data.

		NGS_round initializes with
				Round: <class 'int'> representing the round the samples come from.
				sequence_type: <class 'str'> representing what type of sequence the data was analyzed as.
						Full sequence: sequences
						CDR sequence: CDRs
						only mutated residues: mutations
				samples: <class 'list'> list of the corresponding samples of each sequencing file. The order should be the same as the sequencing files
				sample_of_interest: <class 'str'> contains the primary antigen of interest (e.g. 'ABF')
				path: <class 'str'> path to the frequency and count data
				wild_type: <class 'str'> wild type sequence/CDR_str/mutation of wild type clone if applicable
				mutations_dict: <class 'dict'> dictionary of non-wildtype amino acids muated to. len(mutations_dict.keys()) = len(mutation_str)

		NGS_round will have the following attributes
				self.round: <class 'int'> corresponding round
				self.sample_counts: <class 'dict'> {sample: num_seq_in_sample}
				self.data: <class 'dict'> contains all frequency data for each sample
						{
								sample1: {seq1: freq_sample1_seq1, ..., seqN: freq_sample1_seqN}
								sample2: {seq1: freq_sample2_seq1, ..., seqM: freq_sample2_seqM}
								...
								sampleK: {seq1: freq_sampleK_seq1, ..., seqJ: freq_sampleK_seqJ}
						}
				self.sample_differences: <class 'dict'> contains number of sequences in a sample that aren't observed in sample_of_interest
						{
								sample_of_interest: 0,
								sample2: 5000,
								...
								sampleN: 20942
						}
				self.count_data: <class 'dict'> unused, but contains the distribution of counts for each sample
		'''
		self.round = Round
		if path.endswith('/'):
			self.path = path
		else:
			self.path = path+'/'
		self.samples = samples
		self.sample_counts = {line.split('\n')[0].split(',')[0]: float(line.split('\n')[0].split(',')[1]) for line in open(f'{self.path}R{self.round}/R{self.round}_sort_counts.csv')}
		self.sample_of_interest = sample_of_interest
		self.sequence_type = sequence_type
		self.current_path = os.getcwd()

		if self.sequence_type == 'sequences':
			self.data = self.get_round_dicts('{}R{}/'.format(self.path, self.round), 'Seq*.csv')
		elif self.sequence_type == 'mutations':
			self.data = self.get_round_dicts('{}R{}/'.format(self.path, self.round), 'Mu*.csv')
		elif self.sequence_type == 'CDRs':
			self.data = self.get_round_dicts('{}R{}/'.format(self.path, self.round), 'CD*.csv')

		self.sample_differences = {}
		for sample in self.samples:
			self.sample_differences[sample] = len(self.data[self.sample_of_interest].keys()) - len(set(self.data[self.sample_of_interest].keys()).intersection(set(self.data[sample].keys())))
			
		self.count_data = {}
		for sample in self.samples:
			self.count_data[sample] = []
			for freq in self.data[sample].values():
				self.count_data[sample].append(float(freq*float(self.sample_counts[sample])))

		if wild_type is not None:
			self.wt = wild_type

		if mutations_dict is not None:
			self.mutations_dict = mutations_dict

	def get_seq_dict(self, path):
		'''
				get_seq_dict creates a dictionary of seq: freq from analyzed NGS files.
				input:
						path: path to NGS analyzed csvs.

				output:
						sequence dictionary
						{seq: freq}
		'''
		seq_dict = {}

		for line in open(path, 'r'):
			seq_dict[line.split(',')[0]] = float(line.split(',')[1].split('\n')[0])
		return seq_dict

	def get_round_dicts(self, path, file_type):
		'''
				calls get_seq_dict on all samples within a round and creates a dictionary of all of the individual dictionaries.
		'''
		os.chdir(path)

		muts = glob.glob(file_type)
		muts = sorted(muts, key = lambda f: int(f.split('_')[1].split('.')[0]))
		print(muts)

		if len(muts) != len(self.samples):
			print('The number of samples does not match the number of sample files. Please update samples list.\n')
		if os.cpu_count() > 1:
			p = get_context('fork').Pool(os.cpu_count() - 1)
		else:
			p = get_context('fork').Pool()

		results = p.map(self.get_seq_dict, muts)
		p.close()
		p.join()
		
		r = {}
		for i in range(len(muts)):
			r[self.samples[i]] = results[i]

		os.chdir(self.current_path)
		return r

	def alter_data_dict(self):
		'''
				alter_data_dict adds sequences to sample dictionaries if that sequence is observed in sample_of_interest
				data_added contains dictionaries with increased number of sequences and the frequencies are scalled to sum to unity.
		output:
						self.data_added = {
								Sample1: {seq1: freq, ...}
						}

		'''
		if self.sample_of_interest in self.samples:

			self.data_added = {}
			for sample in self.samples:
				self.data_added[sample] = {}
				# iterates through sequences in sample of interest and another sample
				for seq in set(self.data[self.sample_of_interest].keys()).union(set(self.data[sample].keys())):
					# if seq is only in sample of interest add to data_added[sample]
					if seq not in self.data[sample]:
						self.data_added[sample][seq] = 1/(float(self.sample_counts[sample]) + float(self.sample_differences[sample]))

					else:
						self.data_added[sample][seq] = (self.data[sample][seq] * float(self.sample_counts[sample]))/(float(self.sample_counts[sample]) + float(self.sample_differences[sample]))

		else:
			print('Please check sample_of_interest parameter. {} was not found in sample list'.format(self.sample_of_interest))
			print('Acceptable values are {}'.format(self.samples))

	def condense_data(self, sample_requirements, threshold = 0):
		'''
				condense_data creates another dictionary that only contains sequences that are observed in all of the sample_requirements.
				for non required sorts data is pulled from data_added dictionary even if it is observed in that sort.
				input:
						sample_requirements: <class 'list'> sets which samples sequences need to be observed in.
						threshold: <class 'int'> number of observations of each sequence in each sample. default = 0
				output:
						self.condensed_data = {
								seq: [freq1, freq2, ..., freqN]
								...
						}
						self.condensed_data_order contains the order of how the data was entered.
		'''
		self.alter_data_dict()
		self.condensed_data = {}
		self.condensed_data_order = []
		# checks to ensure sample_requirements is a list and if not populates it correctly.
		if type(sample_requirements) != list:
			print('Enter sample_requirements as list not {}'.format(type(sample_requirements)))
			list_not_done = True
			sample_requirements = []
			while list_not_done:
				required_sample = input('Enter required sample now. Type XXX to stop entry.\n')
				if required_sample != 'XXX':
					if required_sample != '':
						sample_requirements.append(required_sample)
				else:
					list_not_done = False

		seqs = common_clones({self: sample_requirements}, threshold)

		for seq in seqs:
			for sample in self.samples:

				if sample not in self.condensed_data_order:
					self.condensed_data_order.append(sample)

				if sample in sample_requirements:
					if seq in self.data[sample].keys():
						try:
							self.condensed_data[seq].append(self.data[sample][seq])
						except:
							self.condensed_data[seq] = [self.data[sample][seq]]
					else:
						break
				else:
					if seq in self.data_added[sample].keys():
						try:
							self.condensed_data[seq].append(self.data_added[sample][seq])
						except:
							self.condensed_data[seq] = [self.data_added[sample][seq]]
					else:
						break

	def get_enrichment_ratios(self, sample_requirements, threshold = 0):
		'''
				calculates enrichment ratios. log2(freq_sort/freq_input)
				input:
						sample_requirements: <class 'list'> used to create self.condensed_data for calculations
				output:
						self.er = {
								seq: [log2(freq_i/freq_input) ...]
								...
						}
		'''
		self.condense_data(sample_requirements, threshold)

		try:
			input_index = self.condensed_data_order.index('Input')
		except:
			print('Input frequency not found.')

		self.er = {}
		self.er_order = []

		for seq, freq_list in self.condensed_data.items():
			self.er[seq] = []
			for i in range(len(freq_list)):
				if i != input_index:
					if self.condensed_data_order[i] not in self.er_order:
						self.er_order.append(self.condensed_data_order[i])
					self.er[seq].append(np.log2(freq_list[i] / freq_list[input_index]))
	
	def mutational_analysis(self, num_mutations, sample_requirements = None, er_indecies = None, clone_set = None, threshold = 0, recalculate = False):
		'''
		Performs mutational analysis on a NGS_round_data class. 

		Code should be optimized. Specifically creating mut_seqs_dict and wt_seqs_dict. Could probably alter code to find which positions are different from wt
		and only search mutation sets with the same differences. multiprocessing might also be added.
		'''
		
		#get locations of possible mutations
		
		try:
			eval(f'self.mutations_dictionary_{num_mutations}')
			if recalculate:
				raise ValueError
		except:
			mutation_str_length = len(self.wt)
			mut_loc_dict  = {i: loc for i, loc in enumerate(list(itertools.combinations(list(range(mutation_str_length)), num_mutations)))}

			#create possible mutations dictionary -- Fast
			mutation_number = 0
			exec(f'self.mutations_dictionary_{num_mutations} = dict()')
			for loc in tqdm.tqdm(mut_loc_dict.values()):
				mutations_dictionary_subkey = ''
				for i in loc:
					mutations_dictionary_subkey += str(i) + ','
				key_list = [i for i in mutations_dictionary_subkey.split(',') if i != '']
				product_list = [self.mutations_dict[int(j)] for j in key_list]
				possible_mutations = list(itertools.product(*product_list))
				for m in possible_mutations:
					exec(f'self.mutations_dictionary_{num_mutations}[{mutation_number}] = dict()')
					exec(f'self.mutations_dictionary_{num_mutations}[{mutation_number}][{mutations_dictionary_subkey}] = {m}')
					mutation_number += 1
		
		if sample_requirements is not None:
			try:
				eval(f'self.mut_seqs_{num_mutations}')
				eval(f'self.wt_seqs_{num_mutations}')
				if recalculate:
					raise ValueError
			except:
				#calculate enrichment ratios
				try:
					self.er
				except:
					self.get_enrichment_ratios(sample_requirements, threshold)

				#calculate mut_seq_dict and wt_seq_dict  slowest part currently
				mut_seqs = {m: [] for m in eval(f'self.mutations_dictionary_{num_mutations}.keys()')}
				wt_seqs = {m: [] for m in eval(f'self.mutations_dictionary_{num_mutations}.keys()')}
				exec(f'self.mut_seqs_{num_mutations} = {mut_seqs}')
				exec(f'self.wt_seqs_{num_mutations} = {wt_seqs}')

				if clone_set is not None:
					seqs = clone_set
				else:
					seqs = common_clones({self: sample_requirements}, threshold)

				seqs_dict = {i: [aa for aa in s] for i, s in enumerate(seqs)}
				clone_df = pd.DataFrame.from_dict(seqs_dict, orient = 'index')

				for m, mut_dict in tqdm.tqdm(eval(f'self.mutations_dictionary_{num_mutations}.items()')):
					for loc, muts in mut_dict.items():
						mut_sub = clone_df[clone_df.loc[:, loc] == muts].loc[:, loc].dropna()
						wts = [self.wt[i] for i in loc]
						wt_sub = clone_df[clone_df.loc[:, loc] == wts].loc[:, loc].dropna()


						for c in mut_sub.index:
							seq = ''.join(clone_df.loc[c, :])
							eval(f'self.mut_seqs_{num_mutations}[m].append(seq)')
						for c in wt_sub.index:
							seq = ''.join(clone_df.loc[c, :])
							eval(f'self.wt_seqs_{num_mutations}[m].append(seq)')
				
			if er_indecies is not None:
				try:
					eval(f'self.spearman_coefficients_{num_mutations}')
					if recalculate:
						raise ValueError
				except:
					#calculate spearman info
					exec(f'self.spearman_coefficients_{num_mutations} = dict()')		
					for index in tqdm.tqdm(er_indecies):
						exec(f'self.spearman_coefficients_{num_mutations}[self.er_order[index]] = dict()')
						for m in tqdm.tqdm(eval(f'self.mutations_dictionary_{num_mutations}.keys()')): 
							data = {}
							for seq in eval(f'self.mut_seqs_{num_mutations}[m]'):
								try:
									data[self.er[seq][index]][0] += 1
								except:
									data[self.er[seq][index]] = [1, 0]
							for seq in eval(f'self.wt_seqs_{num_mutations}[m]'):
								try:
									data[self.er[seq][index]][1] += 1
								except:
									data[self.er[seq][index]] = [0, 1]

							er_data = []
							freq_data = []
							for e, nums in data.items():
								er_data.append(e)
								freq_data.append(nums[0] / (sum(nums)))
							exec(f'self.spearman_coefficients_{num_mutations}[self.er_order[index]][m] = list(spearmanr(er_data, freq_data))')

	def process_seq_mut_analysis(self, product_entry):
		seq = product_entry[0]
		m = product_entry[1][0]
		mut_dict = product_entry[1][1]
		num_mutations = product_entry[2]
		pbar = product_entry[3]
		
		if hamming_distance(seq, self.wt) < num_mutations:
			return

		for loc, muts in mut_dict.items():
			mut_count = 0
			wt_count = 0
			for i, l in enumerate(list(loc)):
				if seq[int(l)] == muts[i]:
					mut_count += 1
				elif seq[int(l)] == self.wt[int(l)]:
					wt_count += 1
			if mut_count == num_mutations:
				eval(f'self.mut_seqs_{num_mutations}[{m}].append(seq)')
			elif wt_count == num_mutations:
				eval(f'self.wt_seqs_{num_mutations}[{m}].append(seq)')
		pbar.update(1)

	def mutational_analysis_mp(self, num_mutations, sample_requirements = None, er_indecies = None, clone_set = None, threshold = 0, recalculate = False):
		'''
		Performs mutational analysis on a NGS_round_data class. 

		Code should be optimized. Specifically creating mut_seqs_dict and wt_seqs_dict. Could probably alter code to find which positions are different from wt
		and only search mutation sets with the same differences. multiprocessing might also be added.
		'''
		
		#get locations of possible mutations
		
		try:
			eval(f'self.mutations_dictionary_{num_mutations}')
			if recalculate:
				raise ValueError
		except:
			mutation_str_length = len(self.wt)
			mut_loc_dict  = {i: loc for i, loc in enumerate(list(itertools.combinations(list(range(mutation_str_length)), num_mutations)))}

			#create possible mutations dictionary -- Fast
			mutation_number = 0
			exec(f'self.mutations_dictionary_{num_mutations} = dict()')
			for loc in tqdm.tqdm(mut_loc_dict.values()):
				mutations_dictionary_subkey = ''
				for i in loc:
					mutations_dictionary_subkey += str(i) + ','
				key_list = [i for i in mutations_dictionary_subkey.split(',') if i != '']
				product_list = [self.mutations_dict[int(j)] for j in key_list]
				possible_mutations = list(itertools.product(*product_list))
				for m in possible_mutations:
					exec(f'self.mutations_dictionary_{num_mutations}[{mutation_number}] = dict()')
					exec(f'self.mutations_dictionary_{num_mutations}[{mutation_number}][{mutations_dictionary_subkey}] = {m}')
					mutation_number += 1
		
		if sample_requirements is not None:
			try:
				eval(f'self.mut_seqs_{num_mutations}')
				eval(f'self.wt_seqs_{num_mutations}')
				if recalculate:
					raise ValueError
			except:
				#calculate enrichment ratios
				try:
					self.er
				except:
					self.get_enrichment_ratios(sample_requirements, threshold)

				#calculate mut_seq_dict and wt_seq_dict  slowest part currently
				mut_seqs = {m: [] for m in eval(f'self.mutations_dictionary_{num_mutations}.keys()')}
				wt_seqs = {m: [] for m in eval(f'self.mutations_dictionary_{num_mutations}.keys()')}
				exec(f'self.mut_seqs_{num_mutations} = {mut_seqs}')
				exec(f'self.wt_seqs_{num_mutations} = {wt_seqs}')

				if clone_set is not None:
					seqs = clone_set
				else:
					seqs = common_clones({self: sample_requirements}, threshold)

				pbar = tqdm.tqdm()
				entries = list(itertools.product(seqs, eval(f'self.mutations_dictionary_{num_mutations}.items()'), [num_mutations], [pbar]))
				pbar.total = len(entries)
				if os.cpu_count() > 1:
					p = get_context('fork').Pool(os.cpu_count() - 1)
				else:
					p = get_context('fork').Pool()
				p.map(self.process_seq_mut_analysis, entries)
				p.join()
				p.close()
				pbar.close()

			if er_indecies is not None:
				try:
					eval(f'self.spearman_coefficients_{num_mutations}')
					if recalculate:
						raise ValueError
				except:
					#calculate spearman info
					exec(f'self.spearman_coefficients_{num_mutations} = dict()')		
					for index in tqdm.tqdm(er_indecies):
						exec(f'self.spearman_coefficients_{num_mutations}[self.er_order[index]] = dict()')
						for m in tqdm.tqdm(eval(f'self.mutations_dictionary_{num_mutations}.keys()')): 
							data = {}
							for seq in eval(f'self.mut_seqs_{num_mutations}[m]'):
								try:
									data[self.er[seq][index]][0] += 1
								except:
									data[self.er[seq][index]] = [1, 0]
							for seq in eval(f'self.wt_seqs_{num_mutations}[m]'):
								try:
									data[self.er[seq][index]][1] += 1
								except:
									data[self.er[seq][index]] = [0, 1]

							er_data = []
							freq_data = []
							for e, nums in data.items():
								er_data.append(e)
								freq_data.append(nums[0] / (sum(nums)))
							exec(f'self.spearman_coefficients_{num_mutations}[self.er_order[index]][m] = list(spearmanr(er_data, freq_data))')

	def plot(self, data_type, samples, save_fig = False, fig_path = None, font_prop = None, marker_size = None, figsize = None, plot_wt = False):
		'''
		Plots two samples against each other.
		Inputs:
			data_type: <class 'str'> either er or dr sets what value selection is done on
			samples: <class 'list'> two member list of two samples to be plotted against each other samples[0] is x samples[1] is y
			plot_tested_clones: <class 'bool'> plots selected clones
			text_size: <class 'float'> alters figure text size if used.
			marker_size: <class 'float'> alters marker size
		'''
		try:
			eval('self.{}'.format(data_type))
		except:
			print('Data not found. Creating now.\n')
			sample_requirements = ''
			if data_type == 'er':
				self.get_enrichment_ratios(sample_requirements)
			elif data_type == 'dr':
				self.get_display_ratios(sample_requirements)

		xdata = []
		ydata = []

		indecies = [eval('self.{}_order'.format(data_type)).index(samples[0]), eval('self.{}_order'.format(data_type)).index(samples[1])]
		
		# WT Ploting
		if plot_wt:
			wt_dict = {self.wt: []}
		
		else:
			wt_dict = {}

		for seq in eval('self.{}.keys()'.format(data_type)):
			if not seq in wt_dict.keys():
				xdata.append(eval('self.{}'.format(data_type))[seq][indecies[0]])
				ydata.append(eval('self.{}'.format(data_type))[seq][indecies[1]])
			else:
				wt_dict[seq].append(eval('self.{}'.format(data_type))[seq][indecies[0]])
				wt_dict[seq].append(eval('self.{}'.format(data_type))[seq][indecies[1]])

		wt_non_empty = {}
		try:
			if len(wt_dict[self.wt]) == 2:
				wt_non_empty[self.wt] = wt_dict[self.wt]
		except:
			print('WT Not found in both samples')
		
		if figsize is not None:
			fig = plt.figure(figsize = figsize)
		else:
			fig = plt.figure()

		ax = fig.add_subplot(111)
		ax.vlines(0, min(ydata)*0.8, max(ydata)*1.2, linewidth=2, color='black')
		ax.hlines(0, min(xdata)*0.8, max(xdata)*1.2, linewidth=2, color='black')
		ax.set_xlim(min(xdata)*0.8, max(xdata)*1.2)
		ax.set_ylim(min(ydata)*0.8, max(ydata)*1.2)
		
		if marker_size is not None:
			ax.plot(xdata, ydata, '.', color='blue', markersize=marker_size, markerfacecolor = 'none')
		else:
			ax.plot(xdata, ydata, '.', color = 'blue', markerfacecolor = 'none')

		if plot_wt:
			if marker_size is not None:
				try:
					ax.plot(wt_non_empty[self.wt][0], wt_non_empty[self.wt][1], '.', color = 'red', label = plot_wt, markersize = marker_size, markeredgecolor = 'black')
					ax.legend()
				except:
					pass
			else:
				try:
					ax.plot(wt_non_empty[self.wt][0], wt_non_empty[self.wt][1], '.', color = 'red', label = plot_wt, markeredgecolor = 'black')
					ax.legend()
				except:
					pass
			

		if data_type == 'er':
			if font_prop is not None:
				ax.set_xlabel('{} Enrichment Ratio'.format(samples[0]), fontproperties = font_prop)
				ax.set_ylabel('{} Enrichment Ratio'.format(samples[1]), fontproperties = font_prop)

		elif data_type == 'dr':
			ax.set_xlabel('{} Display Ratio'.format(samples[0]))
			ax.set_ylabel('{} Display Ratio'.format(samples[1]))

		#ax.grid()
		if not save_fig:
			plt.show()
		else:
			plt.savefig(f'{fig_path}/R{self.round}_{samples[0]}_{samples[1]}.png', dpi=300)
		plt.close()	

	def get_one_counts(self):
		self.one_counts = {}		
		for sample in self.samples:
			for freq in self.data[sample].values():
				if freq * self.sample_counts[sample] < 1.5:
					try:
						self.one_counts[sample] += 1
					except:
						self.one_counts[sample] = 1

		for s, c in self.one_counts.items():
			self.one_counts[s] = round(c / self.sample_counts[s], 3)
			print(f'{s} \t {round(self.one_counts[s] * 100, 1)}%')
	
	def plot_count_histograms(self, num_rows, num_cols, font_size = 10, figsize = (13,9)):
		fig = plt.figure(figsize = figsize)
		for i, sample in enumerate(self.samples):
			plt.rcParams.update({'font.size': font_size})
			ax = fig.add_subplot(num_rows, num_cols, i + 1)
			data = [self.data[sample][seq] * self.sample_counts[sample] for seq in self.data[sample].keys()] 
			ax.hist(data, bins=range(int(min(data)), 100, 1))
			ax.set_xlabel('Sequence Counts')
			ax.set_ylabel('Count')
			ax.set_title(f'{sample}')

class NGS_DMD(object):

	def __init__(self, sequence_type, samples, path, wild_type, mutations_dict):
		'''
		A class for analyzing time series NGS data for use with DMD

		NGS data should be stored in a sub directory Rx/ where x is the round.
		Inside of Rx/ there should be several files.
				Sequences_XXX.csv
				Mutation_XXX.csv
				CDR_XXX.csv
				RX_sort_counts.csv
				these files will be created upon analyzing the NGS data with NGS_sample classes.

		NGS_round initializes with
				sequence_type: <class 'str'> representing what type of sequence the data was analyzed as.
						Full sequence: sequences
						CDR sequence: CDRs
						only mutated residues: mutations
				samples: <class 'list'> list of the corresponding samples of each sequencing file. The order should be the same as the sequencing files
				path: <class 'str'> path to the frequency and count data
				wild_type: <class 'str'> wild type sequence/CDR_str/mutation of wild type clone if applicable
				mutations_dict: <class 'dict'> dictionary of non-wildtype amino acids muated to. len(mutations_dict.keys()) = len(mutation_str)

		NGS_round will have the following attributes
				self.path
				self.samples
				self.sequence_type
				self.current_path
				self.wt
				self.mutations_dict
		
				self.sample_counts: <class 'dict'> {sample: num_seq_in_sample}
				self.data: <class 'dict'> contains all frequency data for each sample
						{
								sample1: {seq1: freq_sample1_seq1, ..., seqN: freq_sample1_seqN}
								sample2: {seq1: freq_sample2_seq1, ..., seqM: freq_sample2_seqM}
								...
								sampleK: {seq1: freq_sampleK_seq1, ..., seqJ: freq_sampleK_seqJ}
						}
				self.count_data: <class 'dict'> unused, but contains the distribution of counts for each sample
		'''
		if path.endswith('/'):
			self.path = path
		else:
			self.path = path+'/'

		self.samples = samples
		self.sequence_type = sequence_type
		self.current_path = os.getcwd()

		if wild_type is not None:
			self.wt = wild_type

		if mutations_dict is not None:
			self.mutations_dict = mutations_dict

		self.sample_counts = {line.split('\n')[0].split(',')[0]: float(line.split('\n')[0].split(',')[1]) 
								for line in open(f'{self.path}DMD_data/Sort_counts.csv')}		

		if self.sequence_type == 'sequences':
			self.data = self.get_round_dicts(f'{self.path}DMD_data/', 'Seq*.csv')
		elif self.sequence_type == 'mutations':
			self.data = self.get_round_dicts(f'{self.path}DMD_data/', 'Mu*.csv')
		elif self.sequence_type == 'CDRs':
			self.data = self.get_round_dicts(f'{self.path}DMD_data/', 'CD*.csv')

		self.count_data = {}
		for sample in self.samples:
			self.count_data[sample] = []
			for freq in self.data[sample].values():
				self.count_data[sample].append(float(freq*float(self.sample_counts[sample])))

	def get_seq_dict(self, path):
		'''
				get_seq_dict creates a dictionary of seq: freq from analyzed NGS files.
				input:
						path: path to NGS analyzed csvs.

				output:
						sequence dictionary
						{seq: freq}
		'''
		seq_dict = {}

		for line in open(path, 'r'):
			seq_dict[line.split(',')[0]] = float(line.split(',')[1].split('\n')[0])
		return seq_dict

	def get_round_dicts(self, path, file_type):
		'''
				calls get_seq_dict on all samples within a round and creates a dictionary of all of the individual dictionaries.
		'''
		os.chdir(path)

		muts = sorted([g for g in glob.glob(file_type)], key = lambda g: int(g.split('_')[1].split('.')[0]))
		print(muts)

		if len(muts) != len(self.samples):
			print('The number of samples does not match the number of sample files. Please update samples list.\n')
		if os.cpu_count() > 1:
			p = get_context('fork').Pool(os.cpu_count() - 1)
		else:
			p = get_context('fork').Pool()

		results = p.map(self.get_seq_dict, muts)
		p.close()
		p.join()
		
		r = {}
		for i in range(len(muts)):
			r[self.samples[i]] = results[i]

		os.chdir(self.current_path)
		return r

	def DMD(self, X_indecies = None, fillna = None, kernel = False, drop_clones = True):
		'''
		Performs standard DMD on self.data. 
		Custom indecies can be used instead of 0:len(t) - 1 and 1:len(t)
			X_indecies = [I1, I2, I3, I4]
				X = data[I1:I2]
				X' = data[I3:I4]
		'''
		if X_indecies is not None:
			X = pd.DataFrame(self.data)[list(self.data.keys())[X_indecies[0]:X_indecies[1]]]
			Xp = pd.DataFrame(self.data)[list(self.data.keys())[X_indecies[2]:X_indecies[3]]]
		
		else:
			X = pd.DataFrame(self.data)[list(self.data.keys())[0:len(self.data.keys()) - 1]]
			Xp = pd.DataFrame(self.data)[list(self.data.keys())[1:len(self.data.keys())]]
		
		if drop_clones:
			xclones = []
			for i in X.index:
				if not '*' in i:
					xclones.append(i)
			
			xpclones = []
			for i in Xp.index:
				if not '*' in i:
					xpclones.append(i)

			X = X.loc[xclones, :]
			Xp = Xp.loc[xpclones, :]
		
		#Condense to only common clones removing NAN values
		if fillna is None:
			clones = set(X.dropna().index).intersection(set(Xp.dropna().index))
			X = X.loc[clones, :]
			Xp = Xp.loc[clones, :]

		#replaces NAN with specified value
		else:
			X = X.fillna(fillna)
			Xp = Xp.fillna(fillna)
		
		if kernel:
			U, S, VcT = np.linalg.svd(X, full_matrices = False)
			S = np.diag(S)
			V = VcT.conj().T
			
			A = S @ VcT @ X.values.conj().T @ Xp.values @ V @ np.linalg.pinv(S)
			L, W = np.linalg.eig(A)
			pl = Xp @ V @ np.linalg.inv(S) @ L

			self.phi_L = pd.DataFrame(data = pl, index = X.index)

		else:
			U, S, VcT = np.linalg.svd(X, full_matrices = False)
			S = np.diag(S)
			V = VcT.conj().T

			A_tilde = U.conj().T @ Xp.values @ V @ np.linalg.inv(S)
			L, W = np.linalg.eig(A_tilde)

			pl = Xp.values @ V @ np.linalg.inv(S) @ L

			self.phi_L = pd.DataFrame(data = pl, index = X.index)
		
	def get_one_counts(self):
		self.one_counts = {}		
		for sample in self.samples:
			for freq in self.data[sample].values():
				if freq * self.sample_counts[sample] < 1.5:
					try:
						self.one_counts[sample] += 1
					except:
						self.one_counts[sample] = 1

		for s, c in self.one_counts.items():
			self.one_counts[s] = round(c / self.sample_counts[s], 3)
			print(f'{s} \t {round(self.one_counts[s] * 100, 1)}%')
	
	def plot_count_histograms(self, num_rows, num_cols, font_size = 10, figsize = (13,9)):
		fig = plt.figure(figsize = figsize)
		for i, sample in enumerate(self.samples):
			plt.rcParams.update({'font.size': font_size})
			ax = fig.add_subplot(num_rows, num_cols, i + 1)
			data = [self.data[sample][seq] * self.sample_counts[sample] for seq in self.data[sample].keys()] 
			ax.hist(data, bins=range(int(min(data)), 100, 1), label = f'Mean: {np.mean(data)}' + '\n' + f'Median: {np.median(data)}')
			ax.set_xlabel('Sequence Counts')
			ax.set_ylabel('Count')
			ax.set_title(f'{sample}')
			ax.legend()

# Simple functions for file management without need for making NGS_sample object.
# Unzip doesn't seem to work when put through a loop of os.chdir(f'Sample_{i}')-->unzip()-->os.chdir('../')
def unzip():
	for f in glob.glob('*.gz'):
		os.system(f'gunzip {f}')
		print('Unziping of {} complete.'.format(f))

def convert(y_n):
		'''
		converts the Sample_x_merged.fastq file to a .fasta file
		if a .fasta file is found in the directory it is removed.
		'''
		create_fasta = True
		fastas = glob.glob('*.fasta')
		if len(fastas) > 0:
			for a in glob.glob('*.fasta'):
				if y_n.upper() == 'Y':
					os.system('rm {}'.format(a))
					create_fasta = True
				else:
					create_fasta = False
					print('fasta already exists, not overwriting.')
		if create_fasta:
			for f in glob.glob('*_merged*.fastq'):
				os.system(f"paste - - - - < {f} | cut -f 1,2 | sed 's/^@/>/' | tr '\t' '\n' > {f.split('.fastq')[0]}.fasta")
				print('{} successfully converted to .fasta'.format(f))

def bam_convert(sample):
	'''
	converts pacbio .bam files into fastq files using bedtools.
	'''
	files = glob.glob('*.bam')
	for file in sorted(files, key = lambda f: int(f.split('.bc')[1].split('_')[0])):
		os.system(f'bedtools bamtofastq -i {file} -fq Sample_{sample}.fastq')
		print(f'{file} successfully converted to .fastq')

def pacbio_convert():
		'''
		converts the Sample_x_merged.fastq file to a .fasta file
		if a .fasta file is found in the directory it is removed.
		'''
		for a in glob.glob('*.fasta'):
			os.system('rm {}'.format(a))
		for f in glob.glob('*.fastq'):
			#os.system(f"cat {f} | paste - - - - |cut -f 1, 2| sed 's/@/>/'g | tr -s '/t' '/n' > {f.split('.fastq')[0]}.fasta")
			#os.system(f"sudo sed -n '1~4s/^@/>/p;2~4p' {f} > {f.split('.fastq')[0]}.fasta")
			os.system(f"paste - - - - < {f} | cut -f 1,2 | sed 's/^@/>/' | tr '\t' '\n' > {f.split('.fastq')[0]}.fasta")
			print('{} successfully converted to .fasta'.format(f))

def make_pacbio_sample_folders():
	''' moves all files into new separate sample folders'''
	pbi_files = glob.glob('*.pbi')
	bam_files = glob.glob('*.bam')
	xml_files = glob.glob('*.xml')

	assert len(pbi_files) == len(bam_files) and len(bam_files) == len(xml_files)

	for i in range(len(bam_files)):
		os.system(f'mkdir Sample_{i+1}')
		os.system(f"mv {sorted(pbi_files, key = lambda f: int(f.split('.bc')[1].split('_')[0]))[i]} Sample_{i+1}")
		os.system(f"mv {sorted(bam_files, key = lambda f: int(f.split('.bc')[1].split('_')[0]))[i]} Sample_{i+1}")
		os.system(f"mv {sorted(xml_files, key = lambda f: int(f.split('.bc')[1].split('_')[0]))[i]} Sample_{i+1}")

def fix_sample_folders(demux_stats_file):
	'''
	Please update demux stats file to be in the correct order.
	This function will change file folders from barcode ID to given sample ID. This uses a file called DemuxStats_{sample_ID}.csv The format is as follows:
		Project		Sample_ID						Description		Barcode					#Reads		%Reads		%Perfect Index Reads	% One Mismatch Index Reads
		5858-MS		5858-MS-1_CGAGAGTT-ACTATCTG		Ag_pos			CGAGAGTT-ACTATCTG		3355694		15.7		93.51					6.49
		5858-MS		5858-MS-1_CGAGAGTT-ATCGTACG		Input			CGAGAGTT-ATCGTACG		3240806		15.16		93.99					6.01

	An example file name is 5858-MS-1_CGAGAGTT-ATCGTACG_S1_R1_001.fastq. The input to this function is unique characters before the 
	to sort based on sample number.

	run this function in the directory containing the .fastq files.

	input:
		demux_stats_file: <str> file location of demux_stats_file. 
			ex:
				'../DemuxStats_sample_ID.csv'

	'''
	try:
		barcode_key = pd.read_csv(demux_stats_file)
	except AssertionError:
		print('Insure Demux Stats file is in directory.')
		sys.exit()
	
	sample_id = {}
	for i, d in enumerate(barcode_key.loc[:, 'Sample_ID']):
		try:
			os.system(f'mkdir Sample_{i + 1}')
			if f'{d}_S{i + 1}_R1_001.fastq' in os.listdir():
				os.system(f"mv {d}_S{i + 1}_R1_001.fastq ./Sample_{i + 1}/")
				os.system(f"mv {d}_S{i + 1}_R2_001.fastq ./Sample_{i + 1}/")
			elif f'{d}_S{i + 1}_R1_001.fastq.gz' in os.listdir():
				os.system(f"mv {d}_S{i + 1}_R1_001.fastq.gz ./Sample_{i + 1}/")
				os.system(f"mv {d}_S{i + 1}_R2_001.fastq.gz ./Sample_{i + 1}/")
		except:
			print('You might need to fix demux_stats_file to fix the order.')

def generate_round_sample_counts(sample_dict, Round):
	'''
	Input:
		sample_dict: <dict> {sample_id: sample_name} dictionary of ids to sample names
		Round: <int> round that data is from.
	'''
	csv = open(f'R{Round}_sort_counts.csv', 'w')
	
	for sample_id, sample_name in sample_dict.items():
		csv.write(f"{sample_name}, {int(open(f'Count_{sample_id}.csv').read())}\n")
	
	csv.close()
	
	try:
		os.chdir(f'R{Round}')
		os.chdir('../')
	except:
		os.system(f'mkdir R{Round}')
	
	for i in sample_dict.keys():
		files = glob.glob(f'*_{i}.csv')
		for file in files:
			os.system(f'mv {file} R{Round}')
	
	os.system(f'mv R{Round}_sort_counts.csv R{Round}')

def generate_time_series_counts(sample_dict):
	'''
	Input:
		sample_dict: <dict> {sample_id: sample_name} dictionary of ids to sample names
	'''
	csv = open(f'Sort_counts.csv', 'w')
	
	for sample_id, sample_name in sample_dict.items():
		csv.write(f"{sample_name}, {int(open(f'Count_{sample_id}.csv').read())}\n")

	csv.close()
	
	try:
		os.chdir(f'DMD_data')
		os.chdir('../')
	except:
		os.system(f'mkdir DMD_data')
	
	for i in sample_dict.keys():
		files = glob.glob(f'*_{i}.csv')
		for file in files:
			os.system(f'mv {file} DMD_data')
	
	os.system(f'mv Sort_counts.csv DMD_data')

#Useful functions to interact with NGS_round_data objects.
def common_clones(ngs_sample_dict, threshold = 0, data_type = 'clone'):
	'''
	returns a set of all clones common to all NGS_round_data in ngs_list and in all samples of sample_list
	Input: 
		ngs_sample_dict: <class 'dict'>
		nsd = {r1: ['sample1', 'sample2', ... ,'sampleN'], ..., rN: ['sampleA', 'sampleB', ... ,'SampleZ']}
	'''
	if type(ngs_sample_dict) != dict:
		print('Input should be dictionary')
		return
	
	if data_type == 'clone':
		initial = True
		for ngs_round, sample_list in ngs_sample_dict.items():
			for s in sample_list:
				seq_list = []
				for seq, freq in ngs_round.data[s].items():
					if freq * ngs_round.sample_counts[s] > threshold:
						seq_list.append(seq)
				if initial:
					seqs = set(seq_list)
					initial = False
				else:
					seqs = seqs.intersection(set(seq_list))
		return seqs
	
	else:
		initial = True
		for ngs_round, sample_list in ngs_sample_dict.items():
			for s in sample_list:
				seq_list = []
				for seq, freq in ngs_round.framework_data[s].items():
					if freq * ngs_round.sample_counts[s] > threshold:
						seq_list.append(seq)
				if initial:
					seqs = set(seq_list)
					initial = False
				else:
					seqs = seqs.intersection(set(seq_list))
		return seqs

def all_clones_count(ngs, sample, threshold):
	'''
	returns a list of all sequences in a given NGS_round_data sample with a count greater than or equal to the count
	'''
	seqs = []
	for seq, freq in ngs.data[sample].items():
		if freq * ngs.sample_counts[sample] >= threshold:
			seqs.append(seq)
	return seqs

def mut_analysis(ngs_round, mutation_str_length, num_mutations, er = None, er_order = None, sample_requirements = None, indecies = None, clone_set = None):
	'''
	Performs mutational analysis on a non NGS_round_data class. Can be useful for when enrichment ratios are non-standard or for multi concentration NGS data.

	Code should be optimized. Specifically creating mut_seqs_dict and wt_seqs_dict. Could probably alter code to find which positions are different from wt
	and only search mutation sets with the same differences. multiprocessing might also be added.
	'''
	mut_loc_dict  = {i: loc for i, loc in enumerate(list(itertools.combinations(list(range(mutation_str_length)), num_mutations)))}

	#create possible mutations dictionary
	mutation_number = 0
	mutations_dictionary = dict()
	for loc in tqdm.tqdm_notebook(mut_loc_dict.values()):
		mutations_dictionary_subkey = ''
		for i in loc:
			mutations_dictionary_subkey += str(i) + ','
		key_list = [i for i in mutations_dictionary_subkey.split(',') if i != '']
		product_list = [ngs_round.mutations_dict[int(j)] for j in key_list]
		possible_mutations = list(itertools.product(*product_list))
		for m in possible_mutations:
			mutations_dictionary[mutation_number] = {mutations_dictionary_subkey: m}
			mutation_number += 1
	
	if sample_requirements is not None or clone_set is not None:
		#calculate mut_seq_dict and wt_seq_dict  slowest part currently
		mut_seqs = {m: [] for m in mutations_dictionary.keys()}
		wt_seqs = {m: [] for m in mutations_dictionary.keys()}

		if clone_set is not None:
			seqs = clone_set
		else:
			seqs = common_clones({ngs_round: sample_requirements})

		for seq in tqdm.tqdm_notebook(seqs):
			for m, mut_dict in mutations_dictionary.items():
				for loc, muts in mut_dict.items():
					mut_count = 0
					wt_count = 0
					for i, l in enumerate(loc.split(',')):
						if l != '':
							if seq[int(l)] == muts[i]:
								mut_count += 1
							if seq[int(l)] == ngs_round.wt[int(l)]:
								wt_count += 1
					if mut_count == num_mutations:
						mut_seqs[m].append(seq)
					if wt_count == num_mutations:
						wt_seqs[m].append(seq)
	else:
		return mutations_dictionary

	if indecies is not None:
		#calculate spearman info
		spearman_coefficients = dict()
		for index in tqdm.tqdm_notebook(indecies):
			spearman_coefficients[er_order[index]] = dict()
			for m in tqdm.tqdm_notebook(mutations_dictionary.keys()): 
				data = {}
				for seq in mut_seqs[m]:
					try:
						data[er[seq][index]][0] += 1
					except:
						data[er[seq][index]] = [1, 0]
				for seq in wt_seqs[m]:
					try:
						data[er[seq][index]][1] += 1
					except:
						data[er[seq][index]] = [0, 1]

				er_data = []
				freq_data = []
				for e, nums in data.items():
					er_data.append(e)
					freq_data.append(nums[0] / (sum(nums)))
				spearman_coefficients[er_order[index]][m] = list(spearmanr(er_data, freq_data))

		return mutations_dictionary, spearman_coefficients, mut_seqs, wt_seqs
	
	else:
		return mutations_dictionary, mut_seqs, wt_seqs

def hamming_distance(seq1, seq2):
	d = 0
	for i, aa in enumerate(seq1):
		try:
			if aa != seq2[i]:
				d += 1
		except:
			print('Sequences incompatible')
	return  d

def mut_to_cdr(wt_cdr_str, mut_pos, mut_str):
	'''
	converts between mut_str and cdr_str.
	'''
	cdr_str = ''
	for i, aa in enumerate(wt_cdr_str):
		if i in mut_pos:
			cdr_str += mut_str[mut_pos.index(i)]
		else:
			cdr_str += aa
	
	return cdr_str

def cdr_to_mut(cdr_str, mut_pos):
	mutation_str = ''
	for i, aa in enumerate(cdr_str):
		if i in mut_pos:
			mutation_str += aa
	return mutation_str	

def mut_to_seq(wt, mut_pos, mut_str):
    '''
    converts mut_str to full sequence
    '''
    seq = ''
    for i, aa in enumerate(wt):
        if i in mut_pos:
            seq += mut_str[mut_pos.index(i)]
        else:
            seq += aa
    return seq
	
def seq_to_mut(seq, cdr_regions, mut_pos):
	mutation_str = ''
	cdr_str = ''
	for cdr in cdr_regions:
		cdr_str += seq[cdr[0]:cdr[1]]
		cdr_str += '_'

	for i, aa in enumerate(cdr_str):
		if i in mut_pos:
			mutation_str += aa
	
	return mutation_str

if __name__ == '__main__':
	pass

else:
	print('NGS package imported')
