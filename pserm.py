# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 11:33:12 2022

@author: mduncans@umich.edu
"""

import os
import time
import tqdm
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import get_context, Manager

#will make use of common_clones function
from ngs import common_clones

class ngs_analysis(object):
    '''
    create_data_matrix takes the two replicates of data and a list of rounds
    that sequences must be common to to create an averaged dataset of all the 
    sequences to consider and their frequency in each round.

    '''

    def __init__(self, replicates, common_rounds, clone_set = None, mutated_positions = None):
        if clone_set is not None:
            self.clone_set = clone_set
        else:
            cc_dict = {}
            for replicate in replicates:
                cc_dict[replicate] = common_rounds
            self.clone_set = common_clones(cc_dict)
        
        self.library = {i: [aa] for i, aa in enumerate(replicates[0].wt)}
        for i, aa_sampled in replicates[0].mutations_dict.items():
            for aa in aa_sampled:
                self.library[i].append(aa)
    
        self.wt = replicates[0].wt
        trimmed_clones = []
        for c in self.clone_set:
            if '*' not in c:
                trimmed_clones.append(c)
        
        self.clone_set = trimmed_clones
        self.replicates = replicates
        self.common_rounds = common_rounds
        self.AA_order = ['F', 'W', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        self.mutated_positions = mutated_positions
        self.Counts = {}
        self.PPM = {}
        self.PSSM = {}
        self.PSERM = {}
        
    # DATA MATRIX 
    def generate_D(self):

        df = []
        for r in self.replicates:
            d = pd.DataFrame(r.data).fillna(0) * list(r.sample_counts.values())
            df.append(d.loc[:, self.common_rounds])
        

        #find total counts of each sample to return to frequency
        for d in df:
            try:
                D = D.add(d, fill_value = 0)
            except:
                D = d.fillna(0)

        D = D.fillna(0)
        total_counts = D.loc[self.clone_set, :].sum(axis = 0)
        D /= total_counts

        #remove clones that were not considered for analysis
        self.D = D.loc[self.clone_set, :].fillna(0)
        self.total_counts = total_counts
        self.samples = self.D.columns
        self.counts = self.D * self.total_counts

    def generate_MSA(self):
        try:
            self.D
        except:
            self.generate_D()

        seqs_dict = {}
        
        for i, s in enumerate(self.D.index):
            seqs_dict[i] = []
            for aa in s:
                seqs_dict[i].append(aa)
            seqs_dict[i].append(s)
        
        msa_columns = [i for i in range(len(self.wt))]
        msa_columns.append('seq')

        self.MSA = pd.DataFrame.from_dict(seqs_dict, orient = 'index', columns = msa_columns)

    def generate_sequence_sets(self):
        try:
            self.D
        except:
            print('Generate datamatrix D first')
            return

        self.seq_sets = {f'{p}_{aa}': [] for p, aa_sample in self.library.items() for aa in aa_sample}
        for seq in tqdm.tqdm(self.D.index):
            for i, aa in enumerate(seq):
                self.seq_sets[f'{i}_{aa}'].append(seq)

    # PSSM Functions
    def generate_PSSM(self, sample, pseudocount = 1, background_prob = 1/20, show_plots = False, ppm_only = False, name = None, clone_set = None, excluded_sample = None, use_freq = True):
        '''
        generate_PSSM generates several different matrices
            1. self.Counts[name] is a dataframe containing the counts of each amino acid at each sampled position without pseudocounts.
            2. self.PPM[name] is a dataframe containing the relative frequencies of each amino acid at each sampled position with pseudocounts.
            3. self.PSSM[name] is a log2 transform of self.PPM[name] normalized to the background_prob
        
        sample: string, column name of self.D matrix, represents the sample of deep sequencing data to generate the above matrices for.
        '''
        try: 
            self.D
        except:
            print('Generate datam matrix (D) first')
        
        try:
            self.MSA
        except:
            print('Generating MSA.')
            self.generate_MSA()

        if name is None:
            name = sample

        try:
            self.PSSM[name]
        except:
            Count_dict = {pos: {AA: 0 for AA in self.AA_order} for pos in range(len(self.D.iloc[0, :].name))}
            Count_dict = pd.DataFrame.from_dict(Count_dict)
            
            if clone_set is None:
                if excluded_sample is not None:
                    clone_set = self.D.loc[(self.D[sample] > 0) & (self.D[excluded_sample] == 0)].index
                else:
                    clone_set = self.D.loc[self.D[sample] > 0].index
                use_D = True
            
            else:
                if use_freq:
                    use_D = True
                else:
                    use_D = False
            
            if use_D:
                for i, aa_list in tqdm.tqdm(self.library.items()):
                    for aa in aa_list:
                        seqs = list(set(self.MSA[self.MSA[i] == aa].loc[:, 'seq']).intersection(clone_set))
                        Count_dict.loc[aa, i] += self.counts.loc[seqs, sample].sum()

            else: 
                for i, aa_list in self.library.items():
                    for aa in aa_list:
                        total_num_seqs = len(self.MSA[self.MSA[i] == aa].loc[:, 'seq'])
                        Count_dict.loc[aa, i] += total_num_seqs

            self.Counts[name] = Count_dict.fillna(0).reindex(self.AA_order)
            
            N = self.Counts[name][0].sum() #each column will sum to same value.
            if type(pseudocount) == str:
                if pseudocount.lower() == 'proportional':
                    pseudocount = np.sqrt(N) 
            try:
                self.PPM[name] = (self.Counts[name] + pseudocount * background_prob) / (N + pseudocount)
            except:
                self.PPM = {name: (self.Counts[name] + pseudocount * background_prob) / (N + pseudocount)}
            
            self.PSSM[name] = np.log2(self.PPM[name] / background_prob).reindex(self.AA_order)        

    def generate_fixed_resi_PSSM(self, sample, pos, aa, pseudocount = 1, background_prob = 1/20, show_plots = False, ppm_only = False, excluded_sample = None):
        try:
            self.D
        except:
            print('Generate Data matrix first')
            return
        try:
            self.seq_sets
        except:
            self.generate_sequence_sets()
        
        if excluded_sample is not None:
            clone_set = set(self.D.loc[(self.D[sample] > 0) & (self.D[excluded_sample] == 0)].index)
            clone_set = clone_set.intersection(set(self.seq_sets[f'{pos}_{aa}']))
        else:
            clone_set = set(self.D.loc[(self.D[sample] > 0)].index)
            clone_set = clone_set.intersection(set(self.seq_sets[f'{pos}_{aa}']))

        self.generate_PSSM(sample, pseudocount = pseudocount, background_prob = background_prob, show_plots = show_plots, 
                            ppm_only = ppm_only, name = f'{pos}_{aa}_{sample}', clone_set = clone_set, use_freq = True)

    def generate_PSERM(self, In_sample, Out_sample, name = None):
        In_check = False
        Out_check = False
        if not In_sample in self.PSSM.keys():
            print(f'Please compute {In_sample} PSSM First')
        else:
            In_check = True
        
        if not Out_sample in self.PSSM.keys():
            print(f'Please compute {Out_sample} PSSM First')
        else:
            Out_check = True

        if In_check and Out_check:
            if name is None:
                self.PSERM[Out_sample] = self.PSSM[Out_sample] - self.PSSM[In_sample]
            else:
                self.PSERM[name] = self.PSSM[Out_sample] - self.PSSM[In_sample]

    def load_pssm(self, name, fname, excel_or_csv = 'csv'):
        if excel_or_csv == 'csv':
            try:
                self.PSSM[name] = pd.read_csv(f'{fname}', index_col = 0)
                self.PSSM[name].columns = self.PSSM[name].columns.astype('int')
            except:
                self.PSSM = {name: pd.read_csv(f'{fname}', index_col = 0)}
                self.PSSM[name].columns = self.PSSM[name].columns.astype('int')
        else:
            try:
                self.PSSM[name] = pd.read_excel(f'{fname}', index_col = 0)
                self.PSSM[name].columns = self.PSSM[name].columns.astype('int')
            except:
                self.PSSM = {name: pd.read_excel(f'{fname}', index_col = 0)}
                self.PSSM[name].columns = self.PSSM[name].columns.astype('int')
    
    def load_ppm(self, name, fname, excel_or_csv = 'csv'):
        if excel_or_csv == 'csv':
            try:
                self.PPM[name] = pd.read_csv(f'{fname}', index_col = 0)
                self.PPM[name].columns = self.PPM[name].columns.astype('int')
            except:
                self.PPM = {name: pd.read_csv(f'{fname}', index_col = 0)}
                self.PPM[name].columns = self.PPM[name].columns.astype('int')
        else:
            try:
                self.PPM[name] = pd.read_excel(f'{fname}', index_col = 0)
                self.PPM[name].columns = self.PPM[name].columns.astype('int')
            except:
                self.PPM = {name: pd.read_excel(f'{fname}', index_col = 0)}
                self.PPM[name].columns = self.PPM[name].columns.astype('int')

    def load_multi_mut_ppm(self, name, num_muts, fname, excel_or_csv = 'csv'):
        if excel_or_csv == 'csv':
            try:
                self.multi_muts_freq[num_muts][name] = pd.read_csv(f'{fname}', index_col = 0)
            except:
                try:
                    self.multi_muts_freq[num_muts] = {name: pd.read_csv(f'{fname}', index_col = 0)}
                except:
                    self.multi_muts_freq = {num_muts: {name: pd.read_csv(f'{fname}', index_col = 0)}}
        else:
            try:
                self.multi_muts_freq[num_muts][name] = pd.read_excel(f'{fname}', index_col = 0)
            except:
                try:
                    self.multi_muts_freq[num_muts] = {name: pd.read_excel(f'{fname}', index_col = 0)}
                except:
                    self.multi_muts_freq = {num_muts: {name: pd.read_excel(f'{fname}', index_col = 0)}}
                    
    #Create a load_scores function to load in precomputed scores
    
    # INFORMATION FUNCTIONS
    def get_H_pos(self, name):
        try:
            self.PPM[name]
        except:
            print('Create PPM first')
            return 
        try:
            self.H[name] = []
        except:
            self.H = {name: []}

        for i in self.PPM[name].columns:
            h = 0
            for aa in self.library[i]:
                h -= self.PPM[name].loc[aa, i] * np.log2(self.PPM[name].loc[aa, i])
            self.H[name].append((h)/np.log2(len(self.library[i])))

    def get_JS_divergence_pos(self, pos, aa, sample, name = None):
        #Make sure PSERM' and PSERM exist
        try:
            PSERMp = self.PSERM[f'{pos}_{aa}_{sample}']
        except:
            print(f'Create {pos}_{aa}_{sample} PSERM First.')
            return

        try:
            PSERM = self.PSERM[f'{sample}']
        except:
            print(f'Create {sample} PSERM First.')
            return
        
        #Create dictionary/list of positional JS divergences
        if name is None:
            name = f'{pos}_{aa}'
        try:
            self.JS[sample][name] = []
        except:    
            try:
                self.JS[sample] = {name: []}
            except:
                self.JS = {sample: {name: []}}

        for i, aas in self.library.items():
            p = np.array(PSERMp.loc[aas, i])
            q = np.array(PSERM.loc[aas, i])
            
            #Botlzmann like normalization
            p = 2**p / sum(2**p)
            q = 2**q / sum(2**q)
            m = 1/2 * (p + q)

            self.JS[sample][name].append(1/2 * self.KL_divergence(p, m) + self.KL_divergence(q, m))

    def get_all_mut_JS_divergence_matrix(self, sample):
        ''' Still needs to be coded... '''
        try:
            self.JS_data
        except:
            self.JS_data = {}

        for pos, aa_list in self.library.items():
            for aa in aa_list:
                try: 
                    self.PSERM[f'{pos}_{aa}_{sample}']
                except:
                    print('Please compute fixed residue PSERMs first')
                    return

        for pos, aa_list in self.library.items():
            for aa in aa_list:
                self.get_JS_divergence_pos(pos, aa, sample)
        if self.mutated_positions is not None:
            self.JS_data[sample] = pd.DataFrame.from_dict(self.JS[sample], orient = 'index', columns = self.mutated_positions)
        else:
            self.JS_data[sample] = pd.DataFrame.from_dict(self.JS[sample], orient = 'index')

    def KL_divergence(self, p, q):
        kld = 0

        for i, pi in enumerate(p):
            kld += pi * np.log2(pi / q[i])

        return kld

    # PLOTTING FUNCTION
    def plot_matrix(self, matrix_type, data, sample = None):

        fig, ax = plt.subplots(figsize = (7.2, 6))

        if matrix_type == 'Count':
            cax = ax.matshow(np.log10(data), aspect = 'auto')
            fig.colorbar(cax, ax = ax, label = r'$\log_{10}(C_{AA, pos})$')
        
        elif matrix_type == 'PPM':
            cax = ax.matshow(data, aspect = 'auto')
            fig.colorbar(cax, ax = ax, label = r'$\log_{10}(\frac{C_{AA, pos} + \sqrt{N}/20}{N + \sqrt{N}})$')
        
        elif matrix_type == 'PSSM':
            cax = ax.matshow(data, aspect = 'auto')
            fig.colorbar(cax, ax = ax, label = r'$\log_{2}(\frac{PPM_{AA, pos}}{p_aa})$')

        elif matrix_type == 'PSERM':
            cax = ax.matshow(data, aspect = 'auto')
            fig.colorbar(cax, ax = ax, label = f'{sample} ' + r'PSERM')

        ax.xaxis.set_ticks(list(range(data.shape[1])))
        ax.xaxis.set_ticklabels([i for i in range(1, data.shape[1] + 1)])

        ax.yaxis.set_ticks(list(range(20)))
        ax.yaxis.set_ticklabels([AA for AA in data.index])

        plt.show()

    # SCORING FUNCTIONS
    def score_all_samples_seq(self, inputs):
        #scores the input sequence with all PSERMs and stores them in the input dictionary.
        dataset, seq, method = inputs
        
        if method == 'PSERM':
            dataset[seq] = [score_seq(self.PSERM[sample], seq) for sample in self.samples if sample in self.PSERM.keys()]
        
        elif method == 'PSSM':
            dataset[seq] = [score_seq(self.PSSM[sample], seq) for sample in self.samples if sample in self.PSSM.keys()]

    def score_all_samples_seq_uncertainty(self, inputs):
        dataset, seq, In_sample = inputs

        dataset[seq] = [uncertainty_seq(self.Counts[In_sample], self.Counts[sample], seq) for sample in self.samples if sample in self.PSERM.keys()]

    def score_all_clones(self):
        scores = {}
        for s in self.clone_set:
            self.score_all_samples_seq([scores, s])

        self.scores = pd.DataFrame.from_dict(scores, orient = 'index', columns = [f'{s} PSERM Score' for s in self.PSERM.keys()])

    def score_all_clones_smart(self):
        scores = {}
        
    def score_all_clones_mp(self, method = 'PSERM'):
        ''' Scores all clones in self.D.index with all PSERM -- usually should use before context dependent PSERMs are made'''
        try:
            self.D
        except:
            print('Please compute data matrix, D, first.')

        manager = Manager()
        shared_dict = manager.dict()

        if os.cpu_count() > 1:
            p = get_context('fork').Pool(os.cpu_count() - 1)
        else:
            p = get_context('fork').Pool(1)

        
        t1 = time.time()
        result = p.map(self.score_all_samples_seq, [[shared_dict, s, method] for s in self.D.index])
        p.close()
        p.join()

        t2 = time.time()
        if method == 'PSERM':
            method_keys = self.PSERM.keys()
        elif method == 'PSSM':
            method_keys = self.PSSM.keys()

        #If self.scores exists concat new results with existing scores df
        try:    
           self.scores
           self.scores = pd.concat(
                [self.scores, 
                 pd.DataFrame.from_dict(shared_dict, orient = 'index', columns = [f'{s} {method} Score' for s in method_keys])],
                 axis = 1)
        #Else create new scores dataframe with results
        except:
            self.scores = pd.DataFrame.from_dict(shared_dict, orient = 'index', columns = [f'{s} {method} Score' for s in method_keys])
            

        print(f"Done in {round((t2 - t1)/60, 2)} minutes")

    def score_all_clones_uncertainty_mp(self, In_sample = 'Input'):
        ''' Scores all clones in self.D.index with all PSERM -- usually should use before context dependent PSERMs are made'''
        try:
            self.D
        except:
            print('Please compute data matrix, D, first.')

        manager = Manager()
        shared_dict = manager.dict()

        if os.cpu_count() > 1:
            p = get_context('fork').Pool(os.cpu_count() - 1)
        else:
            p = get_context('fork').Pool(1)

        
        t1 = time.time()
        result = p.map(self.score_all_samples_seq_uncertainty, [[shared_dict, s, In_sample] for s in self.D.index])
        p.close()
        p.join()

        t2 = time.time()
        method_keys = self.PSERM.keys()
  
        #If self.scores exists concat new results with existing scores df
        try:    
           self.uncertainties
           self.uncertainties = pd.concat(
                [self.uncertainties, 
                 pd.DataFrame.from_dict(shared_dict, orient = 'index', columns = [f'{s} PSERM Uncertainty' for s in method_keys])],
                 axis = 1)
        #Else create new scores dataframe with results
        except:
            self.uncertainties = pd.DataFrame.from_dict(shared_dict, orient = 'index', columns = [f'{s} PSERM Uncertainty' for s in method_keys])
            

        print(f"Done in {round((t2 - t1)/60, 2)} minutes")

    def score_all_clones_epistatic_wt(self, sample):
        
        try:
            self.D
        except:
            print('Please compute data matrix first')
            return 
            
        try:
            self.context_scores_wt
        except:
            self.context_scores_wt = {sample : {} for sample in self.samples}
        
        t1 = time.time()
        for seq in self.D.index:
            for sample in self.samples:
                score = 0
                for i, aa in enumerate(seq):
                    if aa != self.wt[i]:
                        for j, aap in enumerate(seq):
                            if j != i:
                                score += self.PSERM[f'{i}_{aa}_{sample}'].loc[aap, j]
                    else:
                        for j, aap in enumerate(seq):
                            if j != i:
                                score += self.PSERM[sample].loc[aap, j]
                        
                score /= len(seq)
            
            self.context_scores_wt[sample][seq] = score       
        t2 = time.time()
        print(f'Done in {(t2 - t1)/60:.2e} minutes')

    def score_all_clones_epistatic(self, sample):

        '''
        I want to score a clone using fixed residue PSERMs
            #seq: ACDEFG
            score = 1/6 * ( d(C_2 | A_1) + d(D_3 | A_1) + ... + d(G_6 | A_1) + 
                            d(A_1 | C_2) + d(D_3 | C_2) + ... + d(G_6 | C_2) +
                            d(A_1 | D_3) + d(C_2 | D_3) + ... + d(G_6 | C_2) +
                            ...
                            d(A_1 | G_6) + d(C_2 | G_6) + ... + d(F_5 | G_6) )
        '''
        try:
            self.D
        except:
            print('Please compute data matrix first')
            return 
            
        try:
            self.context_scores
        except:
            self.context_scores = {sample : {} for sample in self.samples}
        
        for seq in self.D.index:
            score = 0
            for i, aa in enumerate(seq):
                for j, aap in enumerate(seq):
                    if j != i:
                        score += self.PSERM[f'{i}_{aa}_{sample}'].loc[aap, j]        
            score /= len(seq)
            
            self.context_scores[sample][seq] = score

    # MULTIMUTATION SCORING
    def generate_multi_mut_freq_data(self, samples, num_muts = 2):
        try:
            self.MSA
        except:
            self.generate_MSA()

        self.single_mutations = [f'{p}_{aa}' for p in self.library.keys() for aa in self.library[p]]

        multi_positions = list(itertools.combinations(self.library.keys(), num_muts))
        
        multi_positions_fixed = []
        for mp in multi_positions:
            pos_id = ''
            for p in mp:
                pos_id += f'{p}_'
            pos_id = pos_id[0:-1]
            multi_positions_fixed.append(pos_id)

        try:
            self.multi_positions[num_muts] = multi_positions_fixed
        except:
            self.multi_positions = {num_muts: multi_positions_fixed}

        all_multi_mutations = list(itertools.combinations(self.single_mutations, num_muts))

        multi_mutations = []
        for mm in all_multi_mutations:
            positions = []
            mut_id = ''
            
            for m in mm:
                positions.append(int(m[0]))
                mut_id += f'{m[-1]}_'
            mut_id = mut_id[0:-1]

            if len(set(positions)) == num_muts:            
                multi_mutations.append(mut_id)

        try:
            self.multi_mutations[num_muts] = list(set(multi_mutations))
        except:
            self.multi_mutations = {num_muts: list(set(multi_mutations))}

        for s in samples:
            try:
                self.multi_muts_freq[num_muts][s] = pd.DataFrame(
                    np.zeros(shape = (len(self.multi_mutations[num_muts]), len(self.multi_positions[num_muts]))), 
                    columns = self.multi_positions[num_muts], 
                    index = self.multi_mutations[num_muts])
            except:
                try:
                    self.multi_muts_freq[num_muts] = {
                        s: pd.DataFrame(
                           np.zeros(shape = (len(self.multi_mutations[num_muts]), len(self.multi_positions[num_muts]))), 
                           columns = self.multi_positions[num_muts], 
                           index = self.multi_mutations[num_muts])
                        }
                except:
                    self.multi_muts_freq = {num_muts: {
                            s: pd.DataFrame(
                               np.zeros(shape = (len(self.multi_mutations[num_muts]), len(self.multi_positions[num_muts]))), 
                               columns = self.multi_positions[num_muts], 
                               index = self.multi_mutations[num_muts])
                            }
                        }

        counts = self.D * self.total_counts
        for multi_m in tqdm.tqdm(self.multi_muts_freq[num_muts][s].index): #Here s was defined previously and is just used to get the shape. they should all be the same for given number of mutations.
            for multi_p in self.multi_muts_freq[num_muts][s].columns:
                positions = np.array(multi_p.split('_'), dtype=int)
                mutations = multi_m.split('_')

                mask = True
                for i, p in enumerate(positions):
                    mask = mask & (self.MSA[p] == mutations[i])

                seqs = self.MSA.loc[mask].loc[:, 'seq']

                for s in samples:
                    try:
                        self.multi_muts_freq[num_muts][s].loc[multi_m, multi_p] = counts.loc[seqs, s].sum()
                    except:
                        pass

        for s in self.multi_muts_freq[num_muts].keys():
            self.multi_muts_freq[num_muts][s] /= self.total_counts[s]

    # LDA ANALYSIS
    def generate_one_hot_encoded(self, clone_set, name):
        '''
        This will create a one hot encoded matrix of all clones observed in any of the clone_set.
        '''
        try:
            self.ohe
        except:
            self.ohe = {}

        try:
            self.ohe[name]
        except:
            features = []
            for p, muts in self.library.items():
                for aa in muts:
                    features.append(f'{p}_{aa}')
            ohe = np.zeros(shape = ( len(clone_set), len(features) ))

            for i, seq in tqdm.tqdm(enumerate(clone_set), total = len(clone_set)):
                for j, aa in enumerate(seq):
                    ohe[i, features.index(f'{j}_{aa}')] = 1

            self.ohe[name] = pd.DataFrame(data = ohe, columns = features)

    def remove_correlated_features(self, name, corr_cutoff = 0.8):
        try:
            self.ohe[name]
        except:
            print(f'Compute {name} OHE matrix first')

        # Create np array of correlations
        corr = abs(self.ohe[name].corr()).values

        # Get indices of upper triangle
        up_tri_ind = np.triu_indices(n = corr.shape[0], k = 1)
        # get indices that need to be dropped. 
        ind_to_drop = (up_tri_ind[0][corr[up_tri_ind] > corr_cutoff], up_tri_ind[1][corr[up_tri_ind] > corr_cutoff])

        # Iterate through indices and drop columns. 
        indices_kept = []
        indices_dropped = []
        for i1, i2 in zip(*ind_to_drop):
            indices_kept.append(i1)
            indices_dropped.append(i2)

        # Check to make sure no columns should be kept AND dropped... Need to update what to do if this isn't true.
        if len(set(indices_kept).intersection(set(indices_dropped))) == 0:
            columns_dropped = []
            for i in indices_dropped:
                columns_dropped.append(self.ohe[name].columns[i])
            try:
                self.dropped_columns[name] = columns_dropped
            except:
                self.dropped_columns = {name: columns_dropped}
            self.ohe[f'{name}_fixed'] = self.ohe[name].drop(columns_dropped, axis = 1)


# NGS ANALYSIS CLASS FOR LIBRARIES WITH DIFFERENT FRAMEWORKS/LENGTHS
class multi_fr_ngs_analysis(ngs_analysis):
    
    def __init__(self, replicates, common_rounds, clone_set = None, mutated_positions = None):
        if clone_set is not None:
            self.clone_set = clone_set
        else:
            cc_dict = {}
            for replicate in replicates:
                cc_dict[replicate] = common_rounds
            self.clone_set = common_clones(cc_dict)

        self.library = {fr: {i:[] for i,_ in enumerate(seq)} for fr, seq in replicates[0].wt.items()}
        for fr, mut_dict in replicates[0].mutations_dict.items():
            for i, aa_sampled in mut_dict.items():
                for aa in aa_sampled:
                    self.library[fr][i].append(aa)
        
        self.wt = replicates[0].wt
        trimmed_clones = []
        for c in self.clone_set:
            if '*' not in c:
                trimmed_clones.append(c)
        
        self.clone_set = trimmed_clones
        self.replicates = replicates
        self.common_rounds = common_rounds
        self.AA_order = ['F', 'W', 'Y', 'P', 'M', 'I', 'L', 'V', 'A', 'G', 'C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R']
        self.mutated_positions = mutated_positions
        self.Counts = {}
        self.PPM = {}
        self.PSSM = {}
        self.PSERM = {}    
    
    # DATA MATRIX 
    def generate_D(self):

        df = []
        for r in self.replicates:
            d = pd.DataFrame(r.data).fillna(0) * list(r.sample_counts.values())
            df.append(d.loc[:, self.common_rounds])
        

        #find total counts of each sample to return to frequency
        for d in df:
            try:
                D = D.add(d, fill_value = 0)
            except:
                D = d.fillna(0)

        D = D.fillna(0)
        total_counts = D.loc[self.clone_set, :].sum(axis = 0)
        D /= total_counts

        #remove clones that were not considered for analysis
        self.D = D.loc[self.clone_set, :].fillna(0)
        self.total_counts = total_counts
        self.samples = self.D.columns
        self.counts = self.D * self.total_counts
        
        FR = []
        for seq in self.D.index:
            for fr, wt_seq in self.wt.items():
                if len(seq) == len(wt_seq):
                    FR.append(fr)
                    break
        self.D.loc[:, 'FR'] = FR

    def generate_MSA(self):
        try:
            self.D
        except:
            self.generate_D()

        seqs_dict = {fr: {} for fr in self.wt.keys()}
        for fr in self.wt.keys():
            for i, s in enumerate(self.D[self.D['FR'] == fr].index):
                seqs_dict[fr][i] = []
                for aa in s:
                    seqs_dict[fr][i].append(aa)
                seqs_dict[fr][i].append(s)
        
        msa_columns = {fr: [i for i in range(len(seq))] for fr, seq in self.wt.items()}
        for fr in self.wt.keys():
            msa_columns[fr].append('seq')

        self.MSA = {fr: pd.DataFrame.from_dict(seqs_dict[fr], orient = 'index', columns = msa_columns[fr]) for fr in self.wt.keys()}

    def generate_PSSM(self, sample, FR, pseudocount = 1, background_prob = 1/20, show_plots = False, ppm_only = False, name = None, clone_set = None, excluded_sample = None, use_freq = True):
        '''
        generate_PSSM creates a Position Specific Scoring Matrix
        based on the input data for a given sample (column of data matrix)
        '''
        try: 
            self.D
        except:
            print('Generate datam matrix (D) first')
        
        try:
            self.MSA
        except:
            print('Generating MSA.')
            self.generate_MSA()

        if name is None:
            name = f'{sample}_{FR}'

        try:
            self.PSSM[name]
        except:
            Count_dict = {pos: {AA: 0 for AA in self.AA_order} for pos in range(len(self.wt[FR]))}
            Count_dict = pd.DataFrame.from_dict(Count_dict)
            
            if clone_set is None:
                if excluded_sample is not None:
                    clone_set = self.D.loc[(self.D[sample] > 0) & (self.D[excluded_sample] == 0) & (self.D['FR'] == FR)].index
                else:
                    clone_set = self.D.loc[(self.D[sample] > 0) & (self.D['FR'] == FR)].index
                use_D = True
            
            else:
                if use_freq:
                    use_D = True
                else:
                    use_D = False
            
            if use_D:
                for i, aa_list in tqdm.tqdm(self.library[FR].items()):
                    for aa in aa_list:
                        seqs = self.MSA[FR][self.MSA[FR][i] == aa].loc[:, 'seq']
                        Count_dict.loc[aa, i] += self.counts.loc[seqs, sample].sum()

            else: 
                for i, aa_list in self.library.items():
                    for aa in aa_list:
                        total_num_seqs = len(self.MSA[FR][self.MSA[FR][i] == aa].loc[:, 'seq'])
                        Count_dict.loc[aa, i] += total_num_seqs

            self.Counts[name] = Count_dict.fillna(0).reindex(self.AA_order)
            
            N = self.Counts[name][0].sum() #each column will sum to same value.
            if type(pseudocount) == str:
                if pseudocount.lower() == 'proportional':
                    pseudocount = np.sqrt(N) 
            
            self.PPM[name] = (self.Counts[name] + pseudocount * background_prob) / (N + pseudocount)
            
            self.PSSM[name] = np.log2(self.PPM[name] / background_prob).reindex(self.AA_order)

    def generate_PSERM(self, In_sample, Out_sample, FR, name = None):
        In_check = False
        Out_check = False
        if not f'{In_sample}_{FR}' in self.PSSM.keys():
            print(f'Please compute {In_sample}_{FR} PSSM First')
        else:
            In_check = True
        
        if not f'{Out_sample}_{FR}' in self.PSSM.keys():
            print(f'Please compute {Out_sample}_{FR} PSSM First')
        else:
            Out_check = True

        if In_check and Out_check:
            if name is None:
                self.PSERM[f'{Out_sample}_{FR}'] = self.PSSM[f'{Out_sample}_{FR}'] - self.PSSM[f'{In_sample}_{FR}']
            else:
                self.PSERM[name] = self.PSSM[f'{Out_sample}_{FR}'] - self.PSSM[f'{In_sample}_{FR}']

    def score_all_samples_fr_seq(self, inputs):
        #scores the input sequence with all PSERMs and stores them in the input dictionary.
        dataset, seq, fr, method = inputs
        
        if method == 'PSERM':
            dataset[seq] = [score_seq(self.PSERM[f'{sample}_{fr}'], seq) for sample in self.samples if f'{sample}_{fr}' in self.PSERM.keys()]
        
        elif method == 'PSSM':
            dataset[seq] = [score_seq(self.PSSM[f'{sample}_{fr}'], seq) for sample in self.samples if f'{sample}_{fr}' in self.PSSM.keys()]

    def score_all_clones_mp(self, method='PSERM'):
        ''' Scores all clones in self.D.index with all PSERM -- usually should use before context dependent PSERMs are made'''
        try:
            self.D
        except:
            print('Please compute data matrix, D, first.')

        manager = Manager()
        shared_dict = manager.dict()

        if os.cpu_count() > 1:
            p = get_context('fork').Pool(os.cpu_count() - 1)
        else:
            p = get_context('fork').Pool(1)

        
        t1 = time.time()
        p.map(self.score_all_samples_fr_seq, [[shared_dict, s, self.D.iloc[i, -1], method] for i, s in enumerate(self.D.index)])
        p.close()
        p.join()

        t2 = time.time()
        method_keys = []
        if method == 'PSERM':
            for s in self.PSERM.keys():
                key = s.split('_')[0]
                if key not in method_keys:
                    method_keys.append(key)
        elif method == 'PSSM':
            for s in self.PSSM.keys():
                key = s.split('_')[0]
                if key not in method_keys:
                    method_keys.append(key)

        #If self.scores exists concat new results with existing scores df
        try:    
           self.scores
           self.scores = pd.concat(
                [self.scores, 
                 pd.DataFrame.from_dict(shared_dict, orient = 'index', columns = [f'{s} {method} Score' for s in method_keys])],
                 axis = 1)
        #Else create new scores dataframe with results
        except:
            self.scores = pd.DataFrame.from_dict(shared_dict, orient = 'index', columns = [f'{s} {method} Score' for s in method_keys])
            

        print(f"Done in {round((t2 - t1)/60, 2)} minutes")        

# GENERAL FUNCTIONS
def score_seq(PSSM, seq, ref_seq = None, ref_locs = None):
    score = 0 
    for i, AA in enumerate(seq):
        if ref_locs is not None:
            if i in ref_locs:
                if ref_seq is not None:
                    score += PSSM.loc[AA, i] - PSSM.loc[ref_seq[i], i]
                else:
                    score += PSSM.loc[AA, i]
        else:
            if ref_seq is not None:
                score += PSSM.loc[AA, i] - PSSM.loc[ref_seq[i], i]
            else:
                score += PSSM.loc[AA, i]

    return score

def uncertainty_seq(In_count, Out_count, seq):
    score = 0 
    for i, AA in enumerate(seq):
        score += ((1 / Out_count.loc[AA, i] + 1 / In_count.loc[AA, i]) / np.log(2))**2
    return np.sqrt(score)

def generate_clone_set(ngs_round_data, samples, excluded_samples = None):
    clone_set = set()
    for r in samples:
        clone_set = clone_set.union(common_clones({ngs_round_data: [r]}))

    if excluded_samples is not None:
        for r in excluded_samples:
            clone_set = clone_set - set(common_clones({ngs_round_data: [r]}))

    clone_set_trimmed = []
    for seq in tqdm.tqdm(clone_set):
        if '*' not in seq:
            if 'X' not in seq:
                if len(seq) == len(ngs_round_data.wt):
                    correct_count = 0
                    for i, aa in enumerate(seq):
                        if aa == ngs_round_data.wt[i]:
                            correct_count += 1
                        elif aa in ngs_round_data.mutations_dict[i]:
                            correct_count += 1
                    if correct_count == len(ngs_round_data.wt):
                        clone_set_trimmed.append(seq)

    return clone_set_trimmed



    