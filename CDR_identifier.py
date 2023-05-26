#!/usr/bin/env python3
#author: mduncans
#180511 23:52

import re
import numpy as np

def cdr_finder(VH_seq = None, vh_loc = None, VL_seq = None, vl_loc = None, definition = 'kabat', which_CDRHs = [True, True, True], which_CDRLs = [True, True, True]):
	'''
	print cmd.get_fastastr(<selection name>) prints sequence to cmd in pymol
	
	input antibody HC or LC protein sequence and obtain dictionary of 
	cdrs and their starting locations. 
	Curently can identify cdr based on kabat or chothia definitions
	CDRH1 might be off. 
	'''
	cdr = {}

	if VL_seq is not None:
		CDRL1 = which_CDRLs[0]
		CDRL2 = which_CDRLs[1]
		CDRL3 = which_CDRLs[2]

		cdrl1_start = re.search(r'[A-Z]C[A-Z]', VL_seq).span()[1] - 1
		cdrl1_stop_seq = ['WYQ', 'WLQ', 'WFQ', 'WYL', r'[A-Z]W[A-Z]']
		for seq in cdrl1_stop_seq:
			cdrl1_stop = re.search(seq, VL_seq)
			try:
				if cdrl1_stop.span()[0] > cdrl1_start + 9:
					cdrl1_stop = cdrl1_stop.span()[0]
					break
			except AttributeError:
				continue
		if CDRL1:
			cdr['L1'] = (VL_seq[cdrl1_start:cdrl1_stop], (cdrl1_start, cdrl1_stop))

		if CDRL2:
			cdrl2_start = cdrl1_stop + 15
			cdrl2_stop = cdrl2_start + 7
			cdr['L2'] = (VL_seq[cdrl2_start:cdrl2_stop], (cdrl2_start, cdrl2_stop))
		if CDRL3:
			cdrl2_start = cdrl1_stop + 15
			cdrl2_stop = cdrl2_start + 7
			cdrl3_start = cdrl2_stop + 32
			cdrl3_stop = re.search(r'[FL][AG][A-Z]G', VL_seq).span()[0]
			cdr['L3'] = (VL_seq[cdrl3_start:cdrl3_stop], (cdrl3_start, cdrl3_stop))
	
	if VH_seq is not None:
		CDRH1, CDRH2, CDRH3 = which_CDRHs

		if definition == 'kabat':
			cdrh1_start_seq = re.findall(r'C[AKSTV][AFGISTVP][STY]', VH_seq)
			if len(cdrh1_start_seq) > 1:
				for seq in cdrh1_start_seq:
					index = VH_seq.index(seq)
					if index > 10:
						if index + 4 < 30:
							cdrh1_start = index + 4
			else:
				cdrh1_start = re.search(r'C[AKSTV][AFGISTVP][STY]', VH_seq).span()[1]

			cdrh1_stop_seq = re.findall(r'[A-Z]W[VIAFYWDL][KQR]', VH_seq)
			for seq in cdrh1_stop_seq:
				index = VH_seq.index(seq)
				if index > cdrh1_start:
					if index - cdrh1_start < 20:
						cdrh1_stop = index + 1
			if CDRH1:	
				cdr['H1'] = (VH_seq[cdrh1_start:cdrh1_stop], (cdrh1_start, cdrh1_stop))
		else:
			cdrh1_start_seq = re.findall(r'C[A-Z][A-Z][A-Z]', VH_seq)
			if len(cdrh1_start_seq) > 1:
				for seq in cdrh1_start_seq:
					index = VH_seq.index(seq)
					if index > 10:
						if index < 30:
							cdrh1_start = index + 4
			else:
				cdrh1_start = re.search(r'C[A-Z][A-Z][A-Z]', VH_seq).span()[1]
			
			cdrh1_stop_seq = re.findall(r'[A-Z]W[VIAFY]', VH_seq)
			for seq in cdrh1_stop_seq:
				index = VH_seq.index(seq)
				if index > cdrh1_start:
					if index - cdrh1_start < 20:
						cdrh1_stop = index -2
			if CDRH1:
				cdr['H1'] = (VH_seq[cdrh1_start:cdrh1_stop], (cdrh1_start, cdrh1_stop))
				
		if CDRH2:
			if definition == 'kabat':
				cdrh2_start = cdrh1_stop + 14
			else:
				cdrh2_start =  cdrh1_stop + 19
			cdrh2_stop_seq = re.findall(r'[KRT][LIVFTAPNGSDF][TSIAQERK]', VH_seq)
			cdr_lengths = []
			for seq in cdrh2_stop_seq:
				index = VH_seq.index(seq)
				if index > cdrh2_start:
					if index - 2 - cdrh2_start < 22:
						cdr_lengths.append(index - 2 - cdrh2_start)
					else:
						cdr_lengths.append(0)
				else:
					cdr_lengths.append(0)
			
			if VH_seq[VH_seq.index(cdrh2_stop_seq[cdr_lengths.index(max(cdr_lengths))]) - 2] in 'RK':
				cdrh2_stop = VH_seq.index(cdrh2_stop_seq[cdr_lengths.index(max(cdr_lengths))]) - 2
			elif VH_seq[VH_seq.index(cdrh2_stop_seq[cdr_lengths.index(max(cdr_lengths))])] in 'RKT':
				cdrh2_stop = VH_seq.index(cdrh2_stop_seq[cdr_lengths.index(max(cdr_lengths))])
			elif VH_seq[VH_seq.index(cdrh2_stop_seq[cdr_lengths.index(max(cdr_lengths))])] in 'KFQ':
				cdrh2_stop = VH_seq.index(cdrh2_stop_seq[cdr_lengths.index(max(cdr_lengths))]) - 4

			if definition == 'chothia':
				cdrh2_stop -7

			cdr['H2'] = (VH_seq[cdrh2_start:cdrh2_stop], (cdrh2_start, cdrh2_stop))

		if CDRH3:
			if definition == 'kabat':
				cdrh2_start = cdrh1_stop + 14
			else:
				cdrh2_start =  cdrh1_stop + 19
			cdrh2_stop_seq = re.findall(r'[KR][LIVFTAPGS][TSIAER]', VH_seq)
			for seq in cdrh2_stop_seq:
				try:
					index = VH_seq[cdrh2_start::].index(seq) + cdrh2_start
					if index > cdrh2_start:
						if index - cdrh2_start < 20:
							cdrh2_stop = index
							if definition == 'chothia':
								cdrh2_stop = cdrh2_stop - 7
				except:
					pass

			cdrh3_start_seq = re.findall(r'C[ABD-Z][A-Z]', VH_seq)
			if len(cdrh3_start_seq) > 1:
				if len(set(cdrh3_start_seq)) > 1:
					distances = []
					for seq in cdrh3_start_seq:
						distances.append(VH_seq.index(seq) - cdrh2_stop)
					index = 0 
					current_difference = np.inf
					
					for i, d in enumerate(distances):
						if definition == 'chothia':
							if d > 0:
								if abs(d - 39) < current_difference:
									current_difference = abs(d - 39)
									index = i
						else:
							if d > 0:
								if abs(d - 32) < current_difference:
									current_difference = abs(d - 32)
									index = i
					cdrh3_start = VH_seq.index(cdrh3_start_seq[index]) + 1
				else:
					for seq in set(cdrh3_start_seq):
						cdrh3_start = cdrh2_stop + 1 + VH_seq[cdrh2_stop::].index(seq)
			else:
				cdrh3_start = VH_seq.index(cdrh3_start_seq[0]) + 1
			cdrh3_stop_seq = re.findall(r'[WR][GIRAPS][A-Z]G', VH_seq)

			for seq in cdrh3_stop_seq:
				try:
					index = VH_seq[cdrh2_stop::].index(seq) + cdrh2_stop
					if index > cdrh2_stop:
						cdrh3_stop = index
				except: 
					pass
			cdr['H3'] = (VH_seq[cdrh3_start:cdrh3_stop], (cdrh3_start, cdrh3_stop))

	return cdr

def split_sequences(VH_seq = None, VL_seq = None):

	split_seq = {}

	if VH_seq is not None:
		if VL_seq is not None:
			cdr = cdr_finder(VH_seq = VH_seq, VL_seq = VL)

			split_seq['VH'] = []
			split_seq['VL'] = []

			split_seq['VH'].append(VH_seq[0:cdr['H1'][1][0] - 1])
			split_seq['VH'].append(cdr['H1'][0])
			split_seq['VH'].append(VH_seq[cdr['H1'][1][1]:cdr['H2'][1][0] - 1])
			split_seq['VH'].append(cdr['H2'][0])
			split_seq['VH'].append(VH_seq[cdr['H2'][1][1]:cdr['H3'][1][0] - 1])
			split_seq['VH'].append(cdr['H3'][0])
			split_seq['VH'].append(VH_seq[cdr['H3'][1][1]::])

			split_seq['VL'].append(VH_seq[0:cdr['L1'][1][0] - 1])
			split_seq['VL'].append(cdr['L1'][0])
			split_seq['VL'].append(VH_seq[cdr['L1'][1][1]:cdr['L2'][1][0] - 1])
			split_seq['VL'].append(cdr['L2'][0])
			split_seq['VL'].append(VH_seq[cdr['L2'][1][1]:cdr['L3'][1][0] - 1])
			split_seq['VL'].append(cdr['L3'][0])
			split_seq['VL'].append(VH_seq[cdr['L3'][1][1]::])

		else:
			cdr = cdr_finder(VH_seq = VH_seq)

			split_seq['VH'] = []

			split_seq['VH'].append(VH_seq[0:cdr['H1'][1][0] - 1])
			split_seq['VH'].append(cdr['H1'][0])
			split_seq['VH'].append(VH_seq[cdr['H1'][1][1]:cdr['H2'][1][0] - 1])
			split_seq['VH'].append(cdr['H2'][0])
			split_seq['VH'].append(VH_seq[cdr['H2'][1][1]:cdr['H3'][1][0] - 1])
			split_seq['VH'].append(cdr['H3'][0])
			split_seq['VH'].append(VH_seq[cdr['H3'][1][1]::])


	elif VL_seq is not None:
		cdr = cdr_finder(VL_seq = VL_seq)

		split_seq['VL'] = []

		split_seq['VL'].append(VH_seq[0:cdr['L1'][1][0] - 1])
		split_seq['VL'].append(cdr['L1'][0])
		split_seq['VL'].append(VH_seq[cdr['L1'][1][1]:cdr['L2'][1][0] - 1])
		split_seq['VL'].append(cdr['L2'][0])
		split_seq['VL'].append(VH_seq[cdr['L2'][1][1]:cdr['L3'][1][0] - 1])
		split_seq['VL'].append(cdr['L3'][0])
		split_seq['VL'].append(VH_seq[cdr['L3'][1][1]::])

	return split_seq
	

VH_N2 = 'QGQLVESGGGSVQAGGSLRLSCAASGIDSSSYCMGWFRQRPGKEREGVARINGLGGVKTAYADSVKDRFTISRDNAENTVYLQMNSLKPEDTAIYYCAAKFSPGYCGGSWSNFGYWGQGTQVTVSSH'
test = 'EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYGMSWVRQAPGKGLELVAYITDCGGATSYPSSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCASGDNWGQGTTVTVSS'
VH_1igt = 'EVKLQESGGGLVQPGGSLKLSCATSGFTFSDYYMYWVRQTPEKRLEWVAYISNGGGSTYYPDTVKGRFTISRDNAKNTLYLQMSRLKSEDTAMYYCARHGGYYAMDYWGQGTTVTVSSAKTTAPSVYPLAPVCGDTTGSSVTLGCLVKGYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPRGPTIKPCPPCKCPAPNLLGGPSVFIFPPKIKDVLMISLSPIVTCVVVDVSEDDPDVQISWFVNNVEVHTAQTQTHREDYNSTLRVVSALPIQHQDWMSGKEFKCKVNNKDLPAPIERTISKPKGSVRAPQVYVLPPPEEEMTKKQVTLTCMVTDFMPEDIYVEWTNNGKTELNYKNTEPVLDSDGSYFMYSKLRVEKKNWVERNSYSCSVVHEGLHNHHTTKSFSR'
VL_1igt = 'DIVLTQSPSSLSASLGDTITITCHASQNINVWLSWYQQKPGNIPKLLIYKASNLHTGVPSRFSGSGSGTGFTLTISSLQPEDIATYYCQQGQSYPLTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC'

VH_3i75 = 'EVKLEESGAELVRPGASVTLSCAASGYTFTDFEIHWVKQPPVGGLEWIGTLDPETGGTAYNQNFKGRATLTADKSSSTAYMELRSLTSEDSAVYYCTRWGKKFYYYGTSYAMDYWGQGTSVTVSSAKTTPPSVYPLAPGATNSMVTLGCLVKGYFPEPVTVTWNSGSLSGGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKIVPRD'
VL_3i75 = 'DIQMTQSPSSLSASLGGKVTITCQSSQDINKYIGWYQHKPGKGPRLLIHYTSILRPDIPSRFSGSGSGRDYSFSISNLEPEDTATYYCLQYDDLLLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSETDQDSKDSTYSMSSTLTLTKDEYERHNTYTCEATHKTSTSPIVKSFNRA'

VH_1yy8 = 'QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSAASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKS'
VL_1yy8 = 'DILLTQSPVILSVSPGERVSFSCRASQSIGTNIHWYQQRTNGSPRLLIKYASESISGIPSRFSGSGSGTDFTLSINSVESEDIADYYCQQNNNWPTTFGAGTKLELKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGA'

VH_3B9V = 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIGWVRRAPGKGEEWVASIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARWGGDGFYAMDYWGQGTLVTVS'
seq =   'TITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTLGQGTKVEIKRTSPNSASHSGSAPQTSSAPGSQEVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTA'
seq2 =  'TITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTSPNSASHSGSAPQTSSAPGSQEVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTA'
works = 'TITCRASQNVANAVSWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHSTTPPTFGQGTKVEIKRTSPNSASHSGSAPQTSSAPGSQEVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPAGGYTRYAASVKGRFTISADTSKNTA'

VH =  'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPAGGYTRYAASVKGRFTISADTSKNTA'
vh =  'EVQLVESGGGWVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPSNGYTRYAGSVKGRFTISADTSKNTA'
#vh2 = 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWGRQAPGKGLEWMARIAPANGTTGYAGSVKGRFTISADTSKDTA'


VL =  'TITCRASQNVANAVSWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHSTTPPTFGQGTKVEIKRT'
vl =  'TITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTLGQGTKVEIKRT'
#vl2 = 'TITRRASQDVDAAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHSTSPPTFGQGTKVEIKRT'

if __name__ == '__main__':
	kabat_cdr = cdr_finder(VH, VL, which_CDRHs = [False, True, False], which_CDRLs = [True, False, True])
	cdr = cdr_finder(vh, vl, which_CDRHs = [False, True, False], which_CDRLs = [True, False, True])
	#cdr2 = cdr_finder(vh2, vl2, which_CDRHs = [False, True, False], which_CDRLs = [True, False, True])
	print(kabat_cdr)
	print(cdr)
	#print(cdr2)
	#kabat_cdr = cdr_finder(VH_N2)

	print(VL[4:15])
	print(VL[69:78])
	print(VH[49:66])


