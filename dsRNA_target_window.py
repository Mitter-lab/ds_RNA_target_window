#Copyright Stephen Fletcher
import get_ref
from collections import Counter
import operator
import csv
import os
import argparse

"""Identifies optimal dsRNA target site of window size n in 
multiple alignment in fasta format.  

At each position in the window, scores:
	Count of most common nucelotide
	Bonus for all seqs containing common nucletide
	Bonus for runs of positions with all seqs containing common nucletide - compounds

Produces results as an ordered score for every window, from best to worst.

Outputs in .csv format
"""

"""
To DO : Include or not include gaps???
"""


def window_analysis(aln_file, out_file, win_size = 300, pos_bonus_mod = 5, run_bonus_mod = 10, twenty_one_bonus_mod = 100):
	win_score = 0
	max_nuc_count = 0

	prev_run = False
	run_bonus_init = run_bonus_mod
	nuc_pos = []

	win_pos = 0
	combined_scores={}
	aln_seq = get_ref.get_ref_f_strand(aln_file) #load seq from aln fasta file
	seq_len = len(aln_seq[aln_seq.keys()[0]])
	while win_pos < seq_len - win_size:
		win_score = 0
		prev_run = False
		nt_21_run = 0
		tot_pos_bonus = 0
		tot_run_bonus = 0
		tot_21_bonus = 0


		for pos in range(win_pos, win_pos+win_size):
			nuc_pos = []	
			for header, seq in aln_seq.iteritems():
				nuc_pos.append(seq[pos]) #make a list of nucs for that position
			nuc_counts =  dict(Counter(nuc_pos)) #count their occurences
			try:
				del nuc_counts['-'] #remove gaps   DO WE ACTUALLY NEED TO DO THIS
			except:
				pass
			try:
				max_nuc_count =  nuc_counts[max(nuc_counts.iteritems(), key=operator.itemgetter(1))[0]]
			except:
				max_nuc_count = 0
			win_score += max_nuc_count
			if max_nuc_count == len(aln_seq): 
				win_score += pos_bonus_mod
				tot_pos_bonus += pos_bonus_mod
				if prev_run == True:
					win_score+=run_bonus_mod
					tot_run_bonus +=run_bonus_mod
					nt_21_run +=1
					run_bonus_mod += run_bonus_init
					if nt_21_run == 20: 
						win_score += twenty_one_bonus_mod
						tot_21_bonus += twenty_one_bonus_mod
						nt_21_run = 0
				prev_run = True
			else:
				prev_run = False
				nt_21_run = 0
				run_bonus_mod = run_bonus_init
			

		combined_scores[win_pos] = win_score
		print "Win score = {0}\tPos_bonus_score = {1}\trun_bonus_score={2}\ttwenty_one_bonus_score={3} at position {4}".format(win_score, tot_pos_bonus, tot_run_bonus, tot_21_bonus, win_pos)	
		win_pos +=1	
	print "\ninput file = {0}".format(aln_file)
	print "output file = {0}".format(out_file)
	print "\nwindow size = {0}".format(win_size)
	print "\ncomplete conservation at position bonus = {0}".format(pos_bonus_mod)
	print "run of conserved positions bonus = {0}".format(run_bonus_mod)
	print '\npos with max score is {0}'.format(max(combined_scores.iteritems(), key=operator.itemgetter(1))[0], win_size)
	calc_optimum_seq(aln_seq, max(combined_scores.iteritems(), key=operator.itemgetter(1))[0],win_size)
	with open(out_file, 'wb') as csvfile:
		mycsv = csv.writer(csvfile)


		for key,value in combined_scores.iteritems():  #prints headers with normalised aligned sRNAs
				mycsv.writerow([key,value])	
	out_pdf = out_file.split('.')[0]
	R_command = 'Rscript ~/Software/dsRNA_target_window/heatmap.R ' + \
		out_file + ' ' + out_pdf
	os.system(R_command)
	sys_command = 'open Rplots.pdf'
	os.system(sys_command)


def calc_optimum_seq(aln_seq, win_pos, win_size):
	"""returns optimum dsRNA construct at the supplied position by caclulating
	the most common nucelotide at each position in the window
	"""
	out_seq=""

	for pos in range(win_pos, win_pos+win_size):
		nuc_pos = []	
		for header, seq in aln_seq.iteritems():
			nuc_pos.append(seq[pos]) #make a list of nucs for that position
		if Counter(nuc_pos).most_common(1)[0][0]!='-':
			out_seq += Counter(nuc_pos).most_common(1)[0][0] #count their occurences
	print "Optimum dsRNA construct seq of length {0} at alignment position {1} is:\n{2}".format(len(out_seq), win_pos, out_seq)	

def comline():
	parser = argparse.ArgumentParser("dsRNA_target_window.py")
	parser.add_argument('input_alignment', type=str)
	parser.add_argument('output_file', type = str)
	parser.add_argument('-win_size','--window_size', default = '300', type=int)
	parser.add_argument('-pos_bonus','--pos_bonus_mod', default = '5', type=int)
	parser.add_argument('-run_bonus','--run_bonus_mod', default = '10', type=int)

	args = parser.parse_args()
	return args

def run_ana():
	a = comline()
	window_analysis(a.input_alignment, a.output_file, a.window_size, a.pos_bonus_mod, a.run_bonus_mod)

run_ana()










#window_analysis('./sge1_cds_alignment.fasta', 'sge1_cds_alignment.csv')