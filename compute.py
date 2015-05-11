import os, cPickle, random, heapq, sys
from utils import *
import numpy as np
import scipy.optimize as op

"""
	usage: compute <k> <sample filename>
"""
k = int(sys.argv[1])
sample_filename = sys.argv[2]

kmer_spectra_filename = 'pickles/kmer_spectra_%d.pickle' % k
with open(kmer_spectra_filename) as f:
	kmer_spectra = cPickle.load(f)

def lnlikes_to_proportions(lnlikes):
	basepoint = max(lnlikes)
	proportions = [0.0]*20
	for index, lnlike in enumerate(lnlikes):
		proportions[index] = np.exp(basepoint - lnlike) * 1.0
	return proportions

def normalize_counts(counts):
	total = sum(counts)
	return map(lambda x: x / float(total), counts)

def predict_alphas(max_num_species=20):
	# get total kmer counts for each genome
	kmer_totals = total_kmers_per_genome()

	# get num reads for progress updating
	reads = open(sample_filename)
	num_reads = sum(1 for line in reads)
	reads.seek(0)

	read_map_counts = [0]*20
	for read in progress(reads, length=num_reads):
		# compute log-likelihoods for each genome for this read
		lnlikes = [0.0]*20
		for kmer in kmers(read, k):
			for index, count in enumerate(kmer_spectra[kmer] if kmer in kmer_spectra else [0]*20):
				if count == 0:
					# smoothing
					like = 1.0 / max(kmer_totals)
				else:
					like = float(count) / kmer_totals[index]
				lnlikes[index] += np.log(like)

		# map to best one
		max_index = 0
		for index, lnlike in enumerate(lnlikes):
			if lnlikes[index] > lnlikes[max_index]:
				max_index = index

		read_map_counts[max_index] += 1

	predicted_alphas = normalize_counts(read_map_counts)
	return predicted_alphas



