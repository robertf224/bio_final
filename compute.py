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
		counts = map(lambda x: 1 if x == 0 else x, kmer_spectra[kmer] if kmer in kmer_spectra else [0]*20)
		for index, count in enumerate(counts):
			lnlikes[index] += np.log(float(count) / kmer_totals[index])
	
	# map to best one
	max_index = 0
	for index, lnlike in enumerate(lnlikes):
		if lnlikes[index] > lnlikes[max_index]:
			max_index = index

	read_map_counts[max_index] += 1

total_reads = sum(read_map_counts)
predicted_alphas = map(lambda x: x / float(total_reads), read_map_counts)
print predicted_alphas

