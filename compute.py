import os, cPickle, random, heapq, sys
from utils import kmer_store, kmers, nucleotides_fna, progress
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

# get totals for each kmer
kmer_totals = [0]*20
for kmer in kmer_spectra:
	for index, count in enumerate(kmer_spectra[kmer]):
		kmer_totals[index] += count

reads = open(sample_filename)
read_map_counts = [0]*20
for read in reads:
	# choose best match in genomes
	for kmer in kmers(read):


estimated_proportions