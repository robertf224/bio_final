import os, cPickle, random, heapq, sys
from utils import *
import numpy as np
import scipy.optimize as op

"""
	usage: compute <k> <sample name>
"""
k = int(sys.argv[1])
sample_name = sys.argv[2]

@memoize
def reads_loglikelihoods(sample_name, k, cutoff=2):
	sample_filename = 'samples/%s.txt' % sample_name

	kmer_spectra_filename = 'pickles/kmer_spectra_%d.pickle' % k
	with open(kmer_spectra_filename) as f:
		kmer_spectra = cPickle.load(f)

	sample_kmers_filename = 'pickles/%s_kmers_%d.pickle' % (sys.argv[1], k)
	with open(sample_kmers_filename) as f:
		sample_kmers = cPickle.load(f)

	# get num reads for progress updating
	reads = open(sample_filename)
	num_reads = sum(1 for line in reads)
	reads.seek(0)

	# denominators of likelihoods for individual kmers
	kmer_totals = total_kmers_per_genome()

	# what we will return
	reads_loglikes = np.empty([num_reads, 20], dtype=float)

	for read_index, read in enumerate(progress(reads, length=num_reads)):
		# get the kmers from this read that are frequent kmers of this sample (and thus, probably not erroneous)
		frequent_kmers = filter(lambda kmer: kmer in sample_kmers and sample_kmers[kmer] >= cutoff, kmers(read, k))

		for kmer in frequent_kmers:
			for genome_index, count in enumerate(kmer_spectra[kmer] if kmer in kmer_spectra else [0]*20):
				if count == 0:
					like = 1.0 / max(kmer_totals)
				else:
					like = float(count) / kmer_totals[index]
				reads_loglikes[read_index][genome_index] += np.log(like)

	return reads_loglikes


def predict_alphas_simple(sample_name, k, cutoff=2):
	"""
		Look at the log likelihoods for each read, map a read to the genome that maximizes its likelihood

		Use proportions of read maps to estimate alphas
	"""
	read_map_counts = [0]*20
	reads_loglikes = read_loglikelihoods(sample_name, k, cutoff)

	for read_loglikes in reads_loglikes:
		read_map_counts[np.argmax(read_loglike)] += 1

	predicted_alphas = normalize_counts(read_map_counts)
	return predicted_alphas



