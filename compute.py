import os, cPickle, random, heapq, sys
from utils import *
import numpy as np
import scipy.optimize as op

"""
	usage: compute <k> <sample name>
"""
k = int(sys.argv[1])
sample_name = sys.argv[2]


def predict_alphas_simple(sample_name, k, cutoff=2):
	"""
		Look at the log likelihoods for each read, map a read to the genome that maximizes its likelihood

		Use proportions of read maps to estimate alphas
	"""
	read_map_counts = [0]*20
	reads_loglikes = reads_loglikelihoods(sample_name, k, cutoff)

	for read_loglikes in reads_loglikes:
		read_map_counts[np.argmax(read_loglike)] += 1

	predicted_alphas = normalize_counts(read_map_counts)
	return predicted_alphas



