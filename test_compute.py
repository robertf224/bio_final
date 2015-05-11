import os, cPickle, random, heapq
from utils import kmer_store, kmers, nucleotides_fna, progress
import numpy as np
import scipy.optimize as op

k = 5
kmer_spectra_filename = 'pickles/kmer_spectra_%d.pickle' % k
with open(kmer_spectra_filename) as f:
	kmer_spectra = cPickle.load(f)

pref = 'S3'
sample_kmers_filename = 'pickles/%s_kmers_%d.pickle' % (pref, k)
with open(sample_kmers_filename) as f:
	sample_kmers = cPickle.load(f)

###############
# simple method (works pwell with k=15)
###############
def predict_simple(cutoff=1, max_num_species=20, percentile=100.0):
	# heuristically filter sample_kmers beforehand
	top_percentile_cutoff = heapq.nlargest(int(len(sample_kmers)/(100/percentile)), sample_kmers.values())[-1] if percentile != 100.0 else 0
	top_sample_kmers = {}
	for kmer in sample_kmers:
		if sample_kmers[kmer] >= top_percentile_cutoff:
			top_sample_kmers[kmer] = sample_kmers[kmer]

	expected_counts = [0.0]*20
	for kmer in top_sample_kmers:
		num = top_sample_kmers[kmer]
		counts_across_spectra = kmer_spectra[kmer] if kmer in kmer_spectra else [0]*20
		num_in_all = sum(counts_across_spectra)
		if num_in_all >= cutoff and len(filter(bool, counts_across_spectra)) <= max_num_species:
			for index, count in enumerate(counts_across_spectra):
				expected_counts[index] += num * (count / float(num_in_all))

	total = sum(expected_counts)
	predicted_alphas = map(lambda x: x / total, expected_counts)
	return predicted_alphas
######################

# some true examples
true = [0.0]*18 + [0.5]*2
true3 = [0.0,0.0,0.2,0.0,0.0,0.0,0.0,0.8]+[0.0]*12

# generate random priors
def gen_priors():
	priors = [0.0] * 20
	points = sorted([random.random() for i in xrange(19)])
	last = 0.0
	for index, point in enumerate(points):
		priors[index] = point - last
		last = point
	priors[-1] = 1 - last
	return priors

kmer_totals = [0]*20
for kmer in kmer_spectra:
	for index, count in enumerate(kmer_spectra[kmer]):
		kmer_totals[index] += count

def lnlike(alphas, cutoff=1, smooth_value=1.0, max_num_species=20):
	lnlike = 0.0
	for kmer in sample_kmers:
		count_in_sample = sample_kmers[kmer]
		counts_across_spectra = kmer_spectra[kmer] if kmer in kmer_spectra else [0]*20
		if sum(counts_across_spectra) < cutoff or len(filter(bool, counts_across_spectra)) > max_num_species:
			continue
		s = 0.0
		for index, count in enumerate(counts_across_spectra):
			freq = float(counts_across_spectra[index])
			if freq == 0.0: freq = smooth_value
			s += alphas[index] * (freq / kmer_totals[index])
		lnlike += count * np.log(s)
	return lnlike

def predict_alphas(priors, maxiter=8, cutoff=1, smooth_value=1.0):
	print 'optimizing...'
	bounds = [(0.0, 1.0) for i in xrange(20)]
	constraints = [{'type': 'eq', 'fun': lambda alphas: sum(alphas) - 1.0}]
	nll = lambda *args: -lnlike(*args, cutoff=cutoff, smooth_value=smooth_value)
	prediction = op.minimize(nll, priors, bounds=bounds, constraints=constraints, options={'maxiter':maxiter})
	return prediction






