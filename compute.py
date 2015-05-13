import sys, json
from utils import *
import numpy as np
import scipy.optimize as op
from gen_alpha_plot import make_bar_plot

"""
	usage: compute <k> <sample name>
"""

def predict_alphas_likelihood(sample_name, k, full=False, cutoff=2, smoothing_function=None):
	"""
		Look at the log likelihoods for each read, map a read to the genome that maximizes its likelihood

		Use proportions of read maps to estimate alphas
	"""
	read_map_counts = [0]*20
	reads_loglikes = reads_loglikelihoods(sample_name, k, full, cutoff, smoothing_function)

	for read_loglikes in reads_loglikes:
		read_map_counts[np.argmax(read_loglikes)] += 1

	predicted_alphas = normalize_counts(read_map_counts)
	return predicted_alphas

def predict_alphas_sum_of_squares(sample_name, k):
	"""
		Compute alphas by minimizing sum of squared differences for each kmer of expected versus observed
	"""
	N = 100000
	L = 150
	kmer_spectra = load_kmer_spectra(k)
	sample_kmers = load_sample_kmers(sample_name, k)
	sizes = total_kmers_per_genome()

	def sum_squared_differences(alphas):
		"""
			Compute sum of squared differences between predicted # kmers and found for recognized kmers
		"""
		alphas = normalize_counts(alphas)
		total_sum_squared_diff = 0.0
		for kmer in sample_kmers:
			if kmer in kmer_spectra:
				counts = kmer_spectra[kmer]
				total_sum_squared_diff += (sum(N*(L-k+1)*alphas[i]*(float(counts[i])/sizes[i]) for i in xrange(20)) - sample_kmers[kmer])**2
		return total_sum_squared_diff

	bounds = [(0, None) for i in xrange(20)]
	objective = lambda alphas: sum_squared_differences(alphas)
	priors = [0.05]*20
	prediction = op.minimize(objective, priors, bounds=bounds)
	return normalize_counts(prediction['x'])

def predict_alphas_uniques_expectation(sample_name, k):
	"""
		Predict alphas using # unique k-mers observed as an estimator of E(# unique k-mers from G_i)
	"""
	sample_filename = 'samples/%s.txt' % sample_name
	reads = open(sample_filename)
	num_reads = sum(1 for line in reads)
	reads.seek(0)

	kmer_spectra = load_kmer_spectra(k)
	genome_sizes = total_kmers_per_genome()

	# calling uniques k-mers that appear only in one genome, and appear at least 4 times in that genome
	uniques = set(filter(lambda kmer: len(filter(bool, kmer_spectra[kmer])) == 1 and sum(kmer_spectra[kmer]) >= 4, kmer_spectra))
	unique_counts = [0]*20
	for kmer in uniques:
		index = np.argmax(kmer_spectra[kmer])
		unique_counts[index] += kmer_spectra[kmer][index]

	unique_frequencies = [0]*20
	for read in progress(reads, length=num_reads):
		for kmer in kmers(read.strip(), k):
			if kmer in uniques:
				index = np.argmax(kmer_spectra[kmer])
				unique_frequencies[index] += 1
	reads.close()

	for index, count in enumerate(unique_frequencies):
		unique_frequencies[index] /= (unique_counts[index] / float(genome_sizes[index]))
		# unique_frequencies[index] /= 100000*(150-k+1) # theoretically what we're doing

	return normalize_counts(unique_frequencies)

if __name__ == '__main__':
	k = int(sys.argv[1])
	sample_name = sys.argv[2]
	version = sys.argv[3]
	method = {
		'l': predict_alphas_likelihood,
		's': predict_alphas_sum_of_squares,
		'u': predict_alphas_uniques_expectation
	}[version]

	# output
	predicted_alphas = method(sample_name, k)
	print 'predicted alpha distribution:' 
	print predicted_alphas
	print ''

	with open('testcases/%s.json' % sample_name) as f:
		expected_alphas = json.load(f)['alphas']
	print 'expected alpha distribution:'
	print expected_alphas
	print ''

	plot_filename = 'plots/%s_%d_%s.png' % (sample_name, k, version)
	chart_title = '%s, k=%d, method=%s' % (sample_name, k, version)
	make_bar_plot(expected_alphas, predicted_alphas, chart_title, plot_filename)
	print 'plot saved to %s' % plot_filename

