import sys, json
from utils import *
import numpy as np
from gen_alpha_plot import make_bar_plot

"""
	usage: compute <k> <sample name>
"""
k = int(sys.argv[1])
sample_name = sys.argv[2]

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

def predict_alphas_expectation(sample_name, k):
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

# output
predicted_alphas = predict_alphas_expectation(sample_name, k)
print 'predicted alpha distribution:' 
print predicted_alphas
print ''

with open('testcases/%s.json' % sample_name) as f:
	expected_alphas = json.load(f)['alphas']
print 'expected alpha distribution:'
print expected_alphas
print ''

plot_filename = 'plots/%s_%d.png' % (sample_name, k)
chart_title = '%s, k=%d' % (sample_name, k)
make_bar_plot(expected_alphas, predicted_alphas, chart_title, plot_filename)
print 'plot saved to %s' % plot_filename

