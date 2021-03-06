import sys, cPickle, os
from collections import defaultdict
from utils import kmer_store, kmers, nucleotides_fna, progress

""""
	gen_kmers <k> sample <sample name>

		or

	gen_kmers <k> genomes <opt:full>
"""
k = int(sys.argv[1])
version = sys.argv[2]

if version == 'sample':
	sample_name = sys.argv[3]
	filename = 'samples/%s.txt' % sample_name
	sample_kmers = kmer_store()
	with open(filename) as f:
		for read in f:
			read_kmers = kmers(read.strip(), k)
			for kmer in read_kmers:
				sample_kmers.update(kmer)

	output_filename = 'pickles/%s_kmers_%d.pickle' % (os.path.basename(os.path.normpath(filename)).replace('.txt',''), k)
	with open(output_filename, 'w') as f:
		cPickle.dump(sample_kmers.kmers, f)

elif version =='genomes':
	full = True if (len(sys.argv) == 4 and sys.argv[3] == 'full') else False
	kmer_spectra = defaultdict(lambda:[0]*20)
	for index, genome_filename in enumerate(progress(filter(lambda x: x.endswith('.fna'), os.listdir('genomes')))):
		kmer_spectrum = {} if full else kmer_store()
		for kmer in kmers(nucleotides_fna('genomes/'+genome_filename), k):
			if full:
				kmer_spectrum[kmer] = kmer_spectrum[kmer]+1 if kmer in kmer_spectrum else 1
			else:
				kmer_spectrum.update(kmer)
		for kmer in kmer_spectrum:
			kmer_spectra[kmer][index] = kmer_spectrum[kmer]

	full_string = 'full_' if full else ''
	with open('pickles/%skmer_spectra_%d.pickle' % (full_string, k), 'w') as f:
		cPickle.dump(dict(kmer_spectra), f)

