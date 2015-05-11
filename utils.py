from pybloom import ScalableBloomFilter
from collections import deque
import random, sys, os

class kmer_store:
	def __init__(self):
		self.bloom_filter = ScalableBloomFilter(initial_capacity=1000000, mode=ScalableBloomFilter.LARGE_SET_GROWTH)
		self.kmers = {}

	def update(self, item):
		if item in self.bloom_filter:
			if item in self.kmers:
				self.kmers[item] += 1
			else:
				self.kmers[item] = 2
		else:
			self.bloom_filter.add(item)

	def __iter__(self):
		for key in self.kmers:
			yield key
	def __getitem__(self, key):
		return self.kmers[key]
	def __repr__(self):
		return str(self.kmers)
	def __str__(self):
		return str(self.kmers)
		
def nucleotides_fna(filename):
	with open(filename) as f:
		f.next()
		for line in f:
			line = line.strip()
			for char in line: 
				yield char

def kmers(nucleotides, k):
	current_kmer = deque()
	for nucleotide in nucleotides:
		current_kmer.append(nucleotide)
		if len(current_kmer) == k:
			yield ''.join(current_kmer)
			current_kmer.popleft()

def reads_fq(filename):
	prev_line = ''
	with open(filename) as f:
		for line in f:
			if line.strip() == '+':
				yield prev_line
			prev_line = line.strip()

def total_kmers_per_genome():
	totals = [0]*20
	if os.path.exists('pickles/totals.txt'):
		with open('pickles/totals.txt') as f:
			for index, line in enumerate(f):
				if index == 20:
					break
				totals[index] = int(line.strip())
		return totals
	for index, genome_filename in enumerate(progress(filter(lambda x: x.endswith('.fna'), os.listdir('genomes')))):
		totals[index] = sum(1 for kmer in kmers(nucleotides_fna('genomes/'+genome_filename), 5))
	with open('pickles/totals.txt', 'w') as f:
		for total in totals:
			f.write('%d\n'%total)
	return totals

def normalize_counts(counts):
	total = sum(counts)
	return map(lambda x: x / float(total), counts)

def memoize(f):
    memo = {}
    def helper(x):
        if x not in memo:            
            memo[x] = f(x)
        return memo[x]
    return helper
	
def reservoir_sample(iterator, size):
    sample = []
    observed = 0

    for item in iterator:
        observed += 1
        if len(sample) < size:
            sample.append(item)
        else:
            index = random.randint(0, observed-1)
            if index < size:
                sample[index] = item

    return sample

def progress(iterable, title=None, length=None):

	# print code (only print when we have to)
	last_num_bars = last_integer_percentage = None
	def print_bar(percentage):
		num_bars, integer_percentage = int(percentage*20), int(percentage*100)
		if num_bars == last_num_bars and integer_percentage == last_integer_percentage:
			return
		sys.stdout.write('\r')
		sys.stdout.write("[%-20s] %d%%" % ('=' * num_bars, integer_percentage))
		sys.stdout.flush()

	# print title regardless
	if title:
		print title

	max_val = len(iterable) if not length and hasattr(iterable, '__len__') else length
	if max_val:
		for index, element in enumerate(iterable):
			# print progress
			percentage = float(index+1) / max_val
			print_bar(percentage)
			yield element
			if index+1 >= max_val:
				break

		print_bar(1)

	else:
		sys.stdout.write('Can\'t calculate progress')
		sys.stdout.flush()
		for element in iterable:
			yield element

	sys.stdout.write('\n')
	sys.stdout.flush()

