import json, sys, os
from utils import reservoir_sample, reads_fq, progress

"""
	args = {
		N: <number of reads total>
		alphas: [a1, a2, ..., a20]
	}
"""
if len(sys.argv) != 2:
	print 'usage: sample_simulation <sample name>'
sample_name = sys.argv[1]
args_filename = 'testcases/%s.json' % sample_name
output_filename = 'samples/%s.txt' % sample_name

with open(args_filename) as f:
	args = json.load(f)

with open(output_filename, 'w') as f:
	i = 1
	for reads_filename in progress(os.listdir('reads')):
		if not reads_filename.endswith('.fq'):
			continue
		reads_filename = 'reads/'+reads_filename

		alpha = args['alphas'][i-1]
		i += 1
		subsample_size = int(args['N'] * alpha)
		if subsample_size == 0:
			continue

		reads = reads_fq(reads_filename)
		subsample = reservoir_sample(reads, subsample_size)
		for read in subsample:
			f.write(read+'\n')





