# k
# sample name
# smoothing function
import subprocess 
import os

def main(k, sample_name):
   # assume reads are present
    if not os.path.exists("samples/%s"%sample_name):
	subprocess.call(["python", "sample_simulation.py", sample_name])

    if not os.path.exists("pickles/kmer_spectra_%d.pickle"%k):
	subprocess.call(["python", "gen_kmers.py", str(k), "genomes"])

    if not os.path.exists("pickles/%s_kmers_%d.pickle"%(sample_name, k)):
	subprocess.call(["python", "gen_kmers.py", str(k), "sample", sample_name])

if __name__ == "__main__":
    if len(sys.argv) != 3:
	print "python run.py <k> <sample name>"

    main(int(sys.argv[1]), sys.argv[2])
