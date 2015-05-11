# Bioinformatics Final

## Sample generation

To generate all reads (do this once):
```bash
python read_simulation.py
```

To generate a sample from a sample definition json file (examples in testcases folder):
```bash
python sample_simulation.py <sample name>
The samples are stored in the samples folder.
```

## kmer generation

To generate kmers for all genomes (only need to do once for a specific k value):
```bash
python gen_kmers.py <k> genomes
```

To generate kmers for a sample:
```bash
python gen_kmers.py <k> sample <sample name>
```

kmer spectra are stored as pickle files in the pickles folder.

## Computation

Preliminary models in compute.py