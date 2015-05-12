# Bioinformatics Final

## How to Run

### Sample generation

To generate all reads (do this once):
```bash
python read_simulation.py
```

To generate a sample from a sample definition json file (examples in testcases folder):
```bash
python sample_simulation.py <sample name>
The samples are stored in the samples folder.
```

### kmer generation

To generate kmers for all genomes (only need to do once for a specific k value):
```bash
python gen_kmers.py <k> genomes
```

To generate kmers for a sample:
```bash
python gen_kmers.py <k> sample <sample name>
```

kmer spectra are stored as pickle files in the pickles folder.

### Computation
```bash
python compute.py <k> <sample name>
```
***

## Writeup

> Compute the k-mer spectra for the different genomes from Sakai. Note that k-mer spectra will
> generally not be disjoint between different species.

- Follow the `kmer generation` instructions above, specifically the second one: `python gen_kmers.py <k> sample <sample name>`

> Create simulated test data by using the genomes from Sakai. Apply a simulator for Illumina
> reads (use 150bp reads and MiSeq simulation, e.g., http://www.niehs.nih.gov/
> research/resources/software/biostatistics/art/). Note, you only need to
> run the sequencing simulation once; use 100x coverage.
> Simulate a (multi-cell) meta-genome experiment by subsampling the reads from the simulation
> output according to the relative proportion of genomes, α1, . . . , αk. Repeat 10 times for
> different choices of α1, . . . , αk (e.g., 1 species present, 3-10 species present with uniform and
> non-uniform αi for the species present). Keep track of what goes into which data set, so that
> you can assess performance. Let S1, . . . , S10 be the test datasets.
> Question: How much sequencing should be done in total for the meta-genome?

- To create the reads and samples, follow the `Sample generation` instructions above.

- We ended up trying 10,000, 100,000 and 1,000,000 reads and found, experimentally, that 100,000 worked the best. 

> From a sequencing data set S, generate all k-mers and identify all k-mers appearing twice or
> more in the data, for example using a Bloom-filter and a dictionary. This yields k-mer frequency
> spectra F(S) = Fk(S).

- Follow the `kmer generation` instructions above, specifically the second one: `python gen_kmers.py <k> genomes`

> Devise and implement a probabilistic model which assigns a likelihood to a species being
> present given the k-mer frequency spectrum Fk(S). Evaluate how well you can identify presence
> of organisms, both when exactly one or more than one species are present.

> Devise a strategy which improves running time of the analysis and/or performance based on the
> observation that there are k-mers which appear in multiple genomes.

> Propose and implement a model for estimating the relative proportions α1, . . . , αk. Evaluate
> how well you can recover the true proportion in your simulated data set.

- Follow all of the above instructions: `sample generation`, `kmer generation`, and then lastly `Computation`. Running `compute.py` will yield alpha values

> Propose a method for deciding when novel species are likely to be present in a sample.

> Describe which information is required or which assumptions are necessary to model the problem
> of estimating the number of novel genomes present.
