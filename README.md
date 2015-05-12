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

> Propose and implement a model for estimating the relative proportions α1, . . . , αk. Evaluate
> how well you can recover the true proportion in your simulated data set.

We're using the observed number of unique k-mers from a specific genome in a sample as an estimator for the expected number of unique k-mers from that genome in the sample given its alpha value.  Specifically, the expected number of uniques from a specific genome is:

![](http://i.imgur.com/IBlQu67.png?1)

The number of unique kmers for a genome should approach this expectation as the sample size increases, so we can use it as an estimator, rearrange the terms, and compute alpha values for every genome, normalizing at the end.

Deciding what are considered unique k-mers is a bit of a heuristic (we chose to consider those k-mers that appear in only one genome, and at least 4 times in that genome.  we also found that k=15 works well), but the model itself has an actual probabilistic motivation, and works fairly well.  Here is an example, generated from S3 in the testcases directory using k=15:

![S3, k=15](http://i.imgur.com/tGhoOSw.png?1)

We also have another method implemented that computes the likelihood of a read coming from a specific genome based on making draws from the k-mer spectrum, and maps that read to the genome that maximizes the likelihood.  The proportions are then estimated by the proportions of reads mapping to different genomes.  This works to some extent, but unfortunately it requires smoothing and is thus difficult to tune.

- Follow all of the above instructions: `sample generation`, `kmer generation`, and then lastly `Computation`. Running `compute.py <k> <sample name>` will yield alpha values and a plot

> Propose a method for deciding when novel species are likely to be present in a sample.

We can easily work novel species detection into our existing framework.  We can essentially add another alpha value to represent the proportion of novel species, and use the number of unrecognized unique kmers in the sample as our estimator for the expected number of unique kmers from novel species in the sample.

> Describe which information is required or which assumptions are necessary to model the problem
> of estimating the number of novel genomes present.

In general, we would need some information regarding the disjoint k-mer content of the novel genomes relative to the original taxa and to each other (i.e. how many unique k-mers we should expect from each, how many common k-mers we should expect from each)
