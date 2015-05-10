import os, subprocess

for genome_filename in os.listdir('genomes'):
	if not genome_filename.endswith('.fna'):
		continue

	args = [
		'./art_bin_ChocolateCherryCake/art_illumina', 
		'-i', 'genomes/'+genome_filename,
		'-l', '150',
		'-f', '100',
		'-ss', 'MS',
		'-o',  'reads/'+genome_filename.replace('.fna','')
	]

	# ./art_bin_ChocolateCherryCake/art_illumina -i genomes/AE013218.fna  -l 150 -f 100 -ss MS -o test/AE013
	subprocess.call(args)


