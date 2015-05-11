import os, subprocess, platform

for genome_filename in os.listdir('genomes'):
	if not genome_filename.endswith('.fna'):
		continue

	if platform.system() == "Linux":
	    args = [
		    './art_bin_ChocolateCherryCake_linux/art_bin_ChocolateCherryCake/art_illumina', 
		    '-i', 'genomes/'+genome_filename,
		    '-l', '150',
		    '-f', '100',
		    '-ss', 'MS',
		    '-o',  'reads/'+genome_filename.replace('.fna',''),
		    '-na'
	    ]
	else:
	    args = [
		    './art_bin_ChocolateCherryCake/art_illumina', 
		    '-i', 'genomes/'+genome_filename,
		    '-l', '150',
		    '-f', '100',
		    '-ss', 'MS',
		    '-o',  'reads/'+genome_filename.replace('.fna',''),
		    '-na'
	    ]

	# ./art_bin_ChocolateCherryCake/art_illumina -i genomes/AE013218.fna  -l 150 -f 100 -ss MS -o test/AE013
	subprocess.call(args)


