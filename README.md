# Megamap
Automates the mapping of multiple sequence reads files against many different reference genomes.

Megamap.sh pulls the names of references and read sets from maplist.txt and sends them, one pair at a time to quickmap.sh which does all the rest of the work calling BWA and Samtools.
The results are packaged in a folder for each read/reference pairing. Additional files include those for use with IGV and statistics on lengths and coverage as well as counts of indels and mutations.

You must modify the maplist file by replacing your read and reference files into the list as depicted in the example maplist.txt file.

Place Megamap, quickmap and maplist into the bwa folder for mapping along with your fasta read and reference files.

See Megaparse for parsing mpile output files.
See sam2circos for parsing sam files for use with Circos mapping of chimreic read locations.
