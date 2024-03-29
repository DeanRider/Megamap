# Megamap
Automates the mapping of multiple sequence reads files against many different reference genomes.

Megamap2.sh pulls the names of references and read sets from maplist.txt and sends them, one pair at a time to quickmap2.sh which does all the rest of the work calling BWA and Samtools. The December 2022 update removed the 8000 read limit on mpileup and added annotations to the log file.
The results are packaged in a folder for each read/reference pairing. Additional files include those for use with IGV and statistics on lengths and coverage as well as counts of indels and mutations.

You must modify the maplist file by replacing your read and reference files into the list as depicted in the example maplist.txt file.

Place Megamap, quickmap and maplist into the bwa folder for mapping along with your fasta read and reference files.  Use bash to run the script.

- Following Megamap2.sh output includes:
   - Index files used for mapping which remain in the bwa folder
   - Data output for each read/reference set in a package containing:
     - log
     - sam
     - bam
     - sorted bam
     - sorted bai
     - mpileup
     - DepthOfCoverage
     - ContigLengthNumreads
     
- Requirements:
   - Mac OS (not likely to be POSIX compliant; originally implemented on Mojave but still working in Monterey)
   - Samtools v1.9 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
   - BWA-MEM v0.7.17-r1188 https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2

See Megaparse for parsing mpileup output files.
See sam2circos for parsing sam files for use with Circos mapping of chimreic read locations.
