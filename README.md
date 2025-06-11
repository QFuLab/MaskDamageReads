# MaskDamageReads

MaskDamageReads.py is a command-line tool to mask the damaged bases in ancient DNA reads, and
it can also extract the reads with damaged bases. This script requires Python 3 and pysam (version >= 0.22.1).

Citation: Fu Q, Cao P, Dai Q, Bennett EA, Feng X, Yang MA, et al. Denisovan mitochondrial DNA from dental calculus of the >146,000-year-old Harbin cranium. Cell. 2025

Code released under GNU General Public License v3.0.

## Overview

This script performs masking of potentially damaged bases in a BAM file, either by a user-specified length from the ends or across the entire read (specified using the '--termini' parameter). It outputs a masked BAM file where damaged bases are replaced with 'N', and the quality scores for these bases are set to 0. Masking is done in a strand-aware manner, depending on whether the library is single-stranded or double-stranded (specified using the --lib parameter). The script can also output the coordinates of each 'N' added by the script into a BED file with the '--output-masked-position' flag.

In default mode, masking is only applied to the potentially damaged aligned read bases that differ from the reference base. This means that masking occurs only when the aligned read base is T or A and the reference base is C or G. If the read base is C or G, or if the read and reference bases are identical (e.g., T or A), no masking will occur.

The default mode may introduce reference bias when the sequencing depth is too low because the script masks all C>T and G>A variations. If a site is covered only once, true C>T and G>A mutations will also be masked. However, if a site is covered more than once, the real alternative read bases (those not at the ends) will be retained.

For low sequencing depth data, we recommend using the '--mode low_coverage' flag. In this mode, the script changes all C and T, G, and A read bases to 'N', regardless of the reference base. Although this mode removes more data, it provides more accurate allele frequency estimates.

No genome reference file is required, as the reference sequence can be reconstructed from the MD tag in the BAM file, so the BAM file **must** contain the MD tag.

With the '--extract-damaged' flag, this script will generate a BAM file containing only reads that have at least one damaged base (T or A, differ from reference base C or G ) within the region specified by the '--termini' parameter.

## Usage

```
MaskDamageReads.py [-h] -i INPUT -t TERMINI -l {SS,DS} -o OUTPUT [-m {normal,low_coverage}] [--extract-damaged] [--output-masked-position]

Masking the terminal damaged bases in BAM file. This script only works with bam file containing MD tag. The command `samtools calmd -b in.bam ref.fa > out.bam` can add MD tag for you.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input BAM file
  -t TERMINI, --termini TERMINI
                        Number bp of the terminal bases to mask. If it is set to 'all', the script will mask the damaged bases in the whole read ('-t all' only works with '-m SS'). If the 5'' and 3'' ends need different base number, you need to seperate them by comma (no space between comma and the second number), and put the 5'' number first. Should be attention that the 5'' and 3'' here refer to the reference direction (Value: an integer, 'all', or 'integer,integer')
  -l {SS,DS}, --lib {SS,DS}
                        The way of library preparation: SS (single-strand library) or DS (double-strand library)
  -o OUTPUT, --output OUTPUT
                        Output BAM file
  -m {normal,low_coverage}, --mode {normal,low_coverage}
                        Mode selection: normal or low_coverage (Default: normal; see README for details)
  --extract-damaged     Flag to extract damaged reads (Default: False)
  --output-masked-position
                        Flag to output the coordinates of each N that this script added into a bed file (Default: False)
```

## Example

* Masking 5 bps from the 5' and 3' ends in a single-strand library BAM file: 

  `python ./MaskDamageReads.py -i SS.bam -t 5 -l SS -o SS.masked.t5.bam` 

  ![Masked SS BAM](/fig/SS.mask.png "SS BAM before and after masking")

* Masking 1 bp from the 5' end and 2 bps from the 3' end in a double-strand library BAM file: 

  `python ./MaskDamageReads.py -i DS.bam -t 1,2 -l DS -o DS.masked.t1-2.bam` 
  
  ![Masked DS BAM](/fig/DS.mask.png "DS BAM before and after masking")

* Extracting reads with damaged 1st base from the 5' end in a double-strand library BAM file: 

  `python ./MaskDamageReads.py -i DS.bam -t 1,0 -l DS -o DS.extracted.t1-0.bam --extract-damaged`
  
  ![Extracted DS BAM](/fig/DS.extract.png "DS BAM before and after extracting")

## Contact

Please report bugs and suggestions on this Github page: [issues](https://github.com/white-sail-dev/MaskDamageReads/issues).
