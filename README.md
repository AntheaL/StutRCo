**How to run**

python correct_errors.py --bam [BAM FILE] --bed [BED FILE] --dst-path [CORRECTED BAM PATH]

**What it does**

A 'UMI family' is a group of reads having the same UMI and covering a specific region in the bam file.

For each site (*s*) in the bed file, for each cell (*c*), we find the UMI families (*umi_families*) corresponding to *(s,c)*.
<br> For each *family* in *umi_families*:
- we select the most frequent indel value across the reads (in case of equality, we select the most frequent indel value across *umi_families*)
- we keep a matching read or build a consensus read out of the reads in *family*

Finally, we save the bam in a new file with one read per corrected UMI family.
