Container to validate denovos by generating a pileup of reads at candidate de novo site in both parents.

The [`dnvval.py`](dnvval.py) script is included (available at `/opt/scripts/dnvval.py` in the container).
It will run samtools pileup for each variant with both mother's and father's BAM file, parse matches and mismatches and annotate VCF records with a true-positive flag if none of the the parents' reads support the variant.

The inputs of `dnvval.py` are:

```
root@3df923e5c1cd:/home# python3 /opt/scripts/dnvval.py --help
usage: dnvval.py [-h] -mbam MBAM -dbam DBAM -r R -v V [-o O] [-d D] [-t T]

optional arguments:
  -h, --help  show this help message and exit
  -mbam MBAM  BAM file (indexed)
  -dbam DBAM  BAM file (indexed)
  -r R        reference FASTA file (indexed)
  -v V        VCF file (can be bgzipped)
  -o O        output (annotated) VCF (will be bgzipped if ending in .gz)
  -d D        output directory
  -t T        number of threads used by bcftools mpileup
```
One new INFO field is added in the output VCF:

- `TP` with a value of '1', if the variant was not found in any reads in both parents.

To build locally and upload to [quay.io](https://quay.io/shnegi/dnvval).
