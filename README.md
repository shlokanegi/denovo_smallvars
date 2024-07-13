# Denovo pipeline for LRS Small Variants
This workflow generates an annotated de novo VCF using trio variant calls. It uses pairwise vcfeval, along with custom-filtering to generate annotated denovo candidates. The anontations include -
* Predicted impact on the genes using SNPeff
* gnomAD's allele frequencies, conservation/CADD scores, and ClinVar's clinical significance.
* Custom de novo read validation flag which helps filter confident de novos.

## Workflow tasks

1. `run_vcfeval` : Runs vcfeval on child and parents' VCFs (trio) to generate total denovo variants per family. Filters denovos by variant QUAL threshold.

2. `subset_denovos_by_region` : Restricts to denovos outside low-complexity regions, i.e. excludes variants in all tandem repeats, homopolymers and satellites. This step could be turned off by not providing this BED file as input.

3. `run_filtering` : Performs effective denovo variant filtering using parent's VCFs to remove possible false-positives. It only keeps variants which are not seen in both parents' VCFs.

4. `annotate_with_gnomad` : Annotates VCF with gnomad allele frequencies (if gnomad VCF is provided).

5. `annotate_with_snpeff` : Annotates VCF with SnpEff annotations (if snpEff database files are provided).

6. `subset_annotate_smallvars_with_db` : 
    - First, subsets to variants with high/moderate impact or with predicted loss of function. 
    - Then, annotates SNPs with presence in ClinVar and some dbNSFP annotations (CADD and GERP++) (if clinvar VCF and dbNSFP database files are provided).
    - This is returned a separate VCF by the workflow.

7. `keep_rare_denovos` : If `KEEP_RARE = true` (default) is set, then subsets to rare de novo snps using gnomad annotations. This keeps both rare variants from gnomad (AF < 0.001), as well as variants not seen in gnomad (not annotated with the AF tag).

8. `validate_denovos` : Validates de novos using parents' BAMs (if provided) by generating a gVCF with variants at candidate sites. Flag de novos as true positives when there's no alternate call in the parents.

## Test locally
```sh
## download GRCh38.105 database
wget https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_GRCh38.105.zip
## Run with miniwdl
miniwdl run --as-me -i test.inputs.json workflow.wdl
```

## Test with Toil
```sh
toil-wdl-runner workflow.wdl --inputs test.inputs.json
```

Workflow wraps rtgtools vcfeval and validation scripts to generate annotated denovo candidates in a trio.