# Denovo pipeline for LRS Small Variants
This workflow generates an annotated de novo VCF using trio variant calls. It uses pairwise vcfeval, along with custom-filtering to generate annotated denovo candidates. The anontations include -
* Predicted impact on the genes using SNPeff
* gnomAD's allele frequencies, conservation/CADD scores, and ClinVar's clinical significance.
* Custom de novo read validation flag which helps filter confident de novos.

## Workflow tasks

1. `run_vcfeval` : Runs vcfeval on child and parents' VCFs (trio) to generate total denovo variants per family. Filters denovos by variant QUAL threshold.

2. `run_filtering` : Performs effective sex-aware denovo variant filtering using parent's VCFs to remove possible false-positives. Removes variants called as 0/0 or ./. in both parents.

3. `annotate_with_gnomad` : Annotates VCF with gnomad allele frequencies (if gnomad VCF is provided).

4. `annotate_with_snpeff` : Annotates VCF with SnpEff annotations (if snpEff database files are provided).

5. `subset_annotate_smallvars_with_db` : 
    - First, subsets to variants with high/moderate impact or with predicted loss of function. 
    - Then, annotates SNPs with presence in ClinVar and some dbNSFP annotations (CADD and GERP++) (if clinvar VCF and dbNSFP database files are provided).
    - This is returned a separate VCF by the workflow.

6. `validate_denovos` : Validates de novos using parents' BAMs (if provided) by generating a gVCF with variants at candidate sites. Flags de novos with based on read-support for alternate call in parents.

7. `keep_rare_denovos` : If `KEEP_RARE = true` (default) is set, then subsets to rare de novos using gnomAD AF. This keeps both rare variants from gnomad (AF < 0.001), as well as variants not seen in gnomad (not annotated with the AF tag).


## Test locally
```sh
## download GRCh38.105 database
wget https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_GRCh38.105.zip
## Run with miniwdl
miniwdl run --as-me -i test.inputs.json wdl/workflow.wdl
```

## Test with Toil
```sh
toil-wdl-runner wdl/workflow.wdl --inputs test.inputs.json
```