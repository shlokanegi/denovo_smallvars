# vcfeval workflow
Run RTG [vcfeval](https://github.com/RealTimeGenomics/rtg-tools) to compare compare any 2 VCFs and output intersection, exclusive VCFs and summary stats.

## Input considerations
* sample1.vcf Gzipped
* sample2.vcf Gzipped
* any unique name
* Reference genome fasta

## Docker container
The `run_vcfeval` task of the workflow uses the [quay.io/shnegi/vcfeval](https://quay.io/repository/shnegi/vcfeval). It contains rtgtools, bcftools and dependencies. It's made with the [Dockerfile](https://github.com/shlokanegi/wdl_workflows/blob/master/vcfeval_wf/Dockerfile) in this repo.

## Test locally
```sh
miniwdl run --as-me -i test.inputs.json workflow.wdl
```

## Test with Toil
```sh
toil-wdl-runner workflow.wdl --inputs test.inputs.json
```
