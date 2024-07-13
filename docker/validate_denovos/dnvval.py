import argparse
from subprocess import run
from cyvcf2 import VCF, Writer
import os
import sys


def alt_present(pileup_output_path, alt):
    with open(pileup_output_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                alt_alleles = line.rstrip().split()[4].split(",")
                # If the alternate allele is present
                if alt in alt_alleles:
                    return True
                else:
                    return False


def evaluate_dnv(dnv_info, ref_fa_path, mom_bam_path, dad_bam_path, output_dir, threads=2):
    # dump = open('/dev/null', 'w')
    region = '{}:{}-{}'.format(dnv_info['chrom'], dnv_info['pos'], dnv_info['pos'])
    depth = '1000'

    # bcftools pileup at variant position in mother's BAM
    pileup_args = ["bcftools", "mpileup", "-r", region, "--threads", str(threads), "-d", depth, "-f", ref_fa_path, mom_bam_path]
    mom_pileup_ouput_path = os.path.join(output_dir, "mom.vcf")
    with open(mom_pileup_ouput_path, 'w') as file:
        run(pileup_args, check=True, stdout=file,
            stderr=sys.stderr, universal_newlines=True)
        
    if alt_present(mom_pileup_ouput_path, dnv_info['alt'][0]):
        return "FP"
    
    # bcftools pileup at variant position in father's BAM
    else:       
        pileup_args = ["bcftools", "mpileup", "-r", region, "--threads", str(threads), "-d", depth, "-f", ref_fa_path, dad_bam_path]
        dad_pileup_ouput_path = os.path.join(output_dir, "dad.vcf")
        with open(dad_pileup_ouput_path, 'w') as file:
            run(pileup_args, check=True, stdout=file,
                stderr=sys.stderr, universal_newlines=True)
        if not alt_present(dad_pileup_ouput_path, dnv_info['alt'][0]):
            return "TP"
        else:
            return "FP"


parser = argparse.ArgumentParser()
parser.add_argument('-mbam', help='BAM file (indexed)', required=True)
parser.add_argument('-dbam', help='BAM file (indexed)', required=True)
parser.add_argument('-r', help='reference FASTA file (indexed)', required=True)
parser.add_argument('-v', help='VCF file (can be bgzipped)',required=True)
parser.add_argument('-o', default='out.vcf', help='output (annotated) VCF (will be bgzipped if ending in .gz)')
parser.add_argument('-d', help='output directory', default='temp_valsv')
parser.add_argument('-t', default=2, help='number of threads used by bcftools mpileup')
args = parser.parse_args()

vcf = VCF(args.v)
vcf.add_info_to_header({'ID': 'TP', 
                        'Description': 'Flagging True-positive denovos identifed with in-silico denovo validation',
                        'Type': 'Float', 'Number': '1'})

vcf_o = Writer(args.o, vcf)

# Read VCF and evaluate each denovo SNP
for variant in vcf:
    dnvinfo = {}
    dnvinfo['chrom'] = variant.CHROM
    dnvinfo['pos'] = variant.POS
    dnvinfo['alt'] = variant.ALT

    flag = evaluate_dnv(dnvinfo, args.r, args.mbam, args.dbam, args.d, threads=args.t)
    if flag == 'TP':
        variant.INFO['TP'] = '1'

    vcf_o.write_record(variant)

vcf_o.close()
vcf.close()