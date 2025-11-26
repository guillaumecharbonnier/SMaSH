#!/usr/bin/env python

from __future__ import print_function

#
# SomaCheck - Somatic Sample Checker
#
# Derived from SMaSH (Sample Matching using SNPs in Humans)
#
# This tool evaluates whether a sequencing sample comes from germline
# (normal) or from a somatic/cancer sample, given a VCF with known
# true somatic variants.
#
# The approach:
# - Read a VCF file containing known somatic variants (e.g., from tumor)
# - For each somatic variant position, count reads in the sample BAM
# - Calculate the variant allele frequency (VAF) at each position
# - Compute probabilities that the sample is somatic vs germline
# - A germline sample should NOT have these somatic mutations
# - A somatic sample SHOULD have these somatic mutations
#

import math
from scipy import stats
import argparse
from time import strftime
import pysam
import os
import sys


def eprint(*args, **kwargs):
    """Function that easily prints to stderr"""
    print(*args, file=sys.stderr, **kwargs)


def calculate_vaf(ref_count, alt_count):
    """Calculate variant allele frequency"""
    total = ref_count + alt_count
    if total == 0:
        return 0.0
    return alt_count / float(total)


def binomial_test_somatic(ref_count, alt_count, min_reads=5):
    """
    Perform a statistical test to determine if the observed VAF
    is significantly different from what we'd expect in a germline sample.
    
    For germline samples at somatic mutation sites:
    - We expect VAF = 0 (no variant reads)
    
    For somatic samples:
    - We expect VAF > 0 (some variant reads present)
    
    Returns a p-value for the hypothesis that this is a germline sample.
    Low p-value suggests this is a somatic sample.
    """
    total = ref_count + alt_count
    
    if total < min_reads:
        return 1.0  # Not enough evidence
    
    # If we see variant reads at a somatic site, that's evidence of somatic origin
    # Binomial test: probability of seeing alt_count or more variant reads
    # under the null hypothesis that VAF should be 0 (germline)
    # We use a small expected rate for sequencing errors (~0.01)
    error_rate = 0.01
    
    if alt_count == 0:
        return 1.0  # No variant reads, consistent with germline
    
    # One-tailed binomial test
    p_value = stats.binom.sf(alt_count - 1, total, error_rate)
    return p_value


def log_likelihood_ratio(ref_count, alt_count, expected_somatic_vaf=0.3):
    """
    Calculate log-likelihood ratio for somatic vs germline hypothesis.
    
    Germline hypothesis: VAF ~ 0 (only sequencing errors)
    Somatic hypothesis: VAF ~ expected_somatic_vaf (typically 0.2-0.5 for heterozygous somatic)
    
    Returns positive value if evidence favors somatic, negative if favors germline.
    """
    total = ref_count + alt_count
    
    if total == 0:
        return 0.0  # No evidence
    
    error_rate = 0.01
    
    # Avoid log(0) issues
    epsilon = 1e-10
    
    # Log-likelihood for germline (VAF ~ error_rate)
    ll_germline = (alt_count * math.log(max(error_rate, epsilon)) + 
                   ref_count * math.log(max(1 - error_rate, epsilon)))
    
    # Log-likelihood for somatic (VAF ~ expected_somatic_vaf)
    ll_somatic = (alt_count * math.log(max(expected_somatic_vaf, epsilon)) + 
                  ref_count * math.log(max(1 - expected_somatic_vaf, epsilon)))
    
    return ll_somatic - ll_germline


def read_vcf(vcf_file, chr_index=0, pos_index=1, ref_index=3, alt_index=4):
    """Read somatic variants from VCF file"""
    variants = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            
            chrom = cols[chr_index]
            pos = int(cols[pos_index])
            ref = cols[ref_index]
            alt = cols[alt_index]
            
            # Parse AF from INFO field if available
            info = cols[7] if len(cols) > 7 else ""
            af = None
            for field in info.split(';'):
                if field.startswith('AF='):
                    try:
                        af = float(field.split('=')[1])
                    except (ValueError, IndexError):
                        pass
                    break
            
            variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'af': af,
                'loc': f"{cols[chr_index]}:{pos}"
            })
    
    return variants


def count_alleles_at_position(samfile, chrom, pos, ref, alt, chrom_refname=""):
    """Count reference and alternate alleles at a given position"""
    
    # Handle chromosome naming convention
    if "chr" not in chrom and chrom_refname == "chr":
        fetch_chrom = "chr" + chrom
    elif "chr" in chrom and chrom_refname == "":
        if chrom == "chrM":
            fetch_chrom = "MT"  # Handle mitochondrial chromosome
        elif chrom in ["chrX", "chrY"]:
            fetch_chrom = chrom[-1]
        else:
            fetch_chrom = ''.join([c for c in chrom if c.isdigit()])
    else:
        fetch_chrom = chrom
    
    ref_count = 0
    alt_count = 0
    other_count = 0
    
    pysam_pos = pos - 1  # pysam uses 0-based coordinates
    
    try:
        for read in samfile.fetch(fetch_chrom, pysam_pos, pos):
            try:
                index = read.positions.index(pysam_pos)
                base = read.query[index]
                
                if base == ref:
                    ref_count += 1
                elif base == alt:
                    alt_count += 1
                else:
                    other_count += 1
            except (ValueError, IndexError):
                continue  # Position not covered by this read
    except ValueError:
        pass  # Chromosome not in BAM
    
    return ref_count, alt_count, other_count


def check_chromosome_convention(samfile, test_chr='1'):
    """Check if BAM uses 'chr' prefix for chromosome names"""
    try:
        samfile.fetch("chr" + test_chr, 1, 1)
        return "chr"
    except ValueError:
        pass
    
    try:
        samfile.fetch(test_chr, 1, 1)
        return ""
    except ValueError:
        pass
    
    return None


def main():
    parser = argparse.ArgumentParser(
        description="SomaCheck - Evaluate if a sequencing sample is germline or somatic/cancer",
        epilog="Example: SomaCheck.py -s somatic_variants.vcf sample.bam"
    )
    
    parser.add_argument('bam', nargs='?', help='BAM/SAM/CRAM file to check')
    parser.add_argument('-s', '--somatic_vcf', action='store', dest='somatic_vcf', required=True,
                        help='VCF file containing known true somatic variants')
    parser.add_argument('-o', '--output', action='store', dest='output', required=False,
                        default='somacheck_out.txt', help='Output file name [somacheck_out.txt]')
    parser.add_argument('-t', '--threshold', action='store', dest='threshold', type=float,
                        required=False, default=0.05,
                        help='P-value threshold for calling a site as having somatic variants [0.05]')
    parser.add_argument('-m', '--min_reads', action='store', dest='min_reads', type=int,
                        required=False, default=10,
                        help='Minimum read depth at a site to consider it [10]')
    parser.add_argument('-v', '--min_variant_sites', action='store', dest='min_variant_sites',
                        type=int, required=False, default=1,
                        help='Minimum number of sites with variant reads to classify as somatic [1]')
    parser.add_argument('-e', '--expected_vaf', action='store', dest='expected_vaf', type=float,
                        required=False, default=0.3,
                        help='Expected VAF for somatic variants [0.3]')
    parser.add_argument('--output_dir', action='store', required=False, default='.',
                        help='Directory to save output files [./]')
    parser.add_argument('--sanity_check_chr', action='store', dest='sanity_check_chr',
                        required=False, default='1',
                        help="Chromosome to use for BAM index sanity check [1]")
    parser.add_argument('--verbose', action='store_true', dest='verbose', required=False,
                        default=False, help='Print detailed output for each variant site')
    
    args = parser.parse_args()
    
    if args.bam is None:
        eprint('ERROR: No BAM file specified')
        parser.print_help()
        sys.exit(1)
    
    bam = args.bam
    somatic_vcf = args.somatic_vcf
    output = args.output
    threshold = args.threshold
    min_reads = args.min_reads
    min_variant_sites = args.min_variant_sites
    expected_vaf = args.expected_vaf
    output_dir = args.output_dir
    verbose = args.verbose
    
    # Validate inputs
    if not os.path.exists(bam):
        eprint(f'ERROR: BAM file not found: {bam}')
        sys.exit(1)
    
    if not os.path.exists(somatic_vcf):
        eprint(f'ERROR: Somatic VCF file not found: {somatic_vcf}')
        sys.exit(1)
    
    # Open BAM file
    print(strftime("[%Y-%m-%d %H:%M:%S]"), f'Opening BAM file: {bam}')
    
    if bam.endswith('.bam'):
        samfile = pysam.AlignmentFile(bam, 'rb')
    elif bam.endswith('.sam'):
        samfile = pysam.AlignmentFile(bam, 'r')
    elif bam.endswith('.cram'):
        samfile = pysam.AlignmentFile(bam, 'rc')
    else:
        eprint(f'ERROR: Cannot determine file type of {bam}. Use .bam, .sam, or .cram extension.')
        sys.exit(1)
    
    # Check chromosome naming convention
    chrom_refname = check_chromosome_convention(samfile, args.sanity_check_chr)
    if chrom_refname is None:
        eprint(f'ERROR: BAM file does not contain chromosome {args.sanity_check_chr} or chr{args.sanity_check_chr}')
        sys.exit(1)
    
    print(strftime("[%Y-%m-%d %H:%M:%S]"), f'Chromosome naming convention: {"chr prefix" if chrom_refname else "no prefix"}')
    
    # Read somatic variants
    print(strftime("[%Y-%m-%d %H:%M:%S]"), f'Reading somatic variants from: {somatic_vcf}')
    variants = read_vcf(somatic_vcf)
    print(strftime("[%Y-%m-%d %H:%M:%S]"), f'Number of somatic variants: {len(variants)}')
    
    if len(variants) == 0:
        eprint('ERROR: No variants found in the somatic VCF file')
        sys.exit(1)
    
    # Process each variant
    results = []
    total_log_lr = 0.0
    sites_with_coverage = 0
    sites_with_variant = 0
    total_ref = 0
    total_alt = 0
    
    print(strftime("[%Y-%m-%d %H:%M:%S]"), 'Analyzing somatic variant sites...')
    
    for var in variants:
        ref_count, alt_count, other_count = count_alleles_at_position(
            samfile, var['chrom'], var['pos'], var['ref'], var['alt'], chrom_refname
        )
        
        total_reads = ref_count + alt_count
        vaf = calculate_vaf(ref_count, alt_count)
        
        if total_reads >= min_reads:
            sites_with_coverage += 1
            p_value = binomial_test_somatic(ref_count, alt_count, min_reads=min_reads)
            log_lr = log_likelihood_ratio(ref_count, alt_count, expected_vaf)
            total_log_lr += log_lr
            
            if alt_count > 0:
                sites_with_variant += 1
            
            total_ref += ref_count
            total_alt += alt_count
            
            result = {
                'loc': var['loc'],
                'ref': var['ref'],
                'alt': var['alt'],
                'ref_count': ref_count,
                'alt_count': alt_count,
                'total': total_reads,
                'vaf': vaf,
                'p_value': p_value,
                'log_lr': log_lr,
                'somatic_evidence': p_value < threshold and alt_count > 0
            }
            results.append(result)
            
            if verbose:
                status = "SOMATIC" if result['somatic_evidence'] else "GERMLINE"
                print(f"  {var['loc']}: REF={ref_count}, ALT={alt_count}, "
                      f"VAF={vaf:.3f}, p={p_value:.4e}, log_LR={log_lr:.2f} [{status}]")
    
    # Calculate final classification
    print(strftime("[%Y-%m-%d %H:%M:%S]"), 'Computing classification...')
    
    # Write detailed output
    output_path = os.path.join(output_dir, output)
    with open(output_path, 'w') as f:
        f.write("# SomaCheck Results\n")
        f.write(f"# BAM: {bam}\n")
        f.write(f"# Somatic VCF: {somatic_vcf}\n")
        f.write(f"# Threshold: {threshold}\n")
        f.write(f"# Min reads: {min_reads}\n")
        f.write("#\n")
        
        header = ["Location", "Ref", "Alt", "RefCount", "AltCount", "Total", 
                  "VAF", "P_Value", "Log_LR", "SomaticEvidence"]
        f.write("\t".join(header) + "\n")
        
        for r in results:
            line = [r['loc'], r['ref'], r['alt'], str(r['ref_count']), 
                    str(r['alt_count']), str(r['total']), f"{r['vaf']:.4f}",
                    f"{r['p_value']:.4e}", f"{r['log_lr']:.4f}",
                    "TRUE" if r['somatic_evidence'] else "FALSE"]
            f.write("\t".join(line) + "\n")
    
    # Print summary
    print()
    print("=" * 60)
    print("SomaCheck Summary")
    print("=" * 60)
    print(f"BAM file: {bam}")
    print(f"Somatic VCF: {somatic_vcf}")
    print(f"Total somatic variant sites in VCF: {len(variants)}")
    print(f"Sites with sufficient coverage (>={min_reads} reads): {sites_with_coverage}")
    print(f"Sites with any variant reads: {sites_with_variant}")
    
    sites_classified_somatic = sum(1 for r in results if r['somatic_evidence'])
    print(f"Sites classified as having somatic variants: {sites_classified_somatic}")
    
    if sites_with_coverage > 0:
        overall_vaf = total_alt / float(total_ref + total_alt) if (total_ref + total_alt) > 0 else 0
        print(f"Overall VAF across covered sites: {overall_vaf:.4f}")
        print(f"Cumulative log-likelihood ratio: {total_log_lr:.2f}")
        
        # Final classification
        print()
        if sites_classified_somatic >= min_variant_sites:
            classification = "SOMATIC"
            confidence = "HIGH" if sites_classified_somatic >= 3 else "MODERATE"
            print(f"CLASSIFICATION: {classification} (confidence: {confidence})")
            print(f"Interpretation: Sample likely comes from a somatic/cancer source.")
            print(f"               Found {sites_classified_somatic} sites with significant variant evidence.")
        else:
            classification = "GERMLINE"
            if sites_with_variant == 0:
                confidence = "HIGH"
            elif sites_with_variant < min_variant_sites:
                confidence = "MODERATE"
            else:
                confidence = "LOW"
            print(f"CLASSIFICATION: {classification} (confidence: {confidence})")
            print(f"Interpretation: Sample likely comes from a germline (normal) source.")
            print(f"               Only {sites_with_variant} sites showed any variant reads.")
    else:
        print("WARNING: No sites had sufficient coverage for analysis.")
        classification = "UNDETERMINED"
        print(f"CLASSIFICATION: {classification}")
    
    print()
    print(f"Detailed results written to: {output_path}")
    print("=" * 60)
    
    # Write summary file
    summary_path = os.path.join(output_dir, "somacheck_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"BAM\t{bam}\n")
        f.write(f"Somatic_VCF\t{somatic_vcf}\n")
        f.write(f"Total_Variants\t{len(variants)}\n")
        f.write(f"Sites_With_Coverage\t{sites_with_coverage}\n")
        f.write(f"Sites_With_Variant\t{sites_with_variant}\n")
        f.write(f"Sites_Somatic\t{sites_classified_somatic}\n")
        if sites_with_coverage > 0:
            f.write(f"Overall_VAF\t{overall_vaf:.4f}\n")
            f.write(f"Log_LR_Sum\t{total_log_lr:.2f}\n")
        f.write(f"Classification\t{classification}\n")
    
    print(strftime("[%Y-%m-%d %H:%M:%S]"), 'SomaCheck finished.')
    
    return 0 if classification != "UNDETERMINED" else 1


if __name__ == "__main__":
    sys.exit(main())
