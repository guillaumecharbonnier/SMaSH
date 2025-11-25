# Liftover Script: hg38 to hs1

This script automates the process of converting SNP positions from the GRCh38/hg38 reference genome to the T2T-CHM13v2.0 (hs1) reference genome.

## Overview

The script `liftover_hg38_to_hs1.sh` performs the following operations:

1. Downloads the hg38 to hs1 chain file from UCSC Genome Browser
2. Downloads the T2T-CHM13v2.0 (hs1) reference genome
3. Applies Picard LiftoverVcf to convert `snps_GRCh38.vcf` to `snps_hs1.vcf`
4. **Fixes INFO field order** for SMaSH compatibility

### Reference Allele Recovery

The script automatically enables `RECOVER_SWAPPED_REF_ALT=true` to handle variants where the REF and ALT alleles are swapped between assemblies. This is common when lifting over between different reference genomes, as the ancestral vs. derived allele designation may differ. This feature significantly improves the liftover success rate by recovering variants that would otherwise be rejected due to mismatched reference alleles.

### INFO Field Order Fix

**Important:** Picard LiftoverVcf reorders INFO fields alphabetically, but SMaSH expects a specific order (AC, AF, AN, NS, DP, ...) because it extracts the allele frequency from a fixed position. The script automatically reorders the INFO fields after liftover to ensure compatibility with SMaSH.

## Prerequisites

### Required Software

- **Bash** (version 4.0 or higher)
- **Java** (version 8 or higher) - for running Picard
- **Picard Tools** - for the LiftoverVcf operation
- **wget** or **curl** - for downloading files
- **gunzip** - for decompressing files

### System Requirements

- **RAM**: At least 8GB of RAM is recommended (default Java heap size)
- **Disk Space**: Approximately 4GB for the hs1 reference genome and chain files

### Optional Software

- **samtools** - for creating FASTA index (recommended)

### Installing Picard

If Picard is not already installed, you can download it from:
https://github.com/broadinstitute/picard/releases

```bash
# Example installation
wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar
sudo mv picard.jar /usr/local/bin/
```

## Usage

### Basic Usage

Simply run the script from the SMaSH directory:

```bash
./liftover_hg38_to_hs1.sh
```

### What the Script Does

1. **Checks for input file**: Verifies that `snps_GRCh38.vcf` exists
2. **Downloads chain file**: Retrieves `hg38ToHs1.over.chain.gz` if not present
3. **Downloads reference**: Retrieves T2T-CHM13v2.0 reference genome if not present
4. **Creates indices**: Generates necessary index files (`.fai` and `.dict`)
5. **Runs liftover**: Converts the VCF file using Picard LiftoverVcf

### Output Files

After successful execution, the following files will be created:

- `snps_hs1.vcf` - The converted VCF file with hs1 coordinates
- `snps_hs1_rejected.vcf` - Variants that could not be lifted over
- `hg38ToHs1.over.chain.gz` - The chain file (retained for future use)
- `hs1.fa` - The hs1 reference genome (retained for future use)
- `hs1.fa.fai` - FASTA index file
- `hs1.dict` - Sequence dictionary file

## Manual Download (if automatic download fails)

If the script cannot download files automatically due to network restrictions:

### Chain File

Download from: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz
```

### hs1 Reference Genome

Download from: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

Or from NCBI:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
mv GCA_009914755.4_T2T-CHM13v2.0_genomic.fna hs1.fa
```

## Using the Lifted VCF with SMaSH

After creating `snps_hs1.vcf`, you can use it with SMaSH for data aligned to the hs1 reference:

```bash
./SMaSH.py -i snps_hs1.vcf sample1.bam sample2.bam
```

## Troubleshooting

### Picard Not Found

If the script cannot find Picard, either:
1. Install Picard and ensure it's in your PATH
2. Set the `PICARD_JAR` environment variable to point to your picard.jar location:
   ```bash
   export PICARD_JAR=/path/to/picard.jar
   ./liftover_hg38_to_hs1.sh
   ```

### Download Failures

If downloads fail:
1. Check your internet connection
2. Verify that the URLs are accessible from your network
3. Download files manually (see "Manual Download" section above)

### Memory Issues

The script defaults to using 8GB of Java heap memory. If you encounter `OutOfMemoryError`, you can increase the memory allocation:

```bash
export JAVA_MEM=16g
./liftover_hg38_to_hs1.sh
```

Alternatively, for systems with limited RAM, you can reduce the memory:
```bash
export JAVA_MEM=4g
./liftover_hg38_to_hs1.sh
```

**Note:** The T2T-CHM13 reference genome is very large (~3GB), so at least 8GB of heap memory is recommended for successful processing.

### VCF Header Issues

If Picard complains about missing INFO field definitions in the VCF header (e.g., "Key AA found in VariantContext field INFO but this key isn't defined in the VCFHeader"), you need to add proper header definitions.

**Option 1: Use bcftools (recommended if available):**
```bash
# bcftools can automatically infer and add missing header definitions
bcftools view -h snps_GRCh38.vcf > header.txt
# Edit header.txt to add missing INFO definitions, then:
bcftools reheader -h header.txt snps_GRCh38.vcf > snps_GRCh38_fixed.vcf
```

**Option 2: Manually add INFO definitions**
The `snps_GRCh38.vcf` file in this repository includes proper VCF 4.3 headers with all necessary INFO field definitions (AC, AF, AN, NS, DP, EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF, AA).

### Missing Dictionary File

If Picard complains about missing dictionary file, create it manually:

```bash
picard CreateSequenceDictionary R=hs1.fa O=hs1.dict
```

### Mismatched Reference Alleles

If you see many variants with "mismatching reference alleles after lift over", the script automatically handles this by enabling `RECOVER_SWAPPED_REF_ALT=true`. This recovers variants where the REF and ALT alleles are swapped between assemblies, which is common when different reference genomes use different conventions for ancestral vs. derived alleles.

**Note:** Some variants may still fail if there are genuine sequence differences between assemblies that go beyond simple REF/ALT swaps. These will be written to the rejected variants file.

## References

- **T2T Consortium**: https://github.com/marbl/CHM13
- **UCSC Genome Browser**: https://genome.ucsc.edu/
- **Picard Tools**: https://broadinstitute.github.io/picard/
- **T2T-CHM13v2.0 Paper**: Nurk et al. (2022) Science, doi: 10.1126/science.abj6987

## License

This script is part of the SMaSH project and follows the same license terms.
