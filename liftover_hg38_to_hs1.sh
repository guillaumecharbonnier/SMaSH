#!/bin/bash

# Script to retrieve hg38 to hs1 chain file and apply Picard liftOverVcf
# to convert snps_GRCh38.vcf to snps_hs1.vcf
#
# Author: SMaSH Team
# Date: 2024

set -euo pipefail

# Variables
# You can set PICARD_JAR environment variable before running this script
# to specify the location of picard.jar, e.g.:
# export PICARD_JAR=/path/to/picard.jar
# You can also set JAVA_MEM to control Java heap memory (default: 8g)
# export JAVA_MEM=16g
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHAIN_FILE="hg38ToHs1.over.chain.gz"
INPUT_VCF="${SCRIPT_DIR}/snps_GRCh38.vcf"
OUTPUT_VCF="${SCRIPT_DIR}/snps_hs1.vcf"
REJECTED_VCF="${SCRIPT_DIR}/snps_hs1_rejected.vcf"
HS1_REFERENCE="hs1.fa"
HS1_REFERENCE_GZ="hs1.fa.gz"
PICARD_JAR="${PICARD_JAR:-}"
JAVA_MEM="${JAVA_MEM:-8g}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

echo_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

echo_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to download chain file
download_chain_file() {
    echo_info "Downloading hg38 to hs1 chain file..."
    
    # Try multiple sources for the chain file
    local chain_urls=(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
        "https://hgdownload.gi.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
        "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
    )
    
    for url in "${chain_urls[@]}"; do
        echo_info "Trying to download from: $url"
        if wget -q --show-progress -O "${CHAIN_FILE}" "$url" 2>/dev/null || curl -L -o "${CHAIN_FILE}" "$url" 2>/dev/null; then
            if [ -f "${CHAIN_FILE}" ] && [ -s "${CHAIN_FILE}" ]; then
                echo_info "Successfully downloaded chain file from $url"
                return 0
            fi
        fi
        echo_warn "Failed to download from $url, trying next source..."
    done
    
    echo_error "Failed to download chain file from all sources."
    echo_error "Please manually download hg38ToHs1.over.chain.gz from UCSC and place it in the script directory."
    echo_error "URL: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz"
    return 1
}

# Function to download hs1 reference genome
download_hs1_reference() {
    echo_info "Downloading hs1 reference genome..."
    
    # Try multiple sources for the reference
    local ref_urls=(
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
    )
    
    for url in "${ref_urls[@]}"; do
        echo_info "Trying to download from: $url"
        if wget -q --show-progress -O "${HS1_REFERENCE_GZ}" "$url" 2>/dev/null || curl -L -o "${HS1_REFERENCE_GZ}" "$url" 2>/dev/null; then
            if [ -f "${HS1_REFERENCE_GZ}" ] && [ -s "${HS1_REFERENCE_GZ}" ]; then
                echo_info "Successfully downloaded reference genome"
                echo_info "Decompressing reference genome..."
                gunzip -f "${HS1_REFERENCE_GZ}"
                
                # Create index if needed
                if command -v samtools &> /dev/null; then
                    echo_info "Creating FASTA index..."
                    samtools faidx "${HS1_REFERENCE}"
                fi
                
                return 0
            fi
        fi
        echo_warn "Failed to download from $url, trying next source..."
    done
    
    echo_error "Failed to download hs1 reference genome from all sources."
    echo_error "Please manually download the T2T-CHM13v2.0 reference genome."
    return 1
}

# Function to create sequence dictionary
create_sequence_dict() {
    local reference=$1
    local dict_file="${reference%.fa}.dict"
    
    if [ ! -f "${dict_file}" ]; then
        echo_info "Creating sequence dictionary..."
        if [ -n "${PICARD_JAR}" ] && [ -f "${PICARD_JAR}" ]; then
            java -Xmx${JAVA_MEM} -jar "${PICARD_JAR}" CreateSequenceDictionary \
                R="${reference}" \
                O="${dict_file}"
        elif command -v picard &> /dev/null; then
            picard -Xmx${JAVA_MEM} CreateSequenceDictionary \
                R="${reference}" \
                O="${dict_file}"
        else
            echo_warn "Picard not found. Dictionary creation skipped."
            echo_warn "You may need to create it manually before running LiftoverVcf."
            return 1
        fi
    fi
    return 0
}

# Function to find Picard
find_picard() {
    echo_info "Looking for Picard..."
    
    # Check if PICARD_JAR environment variable is set
    if [ -n "${PICARD_JAR}" ] && [ -f "${PICARD_JAR}" ]; then
        echo_info "Using Picard from environment variable: ${PICARD_JAR}"
        return 0
    fi
    
    # Check if picard is in PATH
    if command -v picard &> /dev/null; then
        echo_info "Found picard in PATH"
        return 0
    fi
    
    # Common Picard locations
    local picard_locations=(
        "/usr/local/bin/picard.jar"
        "/usr/share/picard/picard.jar"
        "/opt/picard/picard.jar"
        "$HOME/picard.jar"
        "/usr/share/java/picard.jar"
        "/usr/local/share/picard/picard.jar"
    )
    
    for location in "${picard_locations[@]}"; do
        if [ -f "$location" ]; then
            PICARD_JAR="$location"
            echo_info "Found Picard at: $PICARD_JAR"
            return 0
        fi
    done
    
    echo_warn "Picard not found in common locations."
    echo_warn "You can download it from: https://github.com/broadinstitute/picard/releases"
    echo_warn "Or set PICARD_JAR environment variable to the path of picard.jar"
    return 1
}

# Function to run liftover
run_liftover() {
    echo_info "Running Picard LiftoverVcf with ${JAVA_MEM} heap memory..."
    
    local cmd=""
    if [ -n "${PICARD_JAR}" ] && [ -f "${PICARD_JAR}" ]; then
        cmd="java -Xmx${JAVA_MEM} -jar ${PICARD_JAR}"
    elif command -v picard &> /dev/null; then
        cmd="picard -Xmx${JAVA_MEM}"
    else
        echo_error "Picard not found. Cannot run liftover."
        return 1
    fi
    
    # Uncompress chain file if it's compressed
    local chain_uncompressed="${CHAIN_FILE%.gz}"
    if [[ "${CHAIN_FILE}" == *.gz ]]; then
        if [ ! -f "${chain_uncompressed}" ]; then
            echo_info "Decompressing chain file..."
            gunzip -c "${CHAIN_FILE}" > "${chain_uncompressed}"
        fi
    else
        chain_uncompressed="${CHAIN_FILE}"
    fi
    
    if $cmd LiftoverVcf \
        I="${INPUT_VCF}" \
        O="${OUTPUT_VCF}" \
        CHAIN="${chain_uncompressed}" \
        REJECT="${REJECTED_VCF}" \
        R="${HS1_REFERENCE}" \
        WARN_ON_MISSING_CONTIG=true \
        RECOVER_SWAPPED_REF_ALT=true; then
        echo_info "Liftover completed successfully!"
        echo_info "Output VCF: ${OUTPUT_VCF}"
        echo_info "Rejected variants: ${REJECTED_VCF}"
        return 0
    else
        echo_error "Liftover failed!"
        return 1
    fi
}

# Function to fix INFO field order in the output VCF
# Picard LiftoverVcf reorders INFO fields alphabetically, but SMaSH expects
# a specific order (AC, AF, AN, ...) because it uses INFO[1] for AF value
fix_info_field_order() {
    echo_info "Fixing INFO field order for SMaSH compatibility..."
    
    local temp_file="${OUTPUT_VCF}.tmp"
    
    # Use awk to reorder INFO fields
    awk 'BEGIN {FS=OFS="\t"}
    /^#/ {print; next}
    {
        # Parse INFO field
        n = split($8, fields, ";")
        
        # Build associative array and collect flags
        delete info
        flags = ""
        for (i=1; i<=n; i++) {
            if (index(fields[i], "=") > 0) {
                split(fields[i], kv, "=")
                info[kv[1]] = kv[2]
            } else {
                # Flag field (no value)
                flags = flags ";" fields[i]
            }
        }
        
        # Build output in expected order: AC, AF, AN, NS, DP, EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF, AA
        output = ""
        order[1]="AC"; order[2]="AF"; order[3]="AN"; order[4]="NS"; order[5]="DP"
        order[6]="EAS_AF"; order[7]="AMR_AF"; order[8]="AFR_AF"; order[9]="EUR_AF"; order[10]="SAS_AF"
        order[11]="AA"
        
        for (i=1; i<=11; i++) {
            key = order[i]
            if (key in info) {
                if (output != "") output = output ";"
                output = output key "=" info[key]
            }
        }
        
        # Append flags
        if (flags != "") output = output flags
        
        $8 = output
        print
    }' "${OUTPUT_VCF}" > "${temp_file}"
    
    if [ $? -eq 0 ]; then
        mv "${temp_file}" "${OUTPUT_VCF}"
        echo_info "INFO field order fixed successfully."
        return 0
    else
        echo_error "Failed to fix INFO field order."
        rm -f "${temp_file}"
        return 1
    fi
}

# Main execution
main() {
    echo_info "Starting hg38 to hs1 liftover process..."
    echo_info "Working directory: ${SCRIPT_DIR}"
    
    # Check if input VCF exists
    if [ ! -f "${INPUT_VCF}" ]; then
        echo_error "Input VCF file not found: ${INPUT_VCF}"
        exit 1
    fi
    
    # Find Picard
    find_picard
    
    # Download chain file if not present
    if [ ! -f "${CHAIN_FILE}" ]; then
        if ! download_chain_file; then
            echo_error "Cannot proceed without chain file."
            exit 1
        fi
    else
        echo_info "Chain file already exists: ${CHAIN_FILE}"
    fi
    
    # Download reference if not present
    if [ ! -f "${HS1_REFERENCE}" ]; then
        if ! download_hs1_reference; then
            echo_error "Cannot proceed without reference genome."
            exit 1
        fi
    else
        echo_info "Reference genome already exists: ${HS1_REFERENCE}"
    fi
    
    # Create sequence dictionary if needed
    create_sequence_dict "${HS1_REFERENCE}"
    
    # Run liftover
    if run_liftover; then
        # Fix INFO field order for SMaSH compatibility
        fix_info_field_order
        echo_info "Process completed successfully!"
        exit 0
    else
        echo_error "Process failed!"
        exit 1
    fi
}

# Run main function
main
