# YaHS Scaffolding Pipeline

A comprehensive Hi-C scaffolding pipeline using YaHS (Yet another Hi-C scaffolding) for genome assembly improvement. This pipeline integrates preprocessing, scaffolding, and contact map generation in a containerized environment.

## Features

- **Complete Hi-C Scaffolding Workflow**: From raw Hi-C reads to final scaffolded assembly
- **Contact Map Generation**: Supports both Juicer (.hic) and Pretext formats
- **Smart Reruns**: Automatically skips completed steps for efficient pipeline reruns
- **Containerized**: Runs in Singularity container with all dependencies pre-installed
- **Flexible Configuration**: YAML-based configuration for easy customization

## Pipeline Overview

The pipeline consists of three main steps:

1. **Preprocessing & Mapping**: 
   - Index contigs with BWA and samtools
   - Align Hi-C reads to contigs
   - Mark and remove duplicates

2. **Scaffolding with YaHS**:
   - Run YaHS scaffolding algorithm
   - Generate binary files for contact map creation

3. **Contact Map Generation**:
   - Create Juicer Tools compatible .hic files
   - Generate Pretext maps for visualization

## Requirements

- **Singularity** (version 3.0+)
- **Input Files**:
  - Assembled contigs (FASTA format)
  - Hi-C paired-end reads (FASTQ format, can be gzipped)

## Quick Start

### 1. Build the Singularity Container

```bash
# Build the container from the Singularity definition file
singularity build yahs_pipeline.sif Singularity
```

### 2. Prepare Configuration File

Copy and modify the configuration template:

```bash
# Copy the template
cp config_template.yaml my_config.yaml

# Edit the configuration file
nano my_config.yaml
```

### 3. Run the Pipeline

```bash
# Run the complete pipeline
singularity run yahs_pipeline.sif my_config.yaml
```

## Configuration

The pipeline uses a YAML configuration file. Here's an example:

```yaml
# Input files (required)
input:
  contigs: "/path/to/contigs.fa"           # Assembly contigs FASTA file
  
  # Hi-C reads - can specify single pair or multiple pairs
  # Single pair format (for backwards compatibility):
  # hic_r1: "/path/to/hic_R1.fastq.gz"
  # hic_r2: "/path/to/hic_R2.fastq.gz"
  
  # Multiple pairs format (recommended for multiple Hi-C libraries):
  hic_pairs:
    - r1: "/path/to/library1_R1.fastq.gz"
      r2: "/path/to/library1_R2.fastq.gz"
    - r1: "/path/to/library2_R1.fastq.gz"
      r2: "/path/to/library2_R2.fastq.gz"

# Output configuration
output:
  directory: "/path/to/output"             # Output directory
  prefix: "scaffolded"                     # Output files prefix

# Workflow parameters
parameters:
  threads: 16                              # Number of threads for BWA alignment
  java_memory: "32G"                       # Memory for Java applications
  dedup_method: "biobambam2"               # 'biobambam2' or 'picard'
  contact_map_method: "both"               # 'juicer', 'pretext', or 'both'
  
  # YaHS specific parameters
  yahs:
    resolution: "1000,5000,10000"          # Resolution series
    enzyme_motifs: "GATC"                  # Restriction enzyme motifs
    no_contig_ec: false                    # Skip contig error correction
    no_scaffold_ec: false                  # Skip scaffold error correction

# Logging configuration
logging:
  level: "INFO"                            # DEBUG, INFO, WARNING, ERROR
  file: "yahs_workflow.log"                # Log file name
```

## Multiple Hi-C Libraries Support

The pipeline supports processing multiple Hi-C libraries, which is beneficial for:

- **Improved scaffolding quality**: Multiple libraries provide more contact information
- **Different restriction enzymes**: Combining data from different Hi-C protocols
- **Sequencing depth**: Merging multiple sequencing runs for better coverage

### Configuration Formats

**Single Hi-C Library (backwards compatible):**
```yaml
input:
  contigs: "/path/to/contigs.fa"
  hic_r1: "/path/to/hic_R1.fastq.gz"
  hic_r2: "/path/to/hic_R2.fastq.gz"
```

**Multiple Hi-C Libraries (recommended):**
```yaml
input:
  contigs: "/path/to/contigs.fa"
  hic_pairs:
    - r1: "/path/to/DpnII_R1.fastq.gz"
      r2: "/path/to/DpnII_R2.fastq.gz"
    - r1: "/path/to/HindIII_R1.fastq.gz"
      r2: "/path/to/HindIII_R2.fastq.gz"
    - r1: "/path/to/Arima_R1.fastq.gz"
      r2: "/path/to/Arima_R2.fastq.gz"
```

### Processing Workflow

When multiple Hi-C pairs are specified:
1. Each pair is aligned independently using BWA
2. All aligned BAM files are merged with proper name sorting
3. Duplicate marking is performed on the merged dataset
4. Scaffolding uses the combined Hi-C information

## Usage Examples

### Basic Usage

```bash
# Run with default settings
singularity run yahs_pipeline.sif config.yaml
```

### Advanced Usage

```bash
# Run with specific Singularity options
singularity run --bind /data:/data --cleanenv yahs_pipeline.sif config.yaml

# Run individual tools within the container
singularity exec yahs_pipeline.sif yahs --help
singularity exec yahs_pipeline.sif samtools --version
```

### Debugging

```bash
# Run with debug logging
# Set logging level to DEBUG in config.yaml

# Interactive shell in container
singularity shell yahs_pipeline.sif
```

## Output Files

The pipeline generates several output directories:

```
output_directory/
├── preprocessing/
│   ├── hic_name_sorted.bam
│   ├── hic_dedup.bam
│   └── hic_dedup.bed (optional)
├── scaffolding/
│   ├── scaffolded_scaffolds_final.fa
│   ├── scaffolded_scaffolds_final.agp
│   ├── scaffolded.bin
│   └── scaffolded_*.agp (intermediate files)
├── contact_maps/
│   ├── scaffolded.hic (Juicer format)
│   ├── scaffolded.pretext (Pretext format)
│   └── pretext_snapshots/ (PNG images)
└── logs/
    └── yahs_workflow.log
```

## Key Output Files

- **`*_scaffolds_final.fa`**: Final scaffolded assembly
- **`*_scaffolds_final.agp`**: AGP file describing scaffold structure
- **`*.hic`**: Hi-C contact map for Juicebox visualization
- **`*.pretext`**: Pretext map for PretextView visualization

## Performance Optimization

### Memory Requirements

- **Minimum**: 16GB RAM
- **Recommended**: 32GB+ RAM for large genomes
- **Java heap**: Configure `java_memory` parameter based on available RAM

### CPU Usage

- Set `threads` parameter to match available CPU cores
- BWA alignment is the most CPU-intensive step

### Storage

- Temporary files can be large (several GB)
- Ensure sufficient disk space in output directory

## Container Tools

The Singularity container includes:

- **YaHS**: Hi-C scaffolding algorithm
- **BWA**: Sequence alignment
- **samtools**: SAM/BAM manipulation
- **bedtools**: Genome arithmetic
- **biobambam2**: BAM processing
- **Juicer Tools**: Contact map generation
- **Pretext Suite**: Contact map visualization
- **Python 3** with PyYAML

## Troubleshooting

### Common Issues

1. **Missing input files**: Verify file paths in configuration
2. **Insufficient memory**: Increase `java_memory` parameter
3. **Permission errors**: Use `--bind` option for data directories
4. **Contact map generation fails**: Check YaHS output files exist

### Log Analysis

```bash
# Check the main log file
tail -f output_directory/yahs_workflow.log

# Check specific step logs
grep "ERROR" output_directory/yahs_workflow.log
```

### Container Debugging

```bash
# List available tools
singularity exec yahs_pipeline.sif ls /usr/local/bin

# Check tool versions
singularity exec yahs_pipeline.sif yahs --version
singularity exec yahs_pipeline.sif java -version
```

## Smart Reruns

The pipeline automatically detects completed steps and skips them:

- **Step 1**: Checks for `hic_dedup.bam`
- **Step 2**: Checks for `*_scaffolds_final.fa`
- **Step 3**: Checks for contact map files individually

This allows for efficient reruns if a step fails or you want to regenerate only specific outputs.

## Citation

If you use this pipeline, please cite:

- **YaHS**: Zhou, C., McCarthy, S.A., Durbin, R. YaHS: yet another Hi-C scaffolding tool. Bioinformatics 39, 6 (2023).
- **BWA**: Li, H. and Durbin, R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760 (2009).
- **Juicer**: Durand, N.C. et al. Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. Cell Systems 3, 95-98 (2016).
- **Pretext**: https://github.com/wtsi-hpag/PretextView

## License

This pipeline is released under the MIT License. See individual tools for their respective licenses.

## Support

For issues and questions:
1. Check the troubleshooting section
2. Review log files for error messages
3. Consult the YaHS documentation: https://github.com/c-zhou/yahs