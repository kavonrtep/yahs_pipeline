#!/usr/bin/env python3
"""
YaHS Scaffolding Workflow Script
Performs Hi-C scaffolding using YaHS with preprocessing and contact map generation
"""

import os
import sys
import yaml
import subprocess
import logging
import argparse
from pathlib import Path


class YaHSWorkflow:
    def __init__(self, config_file):
        """Initialize workflow with configuration"""
        self.config = self.load_config(config_file)
        self.setup_logging()
        self.validate_config()

    def load_config(self, config_file):
        """Load YAML configuration file"""
        try:
            with open(config_file, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            print(f"Error loading config file {config_file}: {e}")
            sys.exit(1)

    def setup_logging(self):
        """Setup logging configuration"""
        log_level = self.config.get('logging', {}).get('level', 'INFO')
        log_file = self.config.get('logging', {}).get('file', 'yahs_workflow.log')

        # Force DEBUG level if we're troubleshooting contact maps
        if log_level.upper() == 'INFO':
            log_level = 'DEBUG'

        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)]
            )
        self.logger = logging.getLogger(__name__)

    def validate_config(self):
        """Validate required configuration parameters"""
        required_keys = ['input', 'output', 'parameters']
        for key in required_keys:
            if key not in self.config:
                self.logger.error(f"Missing required config section: {key}")
                sys.exit(1)

        # Check required input files
        input_config = self.config['input']
        
        # Check contigs file
        if 'contigs' not in input_config:
            self.logger.error("Missing required input: contigs")
            sys.exit(1)
        if not os.path.exists(input_config['contigs']):
            self.logger.error(f"Contigs file not found: {input_config['contigs']}")
            sys.exit(1)
        
        # Check Hi-C files - support both single pair and multiple pairs format
        if 'hic_pairs' in input_config:
            # Multiple pairs format
            if not isinstance(input_config['hic_pairs'], list) or len(input_config['hic_pairs']) == 0:
                self.logger.error("hic_pairs must be a non-empty list")
                sys.exit(1)
            
            for i, pair in enumerate(input_config['hic_pairs']):
                if 'r1' not in pair or 'r2' not in pair:
                    self.logger.error(f"Hi-C pair {i+1} missing r1 or r2 file")
                    sys.exit(1)
                if not os.path.exists(pair['r1']):
                    self.logger.error(f"Hi-C R1 file not found: {pair['r1']}")
                    sys.exit(1)
                if not os.path.exists(pair['r2']):
                    self.logger.error(f"Hi-C R2 file not found: {pair['r2']}")
                    sys.exit(1)
                    
        elif 'hic_r1' in input_config and 'hic_r2' in input_config:
            # Single pair format (backwards compatibility)
            if not os.path.exists(input_config['hic_r1']):
                self.logger.error(f"Hi-C R1 file not found: {input_config['hic_r1']}")
                sys.exit(1)
            if not os.path.exists(input_config['hic_r2']):
                self.logger.error(f"Hi-C R2 file not found: {input_config['hic_r2']}")
                sys.exit(1)
        else:
            self.logger.error("No Hi-C input files specified. Use either 'hic_pairs' or 'hic_r1'/'hic_r2' format")
            sys.exit(1)

    def run_command(self, cmd, description):
        """Execute shell command with logging"""
        self.logger.info(f"Running: {description}")
        self.logger.debug(f"Command: {cmd}")

        try:
            result = subprocess.run(
                cmd, shell=True, check=True, capture_output=True, text=True
                )
            if result.stdout:
                self.logger.debug(f"STDOUT: {result.stdout}")
            return result
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {description}")
            self.logger.error(f"Exit code: {e.returncode}")
            self.logger.error(f"STDERR: {e.stderr}")
            sys.exit(1)

    def create_output_directory(self):
        """Create output directory structure"""
        output_dir = Path(self.config['output']['directory'])
        output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        subdirs = ['preprocessing', 'scaffolding', 'contact_maps', 'logs']
        for subdir in subdirs:
            (output_dir / subdir).mkdir(exist_ok=True)

        self.output_dir = output_dir
        self.preprocessing_dir = output_dir / 'preprocessing'
        self.scaffolding_dir = output_dir / 'scaffolding'
        self.contact_maps_dir = output_dir / 'contact_maps'

    def step1_preprocessing(self):
        """Step 1: Preprocessing & Mapping"""
        self.logger.info("=== Step 1: Preprocessing & Mapping ===")
        
        # Check if preprocessing is already completed
        hic_dedup_check = self.preprocessing_dir / "hic_dedup.bam"
        if os.path.exists(hic_dedup_check):
            self.logger.info("Preprocessing already completed, skipping step 1")
            self.hic_dedup_bam = hic_dedup_check
            
            # Still log information about Hi-C pairs for transparency
            input_config = self.config['input']
            if 'hic_pairs' in input_config:
                self.logger.info(f"Configuration specifies {len(input_config['hic_pairs'])} Hi-C library pairs")
            else:
                self.logger.info("Configuration specifies single Hi-C library pair")
            return

        contigs_fa = self.config['input']['contigs']
        threads = self.config['parameters'].get('threads', 8)
        input_config = self.config['input']
        
        # Prepare Hi-C file pairs
        hic_pairs = []
        if 'hic_pairs' in input_config:
            # Multiple pairs format
            hic_pairs = input_config['hic_pairs']
            self.logger.info(f"Processing {len(hic_pairs)} Hi-C library pairs")
        else:
            # Single pair format (backwards compatibility)
            hic_pairs = [{'r1': input_config['hic_r1'], 'r2': input_config['hic_r2']}]
            self.logger.info("Processing single Hi-C library pair")

        # Index contigs
        self.logger.info("Indexing contigs...")
        cmd = f"samtools faidx {contigs_fa}"
        self.run_command(cmd, "Index contigs with samtools")

        # create BWA index if not exists
        bwa_index = contigs_fa + '.bwt'
        if not os.path.exists(bwa_index):
            self.logger.info("Creating BWA index for contigs...")
            cmd = f"bwa index {contigs_fa}"
            self.run_command(cmd, "Create BWA index for contigs")

        # Align Hi-C reads (potentially multiple pairs)
        self.logger.info("Aligning Hi-C reads...")
        
        aligned_bams = []
        for i, pair in enumerate(hic_pairs):
            pair_name = f"pair_{i+1:02d}"
            self.logger.info(f"Processing Hi-C pair {i+1}/{len(hic_pairs)}: {pair['r1']}, {pair['r2']}")
            
            hic_aligned = self.preprocessing_dir / f"hic_{pair_name}_aligned.bam"
            cmd = f"""bwa mem -t {threads} {contigs_fa} {pair['r1']} {pair['r2']} | \\
                     samtools view -bSh - | \\
                     samtools sort -n -O BAM -o {hic_aligned}"""
            self.run_command(cmd, f"Align Hi-C pair {i+1} with BWA and sort by name")
            aligned_bams.append(str(hic_aligned))
        
        # Merge multiple BAM files if necessary
        hic_name_sorted = self.preprocessing_dir / "hic_name_sorted.bam"
        if len(aligned_bams) == 1:
            # Single pair - just rename
            cmd = f"mv {aligned_bams[0]} {hic_name_sorted}"
            self.run_command(cmd, "Rename single aligned BAM file")
        else:
            # Multiple pairs - merge them
            self.logger.info(f"Merging {len(aligned_bams)} aligned BAM files...")
            bam_list = " ".join(aligned_bams)
            cmd = f"samtools merge -n {hic_name_sorted} {bam_list}"
            self.run_command(cmd, "Merge multiple Hi-C aligned BAM files")
            
            # Clean up individual BAM files
            for bam in aligned_bams:
                cmd = f"rm {bam}"
                self.run_command(cmd, f"Remove temporary BAM file {bam}")

        # Mark duplicates
        self.logger.info("Marking duplicates...")
        hic_dedup = self.preprocessing_dir / "hic_dedup.bam"
        dedup_method = self.config['parameters'].get('dedup_method', 'biobambam2')

        if dedup_method == 'biobambam2':
            cmd = f"bammarkduplicates2 I={hic_name_sorted} O={hic_dedup}"
            self.run_command(cmd, "Mark duplicates with biobambam2")
        elif dedup_method == 'picard':
            metrics_file = self.preprocessing_dir / "metrics.txt"
            cmd = f"""java -jar picard.jar MarkDuplicates \\
                     I={hic_name_sorted} O={hic_dedup} M={metrics_file}"""
            self.run_command(cmd, "Mark duplicates with Picard")
        elif dedup_method == 'sambamba':
            threads = self.config['parameters'].get('threads', 8)
            cmd = (f"sambamba markdup -t {threads} --show-progress {hic_name_sorted}"
                   f" {hic_dedup}")
            self.run_command(cmd, "Mark duplicates with sambamba")
        else:
            self.logger.error(f"Unknown dedup_method: {dedup_method}. Supported methods: biobambam2, picard, sambamba")
            sys.exit(1)

        # Optional: Convert to BED
        if self.config['parameters'].get('generate_bed', False):
            self.logger.info("Converting BAM to BED...")
            hic_bed = self.preprocessing_dir / "hic_dedup.bed"
            cmd = f"bedtools bamtobed -i {hic_dedup} > {hic_bed}"
            self.run_command(cmd, "Convert BAM to BED format")

        self.hic_dedup_bam = hic_dedup

    def step2_scaffolding(self):
        """Step 2: Scaffolding with YaHS"""
        self.logger.info("=== Step 2: Scaffolding with YaHS ===")
        
        # Check if scaffolding is already completed
        output_prefix = self.scaffolding_dir / self.config['output'].get('prefix', 'scaffolded')
        final_scaffold = f"{output_prefix}_scaffolds_final.fa"
        if os.path.exists(final_scaffold):
            self.logger.info("Scaffolding already completed, skipping step 2")
            self.output_prefix = output_prefix
            return

        contigs_fa = self.config['input']['contigs']
        output_prefix = self.scaffolding_dir / self.config['output'].get(
            'prefix', 'scaffolded'
            )

        # Build YaHS command with --make-bin for contact map generation
        cmd_parts = ["yahs", f"-o {output_prefix}", "--make-bin", contigs_fa, str(self.hic_dedup_bam)]

        # Add optional parameters
        yahs_params = self.config['parameters'].get('yahs', {})

        if 'resolution' in yahs_params:
            cmd_parts.insert(1, f"-r {yahs_params['resolution']}")

        if 'agp_file' in yahs_params:
            cmd_parts.insert(1, f"-a {yahs_params['agp_file']}")

        if yahs_params.get('no_contig_ec', False):
            cmd_parts.insert(1, "--no-contig-ec")

        if yahs_params.get('no_scaffold_ec', False):
            cmd_parts.insert(1, "--no-scaffold-ec")

        if 'enzyme_motifs' in yahs_params:
            cmd_parts.insert(1, f"-e {yahs_params['enzyme_motifs']}")

        cmd = " ".join(cmd_parts)
        self.run_command(cmd, "Run YaHS scaffolding")

        self.output_prefix = output_prefix

    def step3_contact_maps(self):
        """Step 3: Generate Hi-C Contact Maps"""
        self.logger.info("=== Step 3: Generate Hi-C Contact Maps ===")
        
        # Check if contact maps are already generated
        contact_map_method = self.config['parameters'].get('contact_map_method', 'both')
        hic_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.hic"
        pretext_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.pretext"
        
        skip_juicer = False
        skip_pretext = False
        
        if contact_map_method in ['juicer', 'both'] and os.path.exists(hic_file):
            self.logger.info("Juicer .hic file already exists, skipping Juicer generation")
            skip_juicer = True
            
        if contact_map_method in ['pretext', 'both'] and os.path.exists(pretext_file):
            self.logger.info("Pretext map already exists, skipping Pretext generation")
            skip_pretext = True
            
        if skip_juicer and skip_pretext:
            self.logger.info("All contact maps already generated, skipping step 3")
            return

        if contact_map_method in ['juicer', 'both'] and not skip_juicer:
            self.generate_juicer_hic()

        if contact_map_method in ['pretext', 'both'] and not skip_pretext:
            self.generate_pretext_maps()

    def generate_juicer_hic(self):
        """Generate .hic file using Juicer Tools"""
        self.logger.info("Generating .hic file with Juicer Tools...")

        # Required files for Juicer .hic generation
        bin_file = f"{self.output_prefix}.bin"
        agp_file = f"{self.output_prefix}_scaffolds_final.agp"
        contigs_fa = self.config['input']['contigs']
        contigs_fai = f"{contigs_fa}.fai"
        scaffold_fa = f"{self.output_prefix}_scaffolds_final.fa"
        hic_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.hic"
        
        # Check if required files exist
        if not os.path.exists(bin_file):
            self.logger.error(f"YaHS bin file not found: {bin_file}")
            self.logger.error("Cannot generate Juicer .hic file without YaHS bin file")
            return
            
        if not os.path.exists(agp_file):
            self.logger.error(f"YaHS AGP file not found: {agp_file}")
            self.logger.error("Cannot generate Juicer .hic file without AGP file")
            return
            
        if not os.path.exists(contigs_fai):
            self.logger.info(f"Creating contig index file: {contigs_fai}")
            cmd = f"samtools faidx {contigs_fa}"
            self.run_command(cmd, "Create contig index file")
            
        # Create chromosome sizes file for scaffolds
        scaffold_sizes = f"{scaffold_fa}.chrom.sizes"
        if not os.path.exists(scaffold_sizes):
            self.logger.info("Creating scaffold chromosome sizes file...")
            if not os.path.exists(f"{scaffold_fa}.fai"):
                cmd = f"samtools faidx {scaffold_fa}"
                self.run_command(cmd, "Create scaffold index file")
            cmd = f"cut -f1,2 {scaffold_fa}.fai > {scaffold_sizes}"
            self.run_command(cmd, "Create chromosome sizes file")

        # Step 1: Convert YaHS bin file to Juicer format using juicer pre
        alignments_file = self.contact_maps_dir / "alignments_sorted.txt"
        alignments_raw = self.contact_maps_dir / "alignments_raw.txt"
        threads = self.config['parameters'].get('threads', 8)
        memory_sort = self.config['parameters'].get('sort_memory', '32G')
        
        self.logger.info("Converting YaHS bin file to Juicer format...")
        self.logger.debug(f"Using bin file: {bin_file}")
        self.logger.debug(f"Using AGP file: {agp_file}")
        self.logger.debug(f"Using contig index: {contigs_fai}")
        
        # First try to run juicer pre without piping to see raw output
        cmd = f"juicer pre {bin_file} {agp_file} {contigs_fai} > {alignments_raw}"
        self.run_command(cmd, "Run juicer pre to generate raw alignments")
        
        # Check and debug the raw output
        if os.path.exists(alignments_raw):
            file_size = os.path.getsize(alignments_raw)
            self.logger.info(f"Raw alignments file size: {file_size} bytes")
            
            # Read first few lines to debug format
            with open(alignments_raw, 'r') as f:
                lines = [f.readline().strip() for _ in range(10)]
                self.logger.debug("First 10 lines of raw alignments:")
                for i, line in enumerate(lines, 1):
                    self.logger.debug(f"Line {i}: {line}")
            
            # Filter out status messages and keep only data lines
            self.logger.info("Filtering and sorting alignment data...")
            cmd = f"""grep -v "Writing\\|Skipping\\|Start\\|Not including" {alignments_raw} | \\
                     grep -E "^[^\\s]" | \\
                     sort -k2,2d -k6,6d -T {self.contact_maps_dir} --parallel={threads} -S{memory_sort} | \\
                     awk 'NF' > {alignments_file}"""
            
            self.run_command(cmd, "Filter and sort alignment data")
        else:
            self.logger.error(f"Raw alignments file not created: {alignments_raw}")
            return
        
        # Check if final alignments file was created and is not empty
        if not os.path.exists(alignments_file):
            self.logger.error("Failed to create alignments file")
            return
            
        final_size = os.path.getsize(alignments_file)
        self.logger.info(f"Final alignments file size: {final_size} bytes")
        
        if final_size == 0:
            self.logger.error("Alignments file is empty after filtering")
            self.logger.info("Falling back to BAM-based contact map generation...")
            self.generate_juicer_from_bam()
            return

        # Step 2: Generate .hic file with Juicer Tools
        memory = self.config['parameters'].get('java_memory', '32G')
        juicer_jar = self.config['parameters'].get(
            'juicer_tools_jar', '/usr/local/bin/juicer_tools.jar'
            )

        self.logger.info("Generating .hic file with Juicer Tools...")
        cmd = f"""java -Xmx{memory} -jar {juicer_jar} pre \\
                 {alignments_file} {hic_file} {scaffold_sizes}"""

        self.run_command(cmd, "Generate .hic file with Juicer Tools")

    def generate_juicer_from_bam(self):
        """Fallback method: Generate .hic file directly from BAM file"""
        self.logger.info("Generating .hic file from BAM file...")
        
        # Create Juicer-compatible alignment file from BAM
        alignments_file = self.contact_maps_dir / "alignments_from_bam.txt"
        scaffold_fa = f"{self.output_prefix}_scaffolds_final.fa"
        scaffold_sizes = f"{scaffold_fa}.chrom.sizes"
        hic_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.hic"
        
        self.logger.info("Converting BAM to Juicer format...")
        # Convert BAM to Juicer short format
        cmd = f"""samtools view {self.hic_dedup_bam} | \\
                 awk 'BEGIN{{OFS="\\t"}} \\
                 {{ \\
                     if (and($2,0x40)) strand1=0; else strand1=1; \\
                     if (and($2,0x10)) strand1=1-strand1; \\
                     if (and($2,0x20)) strand2=1; else strand2=0; \\
                     if (and($2,0x80)) strand2=1-strand2; \\
                     print $1, strand1, $3, $4, 1, strand2, $7, $8, 1, 1 \\
                 }}' | \\
                 sort -k3,3d -k7,7d -T {self.contact_maps_dir} > {alignments_file}"""
        
        self.run_command(cmd, "Convert BAM to Juicer format")
        
        # Check if alignments file was created
        if not os.path.exists(alignments_file) or os.path.getsize(alignments_file) == 0:
            self.logger.error("Failed to create alignments from BAM file")
            return
            
        # Generate .hic file
        memory = self.config['parameters'].get('java_memory', '32G')
        juicer_jar = self.config['parameters'].get(
            'juicer_tools_jar', '/usr/local/bin/juicer_tools.jar'
            )

        self.logger.info("Generating .hic file from BAM-derived alignments...")
        cmd = f"""java -Xmx{memory} -jar {juicer_jar} pre \\
                 {alignments_file} {hic_file} {scaffold_sizes}"""

        self.run_command(cmd, "Generate .hic file from BAM alignments")

    def generate_pretext_maps(self):
        """Generate Pretext maps"""
        self.logger.info("Generating Pretext maps...")

        pretext_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.pretext"
        
        # Use the dedup BAM file directly for Pretext map generation
        if not os.path.exists(self.hic_dedup_bam):
            self.logger.error(f"Dedup BAM file not found: {self.hic_dedup_bam}")
            self.logger.error("Cannot generate Pretext maps")
            return

        # Generate pretext map from BAM file
        cmd = f"samtools view {self.hic_dedup_bam} | PretextMap -o {pretext_file}"

        self.run_command(cmd, "Generate Pretext map")

        # Generate snapshots
        snapshot_dir = self.contact_maps_dir / "pretext_snapshots"
        snapshot_dir.mkdir(exist_ok=True)

        cmd = f'PretextSnapshot -m {pretext_file} --sequences "=full" -o {snapshot_dir}'
        self.run_command(cmd, "Generate Pretext snapshots")

    def run_workflow(self):
        """Run the complete YaHS workflow"""
        self.logger.info("Starting YaHS scaffolding workflow...")

        try:
            self.create_output_directory()
            self.step1_preprocessing()
            self.step2_scaffolding()
            self.step3_contact_maps()

            self.logger.info("YaHS workflow completed successfully!")
            self.logger.info(f"Results are available in: {self.output_dir}")

        except Exception as e:
            self.logger.error(f"Workflow failed: {e}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description='YaHS Scaffolding Workflow')
    parser.add_argument('config', help='YAML configuration file')
    parser.add_argument('--version', action='version', version='YaHS Workflow v1.0')

    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Configuration file not found: {args.config}")
        sys.exit(1)

    workflow = YaHSWorkflow(args.config)
    workflow.run_workflow()


if __name__ == "__main__":
    main()