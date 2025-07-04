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
        required_inputs = ['contigs', 'hic_r1', 'hic_r2']
        for inp in required_inputs:
            if inp not in input_config:
                self.logger.error(f"Missing required input: {inp}")
                sys.exit(1)
            if not os.path.exists(input_config[inp]):
                self.logger.error(f"Input file not found: {input_config[inp]}")
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

        contigs_fa = self.config['input']['contigs']
        hic_r1 = self.config['input']['hic_r1']
        hic_r2 = self.config['input']['hic_r2']
        threads = self.config['parameters'].get('threads', 8)

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

        # Align Hi-C reads
        self.logger.info("Aligning Hi-C reads...")
        hic_name_sorted = self.preprocessing_dir / "hic_name_sorted.bam"
        cmd = f"""bwa mem -t {threads} {contigs_fa} {hic_r1} {hic_r2} | \\
                 samtools view -bSh - | \\
                 samtools sort -n -O BAM -o {hic_name_sorted}"""
        self.run_command(cmd, "Align Hi-C reads with BWA and sort by name")

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

        contigs_fa = self.config['input']['contigs']
        output_prefix = self.scaffolding_dir / self.config['output'].get(
            'prefix', 'scaffolded'
            )

        # Build YaHS command
        cmd_parts = ["yahs", f"-o {output_prefix}", contigs_fa, str(self.hic_dedup_bam)]

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

        contact_map_method = self.config['parameters'].get('contact_map_method', 'both')

        if contact_map_method in ['juicer', 'both']:
            self.generate_juicer_hic()

        if contact_map_method in ['pretext', 'both']:
            self.generate_pretext_maps()

    def generate_juicer_hic(self):
        """Generate .hic file using Juicer Tools"""
        self.logger.info("Generating .hic file with Juicer Tools...")

        alignments_file = f"{self.output_prefix}_alignments_sorted.txt"
        hic_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.hic"
        chrom_sizes = f"{self.output_prefix}.fa.fai"

        # Check if alignments file exists (should be generated by YaHS)
        if not os.path.exists(alignments_file):
            self.logger.warning(f"Alignments file not found: {alignments_file}")
            self.logger.warning("Skipping Juicer .hic generation")
            return

        memory = self.config['parameters'].get('java_memory', '32G')
        juicer_jar = self.config['parameters'].get(
            'juicer_tools_jar', '/usr/local/bin/juicer_tools.jar'
            )

        cmd = f"""java -Xmx{memory} -jar {juicer_jar} pre \\
                 {alignments_file} {hic_file} {chrom_sizes}"""

        self.run_command(cmd, "Generate .hic file with Juicer Tools")

    def generate_pretext_maps(self):
        """Generate Pretext maps"""
        self.logger.info("Generating Pretext maps...")

        alignments_file = f"{self.output_prefix}_alignments_sorted.txt"
        pretext_file = self.contact_maps_dir / f"{self.config['output'].get('prefix', 'scaffolded')}.pretext"

        if not os.path.exists(alignments_file):
            self.logger.warning(f"Alignments file not found: {alignments_file}")
            self.logger.warning("Skipping Pretext map generation")
            return

        # Generate pretext map
        cmd = f"""PretextMap -o {pretext_file} \\
                 <(awk '{{print ".\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t.\\t."}}' {alignments_file})"""

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