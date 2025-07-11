Bootstrap: docker
From: ubuntu:24.04

%labels
    Author Petr Novak
    Version 0.1
    Description "Singularity container for YaHS scaffolding pipeline with Hi-C tools"

%files
    yahs_workflow.py /opt/yahs-workflow/yahs_workflow.py
    config_template.yaml /opt/yahs-workflow/config_template.yaml

%post
    # Set timezone to avoid interactive prompts
    export DEBIAN_FRONTEND=noninteractive
    ln -fs /usr/share/zoneinfo/Europe/Prague /etc/localtime
    echo 'Europe/Prague' > /etc/timezone

    # Configure alternative Ubuntu mirrors for 24.04 (noble)
    cat > /etc/apt/sources.list << 'EOF'
# Main Ubuntu repositories
deb http://archive.ubuntu.com/ubuntu/ noble main restricted universe multiverse
deb http://archive.ubuntu.com/ubuntu/ noble-updates main restricted universe multiverse
deb http://archive.ubuntu.com/ubuntu/ noble-backports main restricted universe multiverse
deb http://security.ubuntu.com/ubuntu/ noble-security main restricted universe multiverse

# Czech Republic mirrors (closer to Prague)
# deb http://cz.archive.ubuntu.com/ubuntu/ noble main restricted universe multiverse
# deb http://cz.archive.ubuntu.com/ubuntu/ noble-updates main restricted universe multiverse
# deb http://mirror.hosting90.cz/ubuntu/ noble main restricted universe multiverse
EOF

    # Install core build tools and libraries
    apt-get update --fix-missing && apt-get install -y \
        build-essential git cmake zlib1g-dev wget default-jdk \
        bwa samtools bedtools biobambam2 curl unzip \
        python3 python3-yaml tzdata

    # Install YaHS
    cd /opt && git clone https://github.com/c-zhou/yahs.git
    cd yahs && make && cp yahs /usr/local/bin/

    # Download Juicer Tools
    wget -O /usr/local/bin/juicer_tools.jar \
        https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
    
    # Create juicer wrapper script for 'juicer pre' command
    cat > /usr/local/bin/juicer << 'EOF'
#!/bin/bash
# Juicer wrapper script
if [ "$1" = "pre" ]; then
    shift
    java -jar /usr/local/bin/juicer_tools.jar pre "$@"
else
    echo "Usage: juicer pre [options]"
    echo "This is a wrapper for Juicer Tools pre command"
    exit 1
fi
EOF
    chmod +x /usr/local/bin/juicer

    # Install Miniconda
    wget -O /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash /tmp/miniconda.sh -b -p /opt/miniconda
    rm /tmp/miniconda.sh

    # Add conda to PATH and install pretext-suite from bioconda
    export PATH="/opt/miniconda/bin:$PATH"
    /opt/miniconda/bin/conda config --add channels bioconda
    /opt/miniconda/bin/conda config --add channels conda-forge
    /opt/miniconda/bin/conda install -y pretext-suite pyyaml sambamba

    # Copy workflow script to container
    mkdir -p /opt/yahs-workflow
    chmod +x /opt/yahs-workflow/yahs_workflow.py

%environment
    export PATH=/opt/miniconda/bin:/usr/local/bin:/opt/yahs-workflow:$PATH
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%runscript
    if [ $# -eq 0 ]; then
        echo "YaHS Scaffolding Workflow Container"
        echo "Usage: singularity run container.sif config.yaml"
        echo "       singularity exec container.sif yahs_workflow.py config.yaml"
        echo ""
        echo "Available tools: yahs, bwa, samtools, bedtools, sambamba, PretextMap, PretextSnapshot"
        echo "Template config: /opt/yahs-workflow/config_template.yaml"
        exit 1
    fi
    python3 /opt/yahs-workflow/yahs_workflow.py "$@"

%test
    echo "Checking installed tools versions..."
    echo -n "bwa: " && bwa 2>&1 | head -n1
    echo -n "samtools: " && samtools --version | head -n1
    echo -n "bedtools: " && bedtools --version | head -n1
    echo -n "biobambam2 bammarkduplicates2: " && bammarkduplicates2 --version || echo "OK"
    echo -n "sambamba: " && sambamba --version 2>&1 | head -n1
    echo -n "YaHS: " && yahs --version || echo "YaHS installed"
    echo -n "Java: " && java -version 2>&1 | head -n1
    echo -n "JuicerTools: " && java -jar /usr/local/bin/juicer_tools.jar 2>&1 | head -n1
    echo -n "Juicer pre: " && juicer pre 2>&1 | head -n1
    echo -n "PretextMap: " && PretextMap --help | head -n1
    echo -n "PretextSnapshot: " && PretextSnapshot --help | head -n1
    echo -n "Conda: " && conda --version
    python3 - << 'EOF'
import yaml
print("pyyaml OK")
EOF
