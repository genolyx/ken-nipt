# Use Ubuntu 22.04 LTS for security and latest packages
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

# Error 발생시 fix_gpg_keys.sh 로 해결가능
#RUN echo 'Acquire::AllowInsecureRepositories "true";' > /etc/apt/apt.conf.d/99insecure && \
#    echo 'Acquire::AllowDowngradeToInsecureRepositories "true";' >> /etc/apt/apt.conf.d/99insecure

# Install tzdata and set timezone
RUN apt-get update && \
    apt-get install -y tzdata && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Update package lists and install software-properties-common
RUN apt-get update && \
    apt-get install -y --no-install-recommends software-properties-common && \
    add-apt-repository universe && \
    apt-get update && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Install all required packages + ICU libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3 python3-pip \
    python2 python2-dev \
    curl \
    openjdk-11-jdk \
    git wget unzip \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libncurses-dev libreadline-dev \
    gfortran gcc make \
    pkg-config \
    libicu-dev \
    libicu70 \
    bedtools \
    bwa \
    samtools \
    default-jre \
    fastqc \
    seqtk \
    vim \
    zip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Python 2용 pip 설치
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py | python2.7

# Python 2 의존성 패키지 설치 (NumPy 1.14.6으로 고정)
RUN python2.7 -m pip install numpy==1.14.6 matplotlib pysam scipy futures subprocess32 scikit-learn==0.20.4

# Python 3 기본 도구만 설치 (WisecondorX는 나중에 conda 환경에서)
RUN python3 -m pip install --upgrade pip

# Set up work directories
RUN mkdir -p /Work/NIPT/bin \
    /Work/NIPT/data \
    /Work/NIPT/data/bed \
    /Work/NIPT/analysis \
    /Work/NIPT/log \
    /Work/NIPT/fastq

# ----------------------------
# Install external tools
# ----------------------------
# Install bwa-mem2
RUN cd /tmp && \
    wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xvjf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/ && \
    rm -rf /tmp/*

# Install Qualimap
#COPY --from=qualimap /opt/qualimap /opt/qualimap
#RUN ln -s /opt/qualimap/qualimap /usr/local/bin/qualimap
RUN apt-get update && apt-get install -y wget && \
    wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash miniconda.sh -b -p /opt/miniconda && \
    rm miniconda.sh && \
    /opt/miniconda/bin/conda config --set always_yes yes && \
    /opt/miniconda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/miniconda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    /opt/miniconda/bin/conda install -c bioconda qualimap -y && \
    ln -s /opt/miniconda/bin/qualimap /usr/local/bin/qualimap

# Install Picard with wrapper script
RUN mkdir -p /Work/NIPT/bin/picard && \
    wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -O /Work/NIPT/bin/picard/picard.jar && \
    echo '#!/bin/bash\njava "$@" -jar /Work/NIPT/bin/picard/picard.jar' > /usr/local/bin/picard && \
    chmod +x /usr/local/bin/picard

# ----------------------------
# Install Conda and create environment
# ----------------------------
ENV PATH=/opt/conda/bin:$PATH
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    /opt/conda/bin/conda clean -afy

# ----------------------------
# Create conda environment with Python 3.12.2
# ----------------------------
RUN /opt/conda/bin/conda create -n nipt -y python=3.12.2

# ----------------------------
# Install bioinformatics tools via Conda (WisecondorX 제외)
# ----------------------------
RUN /opt/conda/bin/conda run -n nipt \
    /opt/conda/bin/conda install -y -c conda-forge -c bioconda -c defaults \
    hmmcopy \
    r-base=4.4 \
    r-essentials

# ----------------------------
# Install R packages via Conda
# ----------------------------
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate nipt && \
    conda install -y -c conda-forge mamba && \
    mamba install -y -c conda-forge -c bioconda -c r \
    bioconductor-rsamtools \
    bioconductor-genomicalignments \
    bioconductor-genomicranges \
    bioconductor-rtracklayer \
    bioconductor-s4vectors \
    bioconductor-iranges \
    bioconductor-hmmcopy \
    bioconductor-dnacopy \
    r-argparse \
    r-optparse \
    r-r6 \
    r-jsonlite \
    r-data.table \
    r-ggplot2 \
    r-plyr \
    r-mass \
    r-e1071"

# ----------------------------
# Install stringr/stringi (conda 환경에서)
# ----------------------------
RUN /opt/conda/bin/conda run -n nipt \
    /opt/conda/bin/conda install -y -c conda-forge \
    r-stringi r-stringr

# R stringr 검증
RUN /opt/conda/envs/nipt/bin/R --slave -e "\
    library(stringr); \
    result <- str_replace('test', 'test', 'success'); \
    if (result != 'success') stop('StringR test failed'); \
    message('✓ StringR verified successfully');"

# ----------------------------
# seqFF 에러 방지
# ----------------------------
RUN /opt/conda/envs/nipt/bin/R --slave -e "\
    install.packages(c('tibble', 'dplyr'), repos='https://cloud.r-project.org', dependencies=TRUE, force=TRUE)"

# ----------------------------
# Install WisecondorX official recommended dependencies in conda environment
# ----------------------------
RUN /opt/conda/envs/nipt/bin/pip install \
    numpy==1.26.4 \
    scipy==1.13.0 \
    scikit-learn==1.4.2 \
    plotly \
    pandas==2.2.2 \
    matplotlib==3.8.4 \
    pysam==0.22.0

# ----------------------------
# Install WisecondorX in conda environment
# ----------------------------
RUN /opt/conda/envs/nipt/bin/pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX

# Verify WisecondorX installation
RUN /opt/conda/envs/nipt/bin/wisecondorx --help || echo "WisecondorX installed in conda environment"

# Create symbolic link only if not exists
RUN if [ ! -f /usr/local/bin/wisecondorx ]; then \
        ln -s /opt/conda/envs/nipt/bin/wisecondorx /usr/local/bin/wisecondorx; \
    fi && \
    if [ ! -f /usr/local/bin/readCounter ]; then \
        ln -s /opt/conda/envs/nipt/bin/readCounter /usr/local/bin/readCounter; \
    fi

# ----------------------------
# Install Wisecondor (Python 2) - Original version
# ----------------------------
RUN git clone https://github.com/VUmcCGP/wisecondor.git /opt/wisecondor && \
    sed -i '1i #!/usr/bin/env python2' /opt/wisecondor/wisecondor.py && \
    chmod +x /opt/wisecondor/wisecondor.py && \
    ln -s /opt/wisecondor/wisecondor.py /usr/local/bin/wisecondor

# ----------------------------
# Install Wisecondor-FF (from local files)
# ----------------------------
COPY bin/scripts/WisecondorFF /opt/WisecondorFF
RUN chmod +x /opt/WisecondorFF/wisecondorff && \
    sed -i 's|/opt/WisecondorFF/main.py|/opt/WisecondorFF/src/main.py|g' /opt/WisecondorFF/wisecondorff && \
    ln -s /opt/WisecondorFF/wisecondorff /usr/local/bin/wisecondor-ff

# Copy HMMcopy reference WIG files
RUN mkdir -p /opt/conda/envs/nipt/share/hmmcopy
COPY data/wig/hg19.gc.50kb.wig /opt/conda/envs/nipt/share/hmmcopy/
COPY data/wig/hg19.map.50kb.wig /opt/conda/envs/nipt/share/hmmcopy/
COPY data/wig/hg19.gc.10mb.wig /opt/conda/envs/nipt/share/hmmcopy/
COPY data/wig/hg19.map.10mb.wig /opt/conda/envs/nipt/share/hmmcopy/

# Copy BED files for filtering
# In Docker run, there's /Work/NIPT/data. So, it doesn't need to be copied
#COPY resources/bed/ /Work/NIPT/data/bed/

# Copy AWK scripts and R scripts
COPY bin/scripts/fetus.awk /Work/NIPT/bin/
COPY bin/scripts/mom.awk /Work/NIPT/bin/
COPY bin/scripts/HMMcopy.R /Work/NIPT/bin/

# ----------------------------
# SeqFF setup
# ----------------------------
RUN mkdir -p /opt/nipt/bin/scripts/seqFF_R
COPY bin/scripts/seqFF_R /opt/nipt/bin/scripts/seqFF_R

# SeqFF path fixes
RUN if [ -f /opt/nipt/bin/scripts/seqFF_R/seqff.r ]; then \
        sed -i 's|load("pd4615-sup-0008-file1.rdata")|load("/opt/nipt/bin/scripts/seqFF_R/pd4615-sup-0008-file1.rdata")|g' /opt/nipt/bin/scripts/seqFF_R/seqff.r; \
    fi

# SeqFF wrapper script
RUN printf '#!/bin/bash\n\
export PATH=/opt/conda/envs/nipt/bin:$PATH\n\
cd /opt/nipt/bin/scripts/seqFF_R\n\
/opt/conda/envs/nipt/bin/Rscript seqff.r "$@"\n' > /usr/local/bin/seqff && \
    chmod +x /usr/local/bin/seqff

ARG CACHE_BREAK=1

# Copy modules
COPY bin/scripts/modules/prizm_runner.py /Work/NIPT/bin/
COPY bin/scripts/modules/ezd_runner.py /Work/NIPT/bin/
COPY bin/scripts/modules/sca_detector.py /Work/NIPT/bin/
COPY bin/scripts/modules/process_md_result.py /Work/NIPT/bin/
COPY bin/scripts/modules/generate_json_output.py /Work/NIPT/bin/
COPY bin/scripts/modules/html_review_page.py /Work/NIPT/bin/
COPY bin/scripts/modules/gender_detector.py /Work/NIPT/bin/

# Copy and setup pipeline
COPY bin/scripts/nipt_pipeline.py /Work/NIPT/bin/
RUN chmod +x /Work/NIPT/bin/nipt_pipeline.py

# Generate Reference for WC family
COPY bin/scripts/generate_wcfamily_reference.sh /Work/NIPT/bin/
RUN chmod +x /Work/NIPT/bin/generate_wcfamily_reference.sh

# Set environment variables
ENV SRC_DIR=/Work/NIPT/bin \
    ANALYSIS_DIR=/Work/NIPT/analysis \
    LOG_DIR=/Work/NIPT/log \
    FASTQ_DIR=/Work/NIPT/fastq \
    BWA2=bwa-mem2 \
    SAMTools=samtools \
    PICARD=picard \
    qualimap=qualimap \
    WC=/opt/wisecondor \
    WCX=wisecondorx \
    WCFF=wisecondor-ff \
    HMMcopy=/opt/conda/envs/nipt \
    PYTHON2=python2.7 \
    Rscript=/opt/conda/envs/nipt/bin/Rscript \
    R_HOME=/opt/conda/envs/nipt/lib/R \
    R_LIBS_USER=/opt/conda/envs/nipt/lib/R/library \
    PATH=/Work/NIPT/bin:/opt/conda/envs/nipt/bin:$PATH

# Activate conda environment by default
RUN echo "source /opt/conda/bin/activate nipt" >> ~/.bashrc

WORKDIR /Work/NIPT

# Entry point
ENTRYPOINT ["python3", "/Work/NIPT/bin/nipt_pipeline.py"]
CMD ["--help"]
