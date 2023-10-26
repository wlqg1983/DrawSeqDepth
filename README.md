# DrawSeqDepth
The software tool, known as "DrawSeqDepth," appears to be a valuable resource for researchers working with organelle genomes, offering an easy-to-use solution for assessing the quality of genome assemblies. It could potentially streamline the research process and reduce the risk of errors in downstream analyses.

**1 Install the "DrawSeqDepth"**

1.1 Install and activate the conda

$ wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda 

$ echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc

1.2 Set up the conda working environment of DrawSeqDepth

Install the packaged DrawSeqDepth program into the conda working environment.

$ conda create --name DrawSeqDepth python=3.9 -y

$ conda activate DrawSeqDepth

$ wget -c https://github.com/wlqg1983/DrawSeqDepth/archive/refs/heads/main.zip

$ unzip main.zip 

$ chmod 551 DrawSeqDepth-main/* 

$ cp DrawSeqDepth-main/* $CONDA_PREFIX/bin 

$ rm -rf DrawSeqDepth-main main.zip

1.3 Install sra-tools 2.11.0 in offline mode, as the conda network installation of sra-tools does not meet the unzipping requirements for this study.

$ wget -c https://anaconda.org/bioconda/sra-tools/2.11.0/download/linux-64/sra-tools-2.11.0-pl5321ha49a11a_3.tar.bz2 

$ conda install --offline -f sra-tools-2.11.0-pl5321ha49a11a_3.tar.bz2

1.4 Install other dependent programs.

$ conda install -c bioconda samtools=1.17 -y

$ conda install -c bioconda minimap2=2.26 -y

$ conda install -c bioconda bwa=0.7.17 -y

$ conda install -c bioconda bowtie2=2.5.1 -y

$ conda install -c bioconda hisat2=2.2.1 -y

1.5 Install the dependent Python packages.

$ conda install -c conda-forge numpy -y

$ conda install -c conda-forge pandas -y

$ conda install -c conda-forge matplotlib -y

$ conda install -c conda-forge requests -y


**2 Download the sequencing data and reference organelle genomes**

Download and unzip the raw sequencing data:

$ prefetch SRR20647929 SRR12597239 SRR14924549

$ fastq-dump --split-files SRR20647929 SRR12597239 SRR14924549

Download the reference genomes:

$ prefetch_fasta.py ON055287 MT872375 MW553042

If users encounter the error message "This sra toolkit installation has not been configured," please execute the following command: 
$ vdb-config --interactive

The above script will take users to the interactive configuration interface. Enter 'X' to proceed and continue using sratools.


**3 Drawing the genome sequencing depth as a pipeline using default parameters**

The pseudo-code for running DrawSeqDepth is provided below:

$ DrawSeqDepth.py [-h | --help] [-line | -circle] [-alignment <minimap2 | hisat2 | bowtie2 | bwa>] [-reference REFERENCE] [-single FASTQ | -pair FASTQ1 FASTQ2 | -third FASTQ] [-format <png | pdf | jpeg | bmp | tiff | gif | eps | svg>] [-output OUTPUT] [-flag FLAG]

The arguments are described below in detail:

-h, --help: Show the help message of DrawSeqDepth and exit.

-alignment: Optional. Specifies the alignment program to be used, with optional values such as "minimap2", "hisat2", "bowtie2", and "bwa". The default alignment program is minimap2.

-reference: Specifies the reference genome file or a specific genome segment.

-format: Specifies the output image format, with optional values "png", "pdf", "jpeg", "bmp", "tiff", "gif", "eps", and "svg". The default format is "pdf".

-single: If the sequencing type is single-end, use this option to specify the filename of a single FASTQ file.

-pair: If the sequencing type is paired-end, use this option to specify the filenames of a couple of FASTQ files.

-third: If the sequencing type is third-generation, use this option to specify the filename of the third-generation sequencing data.

-output: Specifies the folder name of output files.

-flag: Optional. Specifies the value of the samtools flag parameter; the default is 4.

-line/-circle: Optional. Indicates whether the reference genome is circular (-circle) or linear (-line). If this parameter is missing, the default is linear.

The above parameters do not need to be in order.


**4 Three specific examples**

(1) Align the single-read next-generation sequencing data to a plastome based on minimap2:

$ DrawSeqDepth.py -reference ON055287.fasta -single SRR20647929_1.fastq -output ON055287.minimap2.single.depth

The script above is the most streamlined code for the DrawSeqDepth program. The results obtained are saved in the folder ON055287.minimap2.single.depth_results. The genome is represented as a linear sequence. The sequencing depth and coverage map is shown in Figure 1A. The format of the output figure is “pdf”. The flag parameter for samtools is 4, excluding unmatched reads. 

(2) Align the paired-end next-generation sequencing data to a plastome based on bowtie2:

$ DrawSeqDepth.py -line -alignment bowtie2 -reference MT872375.fasta -pair SRR12597239_1.fastq SRR12597239_2.fastq -format png -output MT872375.bowtie2.single.depth -flag 12

The script above is the most complete code available for the DrawSeqDepth program. The results obtained are saved in the folder MT872375.bowtie2.single.depth_results. The genome is considered to be a linear sequence. The flag value is set to 12 (=4+8), excluding unmatched reads (4) and secondary alignment reads (8). The sequencing depth and coverage map is shown in Figure 1B. 

(3) Align the third-generation sequencing data (Nanopore) to a mitogenome based on minimap2:

$ DrawSeqDepth.py -circle -reference MW553042.fasta -third SRR14924549_1.fastq -output MW553042.minimap2.third.depth -format jpeg

If the mitogenome is a circular sequence, two results for the whole genome and the junction sequence are saved in MW553042.minimap2.third.depth_results and MW553042.minimap2.third.depth_junction_results, respectively. The sequencing depth and coverage map of the whole genome and the junction sequence are shown in Figures 1C and 1D. 


